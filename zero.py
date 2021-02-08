import iio
import cffi
import numpy as np
import matplotlib.pyplot as plt


ffi = cffi.FFI()
ffi.cdef('''
void rgb2luminance(double * input, double * output, int X, int Y, int C);
void compute_grid_votes_per_pixel(double * image, int * votes, int X, int Y);
int detect_global_grids(int * votes, double * lnfa_grids, int X, int Y);
typedef struct {
    int x0, y0, x1, y1;
    int grid;
    double lnfa;
} meaningful_reg;
int detect_forgery(int * votes, int * forgery, int * forgery_e, meaningful_reg * forged_regions,
                   int X, int Y, int main_grid);
''')

libzero = ffi.dlopen('./libzero.so')

cmap1 = plt.get_cmap('tab20')
cmap2 = plt.get_cmap('tab20b')
cmap3 = plt.get_cmap('tab20c')
cmap4 = plt.get_cmap('Set3')

def colormap(v):
    v = v.copy()
    # swap 0 and 4 so that 0 is colored with green
    v0 = v == 0
    v4 = v == 4
    v[v0] = 4
    v[v4] = 0
    v2 = cmap1(v/20)
    v2[v >= 20] = cmap2((v[v >= 20] - 20)/20)
    v2[v >= 40] = cmap3((v[v >= 40] - 40)/20)
    v2[v >= 60] = cmap4((v[v >= 60] - 60)/20)
    return v2

def P(array):
    typestr = 'double*'
    if array.dtype == np.float32:
        typestr = 'float*'
    elif array.dtype == np.bool:
        typestr = 'bool*'
    elif array.dtype == np.int32:
        typestr = 'int*'
    # requires cffi 0.12
    return ffi.from_buffer(typestr, array, require_writable=True)

def main(filename):
    image = iio.read(filename).astype(np.float64)
    h, w, c = image.shape

    print('1. convert to luminance\n')
    if c == 3:
        im = image[:,:,0] * 299/1000 + image[:,:,1] * 587/1000
        + image[:,:,2] * 114/1000
    else:
        im = image[:,:,0]

    iio.write('luminance.png', im)

    print('2. compute vote map\n')
    votes = np.zeros(im.shape, dtype=np.int32)
    libzero.compute_grid_votes_per_pixel(P(im), P(votes), w, h)

    # iio.write('votemap.png', votes)

    print('2bis. color vote map\n')

    colored_votes = 255 * colormap(votes)[...,:3]
    colored_votes[votes == -1] = 0

    colored_votes = colored_votes.astype(np.uint8)

    iio.write('colored_votemap.png', colored_votes)

    print('3. detect global grids\n')
    lnfa_grids = np.zeros((8, 8), dtype=np.float64)
    main_grid = libzero.detect_global_grids(P(votes), P(lnfa_grids), w, h)

    if main_grid == -1:
        print('No overall JPEG grid found') # this means the image has no detectable JPEG traces
    else:
        print("main grid is " + str(main_grid%8) + "," + str(int(main_grid/8)) )

    if main_grid > 0:
        print('The most meaningful JPEG grid origin is not (0,0).\n'
               'This may indicate that the image has been cropped.\n') # this means that the grid is not aligned

    for i in range(64):
        if lnfa_grids[int(i/8)][i%8] < 0.0:
            print("significant grid is "+ str(i%8) + "," + str(int(i/8))
                  + " with log(nfa) = " + str(lnfa_grids[int(i/8)][i%8]))


    print('\n4. detect forgeries\n')
    forgery = np.zeros(im.shape, dtype=np.int32)
    forgery_c = np.zeros(im.shape, dtype=np.int32)

    forged_region = ffi.new('meaningful_reg[]', w*h)

    forgery_found = libzero.detect_forgery(P(votes), P(forgery), P(forgery_c),
                                           forged_region, w, h, main_grid)

    if forgery_found == 0 and main_grid < 1:
        print('No suspicious traces found in the image with the performed analysis\n')

    if forgery_found > 0:
        for i in range(forgery_found):
            print("grid was found here: " + str(forged_region[i].x0) + " "
                  + str(forged_region[i].y0) + " - " + str(forged_region[i].x1) + " "
                  + str(forged_region[i].y1))
            print("\ngrid is " + str(forged_region[i].grid%8) + ","
                  + str(int(forged_region[i].grid/8))
                  + " with log(nfa) = " + str(forged_region[i].lnfa))

    iio.write('forgery.png', forgery_c) # all black if no forgeries

    print('\nok')

if __name__ == '__main__':
    import fire
    fire.Fire(main)
