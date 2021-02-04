import iio
import cffi
import numpy as np

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
        im = image[:,:,0] * 299/1000 + image[:,:,1] * 587/1000 + image[:,:,2] * 114/1000
    else:
        im = image[:,:,0]

    # libzero.rgb2luminance(P(image), P(im), w, h, c)

    iio.write('luminance.png', im)

    print('2. compute vote map\n')
    votes = np.zeros(im.shape, dtype=np.int32)
    libzero.compute_grid_votes_per_pixel(P(im), P(votes), w, h)

    iio.write('votemap.png', votes)

    print('3. detect global grids\n')
    lnfa_grids = np.zeros((8, 8), dtype=np.float64)
    main_grid = libzero.detect_global_grids(P(votes), P(lnfa_grids), w, h)

    if main_grid == -1:
        print('No overall JPEG grid found')
    else:
        print("main grid is " + str(main_grid%8) + "," + str(int(main_grid/8)) )

    if main_grid > 0:
        print('The most meaningful JPEG grid origin is not (0,0).\n'
               'This may indicate that the image has been cropped.\n')

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
            print("\ngrid is " + str(forged_region[i].grid%8) + "," + str(int(forged_region[i].grid/8))
                  + " with log(nfa) = " + str(forged_region[i].lnfa))
        iio.write('forgery.png', forgery_c)

    print('\nok')

if __name__ == '__main__':
    import fire
    fire.Fire(main)
