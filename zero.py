import iio
import cffi
import numpy as np
import matplotlib.pyplot as plt
# from skimage.util import view_as_blocks
import scipy as sp
# import sys
# import scipy.fft
import os

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
int detect_forgeries(int * votes, int * forgery_mask, int * forgery_mask_reg,
                     meaningful_reg * forged_regions,
                     int X, int Y, int grid_to_exclude, int grid_max);
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
    elif array.dtype == bool:
        typestr = 'bool*'
    elif array.dtype == np.int32:
        typestr = 'int*'
    # requires cffi 0.12
    return ffi.from_buffer(typestr, array, require_writable=True)

def main(filename):
    image = iio.read(filename).astype(np.float64)
    h, w, c = image.shape
    image = image.transpose((2, 0, 1))
    image = image.copy(order='C')

    print('1. convert to luminance\n')
    im = np.zeros((h,w), dtype=np.float64)
    im = im.copy(order='C')

    libzero.rgb2luminance(P(image), P(im), w, h, c)
    iio.write('luminance.png', im)

    # intermediate step before statistical validation
    print('2. compute vote map\n')
    votes = np.zeros(im.shape, dtype=np.int32)
    libzero.compute_grid_votes_per_pixel(P(im), P(votes), w, h)

    print('2bis. color vote map\n')
    colored_votes = 255 * colormap(votes)[...,:3]
    colored_votes[votes == -1] = 0
    colored_votes = colored_votes.astype(np.uint8)

    # one color per grid origin + black in case of a tie
    iio.write('colored_votemap.png', colored_votes)

    print('3. detect global grids\n')
    lnfa_grids = np.zeros((8, 8), dtype=np.float64)
    main_grid = libzero.detect_global_grids(P(votes), P(lnfa_grids), w, h)
    significant_grids = np.where(lnfa_grids < 0.0)

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
    forgery_result = np.zeros(im.shape, dtype=np.int32)

    forged_region = ffi.new('meaningful_reg[]', w*h)
    forgery_found = libzero.detect_forgeries(P(votes), P(forgery), P(forgery_c),
                                             forged_region, w, h, main_grid, 63)

    if forgery_found > 0:
        for i in range(forgery_found):
            print("foreign grid was found here: " + str(forged_region[i].x0) + " "
                  + str(forged_region[i].y0) + " - " + str(forged_region[i].x1) + " "
                  + str(forged_region[i].y1))
            print("grid is " + str(forged_region[i].grid%8) + ","
                  + str(int(forged_region[i].grid/8))
                  + " with log(nfa) = " + str(forged_region[i].lnfa))

    forgery_result = forgery_c

    # do the rest only if main grid is detected
    if main_grid > -1:
        print('\n5. create simulated version\n')

        # create JPEG file with PIL
        ######################################################
        from PIL import Image
        pil_image = Image.open(filename)
        pil_image.save('version99.jpg', format='JPEG', quality=99)
        pil_image = iio.read('version99.jpg').astype(np.float64)

        ######################################################

        h, w, c  = pil_image.shape
        pil_image  = pil_image.transpose((2, 0, 1))
        pil_image = pil_image.copy(order='C')

        img = np.zeros((h,w), dtype=np.float64)
        img = img.copy(order='C')

        libzero.rgb2luminance(P(pil_image), P(img), w, h, c)

        votes2 = np.zeros(img.shape, dtype=np.int32)
        libzero.compute_grid_votes_per_pixel(P(img), P(votes2), w, h)

        iio.write('votes2.png', votes2)

        nb_globalgrids = len(np.where(lnfa_grids<0)[0])
        for i in range(nb_globalgrids):
            x = np.where(lnfa_grids<0)[0][i]
            y = np.where(lnfa_grids<0)[1][i]
            coordgrid = x*8+y
            votes2[votes == coordgrid] = -1

        print('2bis. color vote map\n')
        colored_votes = 255 * colormap(votes2)[...,:3]
        colored_votes[votes2 == -1] = 0
        colored_votes = colored_votes.astype(np.uint8)

        # one color per grid origin + black in case of a tie
        iio.write('colored_votemap_new.png', colored_votes)


        print('\n5. detect suspicious areas\n')
        forgery = np.zeros(img.shape, dtype=np.int32)
        forgery_c2 = np.zeros(img.shape, dtype=np.int32)
        forged_region = ffi.new('meaningful_reg[]', w*h)

        forgery_found2 = libzero.detect_forgeries(P(votes2), P(forgery), P(forgery_c2),
                                                  forged_region, w, h, -1, 0)

        if forgery_found2 > 0:
            for i in range(forgery_found2):
                print("an absence of grid was found here: " + str(forged_region[i].x0) + " "
                      + str(forged_region[i].y0) + " - " + str(forged_region[i].x1) + " "
                      + str(forged_region[i].y1))
                print("with log(nfa) = " + str(forged_region[i].lnfa))
                forgery = forgery_c  + 0.5*forgery_c2
                forgery_result = np.clip(forgery, 0, 255)

    iio.write('result_zero.png', forgery_result) # all black if no forgeries

    print('\nok')


if __name__ == '__main__':
    import fire
    fire.Fire(main)
