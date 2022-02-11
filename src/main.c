/*----------------------------------------------------------------------------

  ZERO: JPEG grid detector applied to forgery detection in digital images. This
  code is part of the following publication and was subject to peer review:

    "ZERO: a local JPEG grid origin detector based on the number of DCT zeros
    and its applications in image forensics" by Tina Nikoukhah, Jérémy Anger,
    Miguel Colom, Jean-Michel Morel, Rafael Grompone von Gioi,
    Image Processing On Line, 2021.
    http://dx.doi.org/10.5201/ipol.2021.XXX

  Copyright (c) 2018-2021 Rafael Grompone von Gioi <grompone@gmail.com>
  Copyright (c) 2018-2021 Jérémy Anger <jeremy.anger@ens-paris-saclay.fr>
  Copyright (c) 2018-2021 Tina Nikoukhah <tinanikoukhah@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <https://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "iio.h"
#include "zero.h"

/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
    double * input = NULL;
    double * input_jpeg = NULL;
    double * luminance; double * luminance_jpeg;
    int * votes; int * votes_jpeg;
    int X, Y, C, XX, YY, CC;
    double lnfa_grids[64] = {0.0};
    meaningful_reg * foreign_regions;
    int foreign_regions_n = 0;
    meaningful_reg * missing_regions;
    int missing_regions_n = 0;
    int * mask_f; int * mask_f_reg;
    int * mask_m; int * mask_m_reg;
    int main_grid = -1;
    int global_grids = 0;

    /* check input and usage */
    if (argc != 2 && argc != 3)
        error("usage: zero <image> [image_jpeg99]\n"
              "finds JPEG grids and forgeries");

    /* read input */
    input = iio_read_image_double_split(argv[1], &X, &Y, &C);
    if (argc == 3) {
        input_jpeg = iio_read_image_double_split(argv[2], &XX, &YY, &CC);
        if (X != XX || Y != YY )
            error("image and image_jpeg99 have different size");
    }

    /* allocate memory */
    luminance       = (double *) xcalloc(X*Y, sizeof(double));
    luminance_jpeg  = (double *) xcalloc(X*Y, sizeof(double));
    votes           = (int *) xcalloc(X * Y, sizeof(int));
    votes_jpeg      = (int *) xcalloc(X * Y, sizeof(int));
    foreign_regions = (meaningful_reg *) xcalloc(X*Y, sizeof(meaningful_reg));
    missing_regions = (meaningful_reg *) xcalloc(X*Y, sizeof(meaningful_reg));
    mask_f          = (int *) xcalloc(X * Y, sizeof(int));
    mask_f_reg      = (int *) xcalloc(X * Y, sizeof(int));
    mask_m          = (int *) xcalloc(X * Y, sizeof(int));
    mask_m_reg      = (int *) xcalloc(X * Y, sizeof(int));

    /* run algorithm */
    main_grid = zero(input, input_jpeg, luminance, luminance_jpeg,
                     votes, votes_jpeg, lnfa_grids,
                     foreign_regions, &foreign_regions_n,
                     missing_regions, &missing_regions_n,
                     mask_f, mask_f_reg, mask_m, mask_m_reg,  X, Y, C, CC);

    /* print detection result */

    if (main_grid == -1)
        /* main grid not found */
        printf("No overall JPEG grid found.\n");

    if (main_grid > -1) {
        /* print main grid */
        printf("main grid found: #%d (%d,%d) log(nfa) = %g\n", main_grid,
               main_grid % 8, main_grid / 8, lnfa_grids[main_grid]);
        global_grids++;
    }

    for (int i=0; i<64; i++) {
        /* print list of meaningful grids */
        if (lnfa_grids[i] < 0.0 && i != main_grid) {
            printf("meaningful global grid found: #%d (%d,%d) log(nfa) = %g\n",
                   i, i % 8, i / 8, lnfa_grids[i]);
            global_grids++;
        }
    }

    if (foreign_regions_n != 0) {
        for (int i=0; i<foreign_regions_n; i++) {
            if (main_grid != -1)
                printf("\nA meaningful grid different from the main one "
                       "was found here:\n");
            else
                printf("\nA meaningful grid was found here:\n");
            printf("bounding box: %d %d to %d %d [%dx%d]",
                   foreign_regions[i].x0, foreign_regions[i].y0,
                   foreign_regions[i].x1, foreign_regions[i].y1,
                   foreign_regions[i].x1 - foreign_regions[i].x0+1,
                   foreign_regions[i].y1 - foreign_regions[i].y0+1);
            printf(" grid: #%d (%d,%d)", foreign_regions[i].grid,
                   foreign_regions[i].grid % 8, foreign_regions[i].grid / 8 );
            printf(" log(nfa) = %g\n", foreign_regions[i].lnfa);
        }
    }

    if (main_grid > -1 && missing_regions_n > 0) {
        for (int i=0; i<missing_regions_n; i++) {
            printf("\nA region with missing JPEG grid was found here:\n");
            printf("bounding box: %d %d to %d %d [%dx%d]",
                   missing_regions[i].x0, missing_regions[i].y0,
                   missing_regions[i].x1, missing_regions[i].y1,
                   missing_regions[i].x1 - missing_regions[i].x0+1,
                   missing_regions[i].y1 - missing_regions[i].y0+1);
            printf(" grid: #%d (%d,%d)", missing_regions[i].grid,
                   missing_regions[i].grid % 8, missing_regions[i].grid / 8 );
            printf(" log(nfa) = %g\n", missing_regions[i].lnfa);
        }
    }

    if (foreign_regions_n + missing_regions_n == 0 && main_grid < 1)
        printf("\nNo suspicious traces found in the image "
               "with the performed analysis.\n");

    if (main_grid > 0)
        printf("\nThe most meaningful JPEG grid origin is not (0,0).\n"
               "This may indicate that the image has been cropped.\n");

    if (global_grids > 1)
        printf("\nThere is more than one meaningful grid. "
               "This is suspicious.\n");

    if (foreign_regions_n + missing_regions_n > 0) {
        printf("\nSuspicious traces found in the image.\nThis may be caused "
               "by image manipulations such as resampling, \ncopy-paste, "
               "splicing.  Please examine the deviant meaningful region \n"
               "to make your own opinion about a potential forgery.\n");
    }

    /* store vote map and forgery detection outputs */
    iio_write_image_double("luminance.png", luminance, X, Y);
    iio_write_image_int("votes.png", votes, X, Y);
    iio_write_image_int("votes_jpeg.png", votes_jpeg, X, Y);
    iio_write_image_int("mask_f.png", mask_f_reg, X, Y);
    iio_write_image_int("mask_m.png", mask_m_reg, X, Y);

    /* free memory */
    free((void *) input);
    if (input_jpeg != NULL) free((void *) input_jpeg);
    free((void *) luminance);
    free((void *) luminance_jpeg);
    free((void *) votes);
    free((void *) votes_jpeg);
    free((void *) foreign_regions);
    free((void *) missing_regions);
    free((void *) mask_f);
    free((void *) mask_f_reg);
    free((void *) mask_m);
    free((void *) mask_m_reg);

    return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
