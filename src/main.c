/*----------------------------------------------------------------------------

  Copyright (c) 2018-2020 Rafael Grompone von Gioi <grompone@gmail.com>
  Copyright (c) 2018-2020 Jérémy Anger <anger@cmla.ens-cachan.fr>
  Copyright (c) 2018-2020 Tina Nikoukhah <tinanikoukhah@gmail.com>

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

int main(int argc, char ** argv) {
    double * input;
    double * image;
    int * votes;
    int X,Y,C;
    double lnfa_grids[64] = {0.0};
    int main_grid = -1;
    int global_grids = 0;
    int forgery_found = 0;
    meaningful_reg * forged_regions;
    int * forgery;
    int * forgery_e;


    if (argc != 2) error("use: zero <image>\nfinds JPEG grids and forgeries");

    input = iio_read_image_double_split(argv[1], &X, &Y, &C);

    /* luminance image */
    image = (double *) xcalloc(X*Y, sizeof(double));

    rgb2luminance(input, image, X, Y, C);
    iio_write_image_double("luminance.png", image, X, Y);

    /* compute vote map */
    votes = (int *) xcalloc(X * Y, sizeof(int));

    compute_grid_votes_per_pixel(image, votes, X, Y);
    iio_write_image_int("votes.png", votes, X, Y);

    /* detect global grids and main_grid */
    main_grid = detect_global_grids(votes, lnfa_grids, X, Y);

    if (main_grid > -1) {
        /* print main grid */
        printf("main grid: #%d [%d %d] log(nfa) = %g\n", main_grid,
               main_grid % 8, main_grid / 8, lnfa_grids[main_grid]);
        global_grids++;
    }
    if (main_grid == -1)
        /* main grid not found */
        printf("No overall JPEG grid found.\n");

    if (main_grid > 0)
        printf("The most meaningful JPEG grid origin is not (0,0).\n"
               "This may indicate that the image has been cropped.\n");

    for (int i=0; i<64; i++) {
        /* print list of meaningful grids */
        if (lnfa_grids[i] < 0.0 && i != main_grid) {
            printf("significant global grid: #%d [%d %d] log(nfa) = %g\n", i,
                   i % 8, i / 8, lnfa_grids[i]);
            global_grids++;
        }
    }

    if (global_grids > 1)
        printf("There is more than one meaningful grid.\n"
               "This is suspicious.\n");

    /* compute forged regions, TODO: how to fill it?  */
    forged_regions = (meaningful_reg *) xcalloc(X*Y, sizeof(meaningful_reg));

    /* compute forgery masks */
    forgery = (int *) xcalloc(X * Y, sizeof(int));
    forgery_e = (int *) xcalloc(X * Y, sizeof(int));

    /* look for local significant grids different from the main one */
    forgery_found = detect_forgery(votes, forgery, forgery_e, forged_regions,
                                   X, Y, main_grid);

    if (forgery_found == 0 && main_grid < 1)
        printf("\nNo suspicious traces found in the image "
               "with the performed analysis.\n");

    if (forgery_found != 0) {
        for (int i=0; i<forgery_found; i++) {
            if (main_grid != -1)
                printf("\nA meaningful grid different from the main one was found here: ");
            else
                printf("\nA grid was found here: ");
            printf("%d %d - %d %d [%dx%d]", forged_regions[i].x0, forged_regions[i].y0,
                   forged_regions[i].x1, forged_regions[i].y1,
                   forged_regions[i].x1-forged_regions[i].x0+1,
                   forged_regions[i].y1-forged_regions[i].y0+1);
            printf("\ngrid: #%d [%d %d] ", forged_regions[i].grid,
                   forged_regions[i].grid % 8, forged_regions[i].grid / 8 );
            printf("log(nfa) = %g\n", forged_regions[i].lnfa);
        }

        /* store forgery detection outputs */
        iio_write_image_int("forgery.png", forgery, X, Y);
        iio_write_image_int("forgery_c.png", forgery_e, X, Y);

        printf("\nSuspicious traces found in the image.\n"
               "This may be caused by image manipulations such as resampling, \n"
               "copy-paste, splicing. Please examine the deviant meaningful region \n"
               "to make your own opinion about a potential forgery.\n");
    }

    /* free memory */
    free((void *) input);
    free((void *) image);
    free((void *) votes);
    free((void *) forged_regions);
    free((void *) forgery);
    free((void *) forgery_e);



    return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
