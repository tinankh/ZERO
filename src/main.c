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
    double * input;  //double * input2;
    double * image;  //double * image2;
    int X, Y, C;  //int X2, Y2, C2;
    int * votes;  //int * votes2;
    double lnfa_grids[64] = {0.0};
    meaningful_reg * forged_regions;  //meaningful_reg * forged_regions2;
    int * forgery;  //int * forgery2;
    int * forgery_e;  //int * forgery_e2;
    int main_grid = -1;
    FILE * maingrid_file;
    maingrid_file = fopen("main_grid.txt", "w");

    if (argc < 2) error("use: zero <image>\nfinds JPEG grids and forgeries");

    input = iio_read_image_double_split(argv[1], &X, &Y, &C);
    /* input2 = iio_read_image_double_split(argv[2], &X2, &Y2, &C2); */

    /* luminance image */
    image = (double *) xcalloc(X*Y, sizeof(double));
    /* image2 = (double *) xcalloc(X2*Y2, sizeof(double)); */

    /* vote map */
    votes = (int *) xcalloc(X * Y, sizeof(int));
    /* votes2 = (int *) xcalloc(X2 * Y2, sizeof(int)); */

    /* compute forged regions */
    forged_regions = (meaningful_reg *) xcalloc(X*Y, sizeof(meaningful_reg));
    /* forged_regions2 = (meaningful_reg *) xcalloc(X2*Y2, sizeof(meaningful_reg)); */

    /* compute forgery masks */
    forgery = (int *) xcalloc(X * Y, sizeof(int));
    forgery_e = (int *) xcalloc(X * Y, sizeof(int));
    /* forgery2 = (int *) xcalloc(X2 * Y2, sizeof(int)); */
    /* forgery_e2 = (int *) xcalloc(X2 * Y2, sizeof(int)); */

    /* run algorithm */
    main_grid = zero(input, image, votes, lnfa_grids, forged_regions, forgery, forgery_e, X, Y, C);
    /* zero2(input2, image2, votes2, forged_regions2, forgery2, forgery_e2, X2, Y2, C2); */

    /* store vote map and forgery detection outputs */
    iio_write_image_double("luminance.png", image, X, Y);

    fprintf(maingrid_file, "%d", main_grid);
    iio_write_image_int("votes.png", votes, X, Y);
    iio_write_image_int("forgery.png", forgery, X, Y);
    iio_write_image_int("forgery_c.png", forgery_e, X, Y);


    fclose(maingrid_file);

    /* free memory */
    free((void *) input);   //free((void *) input2);
    free((void *) image);   //free((void *) image2);
    free((void *) votes);   //free((void *) votes2);
    free((void *) forged_regions);  //free((void *) forged_regions2);
    free((void *) forgery);  //free((void *) forgery2);
    free((void *) forgery_e);   //free((void *) forgery_e2);

    return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
