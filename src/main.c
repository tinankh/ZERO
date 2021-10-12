/*----------------------------------------------------------------------------

  Copyright (c) 2018-2021 Rafael Grompone von Gioi <grompone@gmail.com>
  Copyright (c) 2018-2021 Jérémy Anger <anger@cmla.ens-cachan.fr>
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

int main(int argc, char ** argv) {
    double * input; double * input_compressed;
    double * image; double * image_compressed;
    int X, Y, C;
    int * votes; int * votes_compressed;
    double lnfa_grids[64] = {0.0};
    meaningful_reg * foreign_regions;
    meaningful_reg * missing_regions;
    int * f; int * mask_f;
    int * m; int * mask_m;
    int main_grid = -1;

    if (argc < 3) error("use: zero <image> <imagecompressed>\nfinds JPEG grids and forgeries");

    input = iio_read_image_double_split(argv[1], &X, &Y, &C);
    input_compressed = iio_read_image_double_split(argv[2], &X, &Y, &C);

    /* luminance image */
    image = (double *) xcalloc(X*Y, sizeof(double));
    image_compressed = (double *) xcalloc(X*Y, sizeof(double));

    /* vote map */
    votes = (int *) xcalloc(X * Y, sizeof(int));
    votes_compressed = (int *) xcalloc(X * Y, sizeof(int));

    /* compute forged regions */
    foreign_regions = (meaningful_reg *) xcalloc(X*Y, sizeof(meaningful_reg));
    missing_regions = (meaningful_reg *) xcalloc(X*Y, sizeof(meaningful_reg));

    /* compute forgery masks */
    f = (int *) xcalloc(X * Y, sizeof(int));
    mask_f = (int *) xcalloc(X * Y, sizeof(int));

    m = (int *) xcalloc(X * Y, sizeof(int));
    mask_m = (int *) xcalloc(X * Y, sizeof(int));

    /* run algorithm */
    main_grid = zero(input, input_compressed, image, image_compressed, votes, votes_compressed,
                     lnfa_grids, foreign_regions, missing_regions, f, mask_f, m, mask_m,  X, Y, C);

    /* store vote map and forgery detection outputs */
    iio_write_image_double("luminance.png", image, X, Y);
    iio_write_image_int("votes.png", votes, X, Y);
    iio_write_image_int("votes_compressed.png", votes_compressed, X, Y);
    iio_write_image_int("mask_f.png", mask_f, X, Y);
    iio_write_image_int("mask_m.png", mask_m, X, Y);

    /* free memory */
    free((void *) input);
    free((void *) input_compressed);

    free((void *) image);
    free((void *) image_compressed);

    free((void *) votes);
    free((void *) votes_compressed);

    free((void *) foreign_regions);
    free((void *) missing_regions);

    free((void *) f);
    free((void *) mask_f);
    free((void *) m);
    free((void *) mask_m);

    return EXIT_SUCCESS;
}
/* ---------------------------------------------------------------------------- */
