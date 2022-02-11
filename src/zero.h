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

#ifndef ZERO_H
#define ZERO_H

/* Structure for local forgeries */
typedef struct {
    int x0, y0, x1, y1;
    int grid;
    double lnfa;
}  meaningful_reg;

void error(char * msg);

void * xcalloc(size_t n_items, size_t size);

void rgb2luminance(double * input, double * output, int X, int Y, int C);

double log_nfa(int n, int k, double p, double logNT);

void compute_grid_votes_per_pixel(double * image, int * votes, int X, int Y);

int detect_global_grids(int * votes, double * lnfa_grids, int X, int Y);

int detect_forgeries(int * votes, int * forgery_mask, int * forgery_mask_reg,
                     meaningful_reg * forged_regions,
                     int X, int Y, int grid_to_exclude, int grid_max);

int zero(double * input, double * input_jpeg,
         double * luminance, double * luminance_jpeg,
         int * votes, int * votes_jpeg,
         double * lnfa_grids,
         meaningful_reg * foreign_regions, int * foreign_regions_n,
         meaningful_reg * missing_regions, int * missing_regions_n,
         int * mask_f, int * mask_f_reg, int * mask_m, int * mask_m_reg,
         int X, int Y, int C, int C_jpeg);

#endif
/*----------------------------------------------------------------------------*/
