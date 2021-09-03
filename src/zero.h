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

int detect_foreign_grids(int * votes, int * forgery, int * forgery_e,
                   meaningful_reg * forged_regions,
                   int X, int Y, int main_grid);

int detect_missing_grid(int * votes, int * forgery, int * forgery_ext,
                   meaningful_reg * forged_regions, int X, int Y);

void region_growing(int * votes, int X, int Y, int x0, int x1,
                    int y0, int y1);


int zero(double * input, double * image, int * votes, double * lnfa_grids,
         meaningful_reg * forged_regions, int * forgery, int *forgery_e,
         int X, int Y, int C);

int zero_bis(double * input, double * image, int * votes,
         meaningful_reg * forged_regions, int * forgery, int *forgery_e,
         int X, int Y, int C);


#endif
