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
#include <float.h>

#include "zero.h"

/*----------------------------------------------------------------------------*/
/* LN10 */
#ifndef M_LN10
#define M_LN10      2.30258509299404568401799145468436421
#endif /* !M_LN10 */

/*----------------------------------------------------------------------------*/
/* PI */
#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif /* !M_PI */

/*----------------------------------------------------------------------------*/
#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/* fatal error, print a message to standard-error output and exit.
 */
void error(char * msg) {
    fprintf(stderr,"error: %s\n",msg);
    exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* memory allocation, initialization to 0, print an error and exit if fail.
 */
void * xcalloc(size_t n_items, size_t size) {
    void * p;
    if (size == 0)
        error("xcalloc: zero size");
    p = calloc(n_items, size);
    if (p == NULL)
        error("xcalloc: out of memory");
    return p;
}

/*----------------------------------------------------------------------------*/
/* convert rgb image to luminance.
 */
void rgb2luminance(double * input, double * output, int X, int Y, int C) {
    if (C >= 3) {
        for (int x=0; x<X; x++)
            for(int y=0; y<Y; y++)
                output[x+y*X] = round(0.299 * input[x + y*X + 0*X*Y]
                                    + 0.587 * input[x + y*X + 1*X*Y]
                                    + 0.114 * input[x + y*X + 2*X*Y]);
    }
    else
        memcpy(output, input, X*Y*sizeof(double));
}

/*----------------------------------------------------------------------------*/
/* computes the logarithm of NFA to base 10.

   NFA = NT.b(n,k,p)
   the return value is log10(NFA)

   n,k,p - binomial parameters.
   logNT - logarithm of Number of Tests
*/
#define TABSIZE 100000
double log_nfa(int n, int k, double p, double logNT) {
    static double inv[TABSIZE];  /* table to keep computed inverse values */
    double tolerance = 0.1;      /* an error of 10% in the result is accepted */
    double log1term, term, bin_term, mult_term, bin_tail, err;
    double p_term = p / (1.0-p);

    if (n<0 || k<0 || k>n || p<0.0 || p>1.0)
        error("wrong n, k or p values in nfa()");

    if (n==0 || k==0) return logNT;
    if (n==k) return logNT + (double)n * log10(p);

    log1term = lgamma((double)n+1.0) - lgamma((double)k+1.0)
        - lgamma((double)(n-k)+1.0)
        + (double)k * log(p) + (double)(n-k) * log(1.0-p);

    term = exp(log1term);
    if (term == 0.0) {                       /* the first term is almost zero */
        if ((double)k > (double)n * p)    /* at beginning or end of the tail? */
            return log1term / M_LN10 + logNT; /* end: use just the first term */
        else
            return logNT;                     /* begin: the tail is roughly 1 */
    }

    bin_tail = term;
    for (int i=k+1; i<=n; i++) {
        bin_term = (double)(n-i+1) * (i<TABSIZE ?
                                      (inv[i] ? inv[i] : (inv[i]=1.0/(double)i))
                                      : 1.0/(double)i);
        mult_term = bin_term * p_term;
        term *= mult_term;
        bin_tail += term;
        if (bin_term<1.0) {
            /* when bin_term<1 then mult_term_j<mult_term_i for j>i.
               then, the error on the binomial tail when truncated at
               the i term can be bounded by a geometric serie of form
               term_i * sum mult_term_i^j.                            */
            err = term * (( 1.0 - pow(mult_term, (double)(n-i+1))) /
                          (1.0 - mult_term) - 1.0 );

            /* one wants an error at most of tolerance*final_result, or:
               tolerance * abs(-log10(bin_tail)-logNT).
               now, the error that can be accepted on bin_tail is
               given by tolerance*final_result divided by the derivative
               of -log10(x) when x=bin_tail. that is:
               tolerance * abs(-log10(bin_tail)-logNT) / (1/bin_tail)
               finally, we truncate the tail if the error is less than:
               tolerance * abs(-log10(bin_tail)-logNT) * bin_tail        */
            if (err < tolerance * fabs(-log10(bin_tail)-logNT) * bin_tail)
                break;
        }
    }

    return log10(bin_tail) + logNT;
}

/*----------------------------------------------------------------------------*/
/* computes the vote map.
 */
void compute_grid_votes_per_pixel(double * image, int * votes, int X, int Y) {
    double cos_t[8][8];
    int * zeros;  /* maximal number of zeros found for a given pixel */

    /* compute cosine table */
    for (int k=0; k<8; k++)
        for (int l=0; l<8; l++)
            cos_t[k][l] = cos((2.0 * k + 1.0) * l * M_PI / 16.0);

    /* initialize zero to 0 with calloc */
    zeros = (int *) xcalloc(X * Y, sizeof(int));

    /* initialize votes to -1 */
    for (int n=0; n<X*Y; n++) votes[n] = -1;

    /* compute DCT by 8x8 blocks */
#pragma omp parallel for
    for (int x=0; x<X-7; x++)
        for (int y=0; y<Y-7; y++) {
            int z = 0; /* number of zeros */
            int const_along_x = TRUE;
            int const_along_y = TRUE;

            /* check whether the block is constant along x or y axis */
            for (int xx=0; xx<8 && (const_along_x || const_along_y); xx++)
                for (int yy=0; yy<8 && (const_along_x || const_along_y); yy++) {
                    if (image[x+xx + (y+yy) * X] != image[x+0 + (y+yy) * X])
                        const_along_x = FALSE;
                    if (image[x+xx + (y+yy) * X] != image[x+xx + (y+0) * X])
                        const_along_y = FALSE;
                }

            /* compute DCT for 8x8 blocks staring at x,y and count its zeros */
            for (int i=0; i<8; i++)
                for (int j=0; j<8; j++)
                    if (i > 0 || j > 0) { /* coefficient 0 is not be counted */
                        double dct_ij = 0.0;

                        for (int xx=0; xx<8; xx++)
                            for (int yy=0; yy<8; yy++)
                                dct_ij += image[ x+xx + (y+yy) * X ]
                                    * cos_t[xx][i] * cos_t[yy][j];
                        dct_ij *= 0.25 * (i==0 ? 1.0/sqrt(2.0) : 1.0)
                            * (j==0 ? 1.0/sqrt(2.0) : 1.0);

                        /* the finest quantization in JPEG is to integer values.
                           in such case, the optimal threshold to decide if a
                           coefficient is zero or not is the midpoint between
                           0 and 1, thus 0.5 */
                        if (fabs(dct_ij) < 0.5) {
                            z++;
                        }
                    }

            /* check all pixels in the block and update votes */
#pragma omp critical
            for (int xx=x; xx<x+8; xx++)
                for (int yy=y; yy<y+8; yy++) {
                    /* if two grids are tied in number of zeros, do not vote. */
                    if (z == zeros[xx+yy*X])
                        votes[xx+yy*X] = -1;

                    /* update votes when the current grid has more zeros. */
                    if (z > zeros[xx+yy*X]) {
                        zeros[xx+yy*X] = z;
                        votes[xx+yy*X] = const_along_x || const_along_y ? -1
                            : (x % 8) + (y % 8) * 8;
                    }
                }
        }

    /* set pixels on the border to non valid votes - only pixels that
       belong to the 64 full 8x8 blocks inside the image can vote */
    for (int xx=0; xx<X; xx++) {
        for (int yy=0; yy<7; yy++)
            votes[xx+yy*X] = -1;
        for (int yy=Y-7; yy<Y; yy++)
            votes[xx+yy*X] = -1;
    }
    for (int yy=0; yy<Y; yy++) {
        for (int xx=0; xx<7; xx++)
            votes[xx+yy*X] = -1;
        for (int xx=X-7; xx<X; xx++)
            votes[xx+yy*X] = -1;
    }

    /* free memory */
    free((void *) zeros);
}

/*----------------------------------------------------------------------------*/
/* detect the main grid of the image from the vote map.
 */
int detect_global_grids(int * votes, double * lnfa_grids, int X, int Y) {
    double logNT = 2.0 * log10(64.0) + 2.0 * log10(X) + 2.0 * log10(Y);
    int grid_votes[64];
    int max_votes = 0;
    int most_voted_grid = -1;
    double p = 1.0 / 64.0;

    /* initialize counts */
    for (int i=0; i<64; i++) grid_votes[i] = 0;

    /* count votes per possible grid origin */
    for (int x=0; x<X; x++)
        for (int y=0; y<Y; y++)
            if (votes[x+y*X] >= 0 && votes[x+y*X] < 64) {
                int grid = votes[x+y*X];

                grid_votes[grid]++; /* increment votes */

                /* keep track of maximum of votes and the associated grid */
                if (grid_votes[grid] > max_votes) {
                    max_votes = grid_votes[grid];
                    most_voted_grid = grid;
                }
            }

    /* compute the NFA value for all the significant grids.  votes are
       correlated by irregular 8x8 blocks dividing by 64 gives a rough
       count of the number of independent votes */
    for (int i=0; i<64; i++) {
        int n = X * Y / 64;
        int k = grid_votes[i] / 64;
        lnfa_grids[i] = log_nfa(n, k, p, logNT);
    }

    /* meaningful grid -> main grid found! */
    if (most_voted_grid >= 0 && most_voted_grid < 64
        && lnfa_grids[most_voted_grid] < 0.0)
        return most_voted_grid;

    return -1;
}

/*----------------------------------------------------------------------------*/
/* detects zones which are inconsistent with a given grid.
 */
int detect_forgeries(int * votes, int * forgery_mask, int * forgery_mask_reg,
                     meaningful_reg * forged_regions,
                     int X, int Y, int grid_to_exclude, int grid_max) {
    double logNT = 2.0 * log10(64.0) + 2.0 * log10(X) + 2.0 * log10(Y);
    double p = 1.0 / 64.0;
    int forgery_n = 0;  /* number of forgeries found */
    int * mask_aux;
    int * used;
    int * reg_x;
    int * reg_y;
    int W = 9; /* Distance to look for neighbors in the region growing process.
                  A meaningful forgery must have a density of votes of at least
                  1/64. Thus, its votes should not be in mean further away one
                  from another than a distance of 8. One could use a little
                  more to allow for some variation in the distribution. */

    /* minimal block size that can lead to a meaningful detection */
    int min_size = ceil( 64.0 * logNT / log10( 64.0 ) );

    /* initialize global data */
    mask_aux = (int *) xcalloc(X * Y, sizeof(int));
    used = (int *) xcalloc(X * Y, sizeof(int));
    reg_x = (int *) xcalloc(X * Y, sizeof(int));
    reg_y = (int *) xcalloc(X * Y, sizeof(int));

    for (int i=0; i<X*Y; i++) used[i] = FALSE;

    /* region growing of zones that voted for other than the main grid */
    for (int x=0; x<X; x++)
        for (int y=0; y<Y; y++)
            if (used[x+y*X] == FALSE && votes[x+y*X] != grid_to_exclude
                && votes[x+y*X] >= 0 && votes[x+y*X] <= grid_max) {
                /* initialize region with the seed pixel */
                int reg_size = 0;
                int grid = votes[x+y*X];
                int x0 = x; /* region bounding box */
                int y0 = y;
                int x1 = x;
                int y1 = y;
                used[x+y*X] = TRUE;
                reg_x[reg_size] = x;
                reg_y[reg_size] = y;
                ++reg_size;

                /* iteratively add neighbor pixel of pixels in the region */
                for (int i=0; i<reg_size; i++)
                    for (int xx=reg_x[i]-W; xx<=reg_x[i]+W; xx++)
                        for (int yy=reg_y[i]-W; yy<=reg_y[i]+W; yy++)
                            if (xx >=0 && xx < X && yy >= 0 && yy < Y)
                                if (used[xx+yy*X] == FALSE
                                    && votes[xx+yy*X] == grid) {
                                    used[xx+yy*X] = TRUE;
                                    reg_x[reg_size] = xx;
                                    reg_y[reg_size] = yy;
                                    reg_size++;
                                    if (xx < x0) x0 = xx; /* update       */
                                    if (yy < y0) y0 = yy; /* bounding box */
                                    if (xx > x1) x1 = xx;
                                    if (yy > y1) y1 = yy;
                                }

                /* compute NFA for regions with at least the minimal size */
                if (reg_size >= min_size) {
                    int n = (x1 - x0 + 1) * (y1 - y0 + 1) / 64;
                    int k = reg_size / 64;
                    double lnfa = log_nfa( n, k, p, logNT );

                    if (lnfa < 0.0) { /* meaningful different grid found */
                        forged_regions[forgery_n].x0 = x0;
                        forged_regions[forgery_n].x1 = x1;
                        forged_regions[forgery_n].y0 = y0;
                        forged_regions[forgery_n].y1 = y1;
                        forged_regions[forgery_n].grid = grid;
                        forged_regions[forgery_n].lnfa = lnfa;

                        ++forgery_n;

                        /* mark points of the region in the forgery mask */
                        for (int i=0; i<reg_size; i++)
                            forgery_mask[reg_x[i] + reg_y[i]*X] = 255;
                    }
                }
            }

    /* regularized forgery mask by a morphologic closing operator */
    for (int x=W; x<X-W; x++)
        for (int y=W; y<Y-W; y++)
            if (forgery_mask[x+y*X] != 0)
                for (int xx=x-W; xx<=x+W; xx++)
                    for (int yy=y-W; yy<=y+W; yy++)
                        mask_aux[xx+yy*X] = forgery_mask_reg[xx+yy*X] = 255;
    for (int x=W; x<X-W; x++)
        for (int y=W; y<Y-W; y++)
            if (mask_aux[x+y*X] == 0)
                for (int xx=x-W; xx<=x+W; xx++)
                    for (int yy=y-W; yy<=y+W; yy++)
                        forgery_mask_reg[xx+yy*X] = 0;

    /* free memory */
    free((void *) mask_aux);
    free((void *) used);
    free((void *) reg_x);
    free((void *) reg_y);

    return forgery_n;
}

/*----------------------------------------------------------------------------*/
/* Zero algorithm
 */
int zero(double * input, double * input_jpeg,
         double * luminance, double * luminance_jpeg,
         int * votes, int * votes_jpeg,
         double * lnfa_grids,
         meaningful_reg * foreign_regions, int * foreign_regions_n,
         meaningful_reg * missing_regions, int * missing_regions_n,
         int * mask_f, int * mask_f_reg, int * mask_m, int * mask_m_reg,
         int X, int Y, int C, int C_jpeg) {

    int main_grid = -1;

    /* luminance image */
    rgb2luminance(input, luminance, X, Y, C);

    /* compute vote map */
    compute_grid_votes_per_pixel(luminance, votes, X, Y);

    /* detect global grids and main_grid */
    main_grid = detect_global_grids(votes, lnfa_grids, X, Y);

    /* compute forged regions */
    *foreign_regions_n = detect_forgeries(votes, mask_f, mask_f_reg,
                                          foreign_regions, X, Y,
                                          main_grid, 63);

    /* if a global grid is found and the JPEG QF-99 version is provided,
       try to found regions with missing JPEG grid */
    if (main_grid > -1 && input_jpeg != NULL) {
        /* luminance image */
        rgb2luminance(input_jpeg, luminance_jpeg, X, Y, C_jpeg);

        /* compute vote map */
        compute_grid_votes_per_pixel(luminance_jpeg, votes_jpeg, X, Y);

        /* update votemap by avoiding the votes for the main grid */
        for (int x=0; x<X; x++)
            for (int y=0; y<Y; y++)
                if (votes[x+y*X] == main_grid)
                    votes_jpeg[x+y*X] = -1;

        /* Try to detect an imposed JPEG grid.  No grid is to be excluded
           and we are interested only in grid with origin (0,0), so:
           grid_to_exclude = -1 and grid_max = 0 */
        *missing_regions_n = detect_forgeries(votes_jpeg, mask_m, mask_m_reg,
                                              missing_regions, X, Y, -1, 0);
    }

    return main_grid;
}
/*----------------------------------------------------------------------------*/
