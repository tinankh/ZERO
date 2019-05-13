/*----------------------------------------------------------------------------

  Copyright (c) 2018-2019 Rafael Grompone von Gioi <grompone@gmail.com>
  Copyright (c) 2018-2019 Tina Nikoukhah <nikoukhah@cmla.ens-cachan.fr>
  Copyright (c) 2018-2019 Jérémy Anger <anger@cmla.ens-cachan.fr>
  Copyright (c) 2018-2019 Thibaud Ehret <ehret@cmla.ens-cachan.fr>

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
#include "iio.h"

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
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

/*----------------------------------------------------------------------------*/
/* fatal error, print a message to standard-error output and exit.
 */
void error(char * msg)
{
  fprintf(stderr,"error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* memory allocation, print an error and exit if fail.
 */
void * xmalloc(size_t size)
{
  void * p;
  if( size == 0 ) error("xmalloc: zero size");
  p = malloc(size);
  if( p == NULL ) error("xmalloc: out of memory");
  return p;
}

/*----------------------------------------------------------------------------*/
double * rgb2luminance(double * input, int X, int Y, int C)
{
  double * output;
  int x,y;

  if( C >= 3 )
    {
      output = xmalloc( X * Y * sizeof(double) );
      for(x=0; x<X; x++)
      for(y=0; y<Y; y++)
        output[x+y*X] = 0.299 * input[x + y*X + 0*X*Y]
                      + 0.587 * input[x + y*X + 1*X*Y]
                      + 0.114 * input[x + y*X + 2*X*Y];
    }
  else output = input;

  return output;
}

/*----------------------------------------------------------------------------*/
/* computes the logarithm of NFA to base 10.

   NFA = NT.b(n,k,p)
   the return value is log10(NFA)

   n,k,p - binomial parameters.
   logNT - logarithm of Number of Tests
 */
#define TABSIZE 100000
double log_nfa(int n, int k, double p, double logNT)
{
  static double inv[TABSIZE];   /* table to keep computed inverse values */
  double tolerance = 0.1;       /* an error of 10% in the result is accepted */
  double log1term,term,bin_term,mult_term,bin_tail,err;
  double p_term = p / (1.0-p);
  int i;

  if( n<0 || k<0 || k>n || p<0.0 || p>1.0 )
    error("wrong n, k or p values in nfa()");

  if( n==0 || k==0 ) return logNT;
  if( n==k ) return logNT + (double)n * log10(p);

  log1term = lgamma((double)n+1.0) - lgamma((double)k+1.0)
           - lgamma((double)(n-k)+1.0)
           + (double)k * log(p) + (double)(n-k) * log(1.0-p);

  term = exp(log1term);
  if( term == 0.0 )                        /* the first term is almost zero */
    {
      if( (double)k > (double)n * p )      /* at begining or end of the tail? */
        return log1term / M_LN10 + logNT;  /* end: use just the first term */
      else
        return logNT;                      /* begin: the tail is roughly 1 */
    }

  bin_tail = term;
  for(i=k+1;i<=n;i++)
    {
      bin_term = (double)(n-i+1) * ( i<TABSIZE ?
                   (inv[i] ? inv[i] : (inv[i]=1.0/(double)i)) : 1.0/(double)i );
      mult_term = bin_term * p_term;
      term *= mult_term;
      bin_tail += term;
      if(bin_term<1.0)
        {
          /* when bin_term<1 then mult_term_j<mult_term_i for j>i.
             then, the error on the binomial tail when truncated at
             the i term can be bounded by a geometric serie of form
             term_i * sum mult_term_i^j.                            */
          err = term * ( ( 1.0 - pow( mult_term, (double)(n-i+1) ) ) /
                         (1.0-mult_term) - 1.0 );

          /* one wants an error at most of tolerance*final_result, or:
             tolerance * abs(-log10(bin_tail)-logNT).
             now, the error that can be accepted on bin_tail is
             given by tolerance*final_result divided by the derivative
             of -log10(x) when x=bin_tail. that is:
             tolerance * abs(-log10(bin_tail)-logNT) / (1/bin_tail)
             finally, we truncate the tail if the error is less than:
             tolerance * abs(-log10(bin_tail)-logNT) * bin_tail        */
          if( err < tolerance * fabs(-log10(bin_tail)-logNT) * bin_tail ) break;
        }
    }
  return log10(bin_tail) + logNT;
}

/*----------------------------------------------------------------------------*/
int * compute_grid_votes_per_pixel(double * image, int X, int Y)
{
  double cos_t[8][8];
  int * zeros;
  int * votes;
  int x,y,k,l,n;

  /* compute cosine table */
  for(k=0; k<8; k++)
  for(l=0; l<8; l++)
    cos_t[k][l] = cos( (2.0 * k + 1.0) * l * M_PI / 16.0 );

  /* initialize zeros and votes */
  zeros = (int *) xmalloc( X * Y * sizeof(int) );
  votes = (int *) xmalloc( X * Y * sizeof(int) );
  for(n=0; n<X*Y; n++) zeros[n] = votes[n] = -1;

  /* compute DCT by 8x8 blocks */
#pragma omp parallel for private(x,y)
  for(x=0; x<X-7; x++)
  for(y=0; y<Y-7; y++)
    {
      int z = 0; /* number of zeros */
      int xx,yy,i,j;

      /* compute DCT for the 8x8 block staring at x,y and count its zeros */
      for(i=0; i<8; i++)
      for(j=0; j<8; j++)
      if( i > 0 || j > 0 ) /* the coefficient 0 should not be counted */
        {
          double dct_ij = 0.0;

          for(xx=0; xx<8; xx++)
          for(yy=0; yy<8; yy++)
            dct_ij += image[ x+xx + (y+yy) * X ] * cos_t[xx][i] * cos_t[yy][j];
          dct_ij *= 0.25 * ( i==0 ? 1.0/sqrt(2.0) : 1.0 )
                         * ( j==0 ? 1.0/sqrt(2.0) : 1.0 );

          /* the finest quantization in JPEG is to integer values.
             in such case, the optimal threshold to decide if a
             coefficient is zero or not is the midpoint between
             0 and 1, thus 0.5 */
          if( fabs(dct_ij) < 0.5 ) ++z;
        }

      /* check all pixels in the block and update votes */
#pragma omp critical
      for(xx=x; xx<x+8; xx++)
      for(yy=y; yy<y+8; yy++)
        {
          /* if two grids are tied in number of zeros, do not vote */
          if( z == zeros[xx+yy*X] ) votes[xx+yy*X] = -1;

          /* update votes when the current grid has more zeros */
          if( z > zeros[xx+yy*X] )
            {
              zeros[xx+yy*X] = z;
              votes[xx+yy*X] = (x % 8) + (y % 8) * 8;
            }
        }
    }

  /* store zeros and votes */
  iio_write_image_int("votes.png",votes,X,Y);

  /* free memory */
  free( (void *) zeros );

  return votes;
}

/*----------------------------------------------------------------------------*/
int detect_main_grid(int * votes, int X, int Y)
{
  double logNT = log10(64.0) + 1.5 * log10(X) + 1.5 * log10(Y);
  int grid_votes[64];
  int max_votes = 0;
  int most_voted_grid = -1;
  double p = 1.0 / 64.0;
  double lnfa;
  int x,y,i,n,k;

  /* initialize counts */
  for(i=0; i<64; i++) grid_votes[i] = 0;

  /* count votes per possible grid origine */
  for(x=0; x<X; x++)
  for(y=0; y<Y; y++)
    if( votes[x+y*X] >= 0 && votes[x+y*X] < 64 )
      {
        int grid = votes[x+y*X];

        ++grid_votes[grid]; /* increment votes */

        /* keep track of maximum of votes and the associated grid */
        if( grid_votes[grid] > max_votes )
          {
            max_votes = grid_votes[grid];
            most_voted_grid = grid;
          }
      }

  /* compute the NFA value of the most voted grid.  votes are
     correlated by irregular 8x8 blocks dividing by 64 gives a rough
     count of the number of independent votes */
  n = X * Y / 64;
  k = grid_votes[most_voted_grid] / 64;
  lnfa = log_nfa( n, k, p, logNT );

  /* meaningful grid -> main grid found! */
  if( lnfa < 0.0 )
    {
      printf("main grid: #%d [%d %d] log(nfa) = %g\n", most_voted_grid,
             most_voted_grid % 8, most_voted_grid / 8, lnfa);
      return most_voted_grid;
    }

  /* main grid not found */
  printf("no overall JPEG grid found\n");
  return -1;
}

/*----------------------------------------------------------------------------*/
void detect_forgery(int * votes, int X, int Y, int main_grid)
{
  double logNT = log10(64.0) + 1.5 * log10(X) + 1.5 * log10(Y);
  double p = 1.0 / 64.0;
  int * forgery;
  int * forgery_d;
  int * forgery_e;
  int * used;
  int * reg_x;
  int * reg_y;
  int x,y,xx,yy,i;
  int W = 12; /* distance to look for neighbors in the region growing process.
                 A meaningful forgery must have a density of votes of at least
                 1/64. thus, its votes should not be in mean further away one
                 from another than a distance of 8.
                 we use a little more to allow for some variation in the
                 distribution. */

  /* miminal block size that can lead to a meaningful detection */
  int min_size = ceil( 64.0 * logNT / log10( 64.0 ) );

  /* initialize global data */
  forgery = (int *) xmalloc( X * Y * sizeof(int) );
  forgery_d = (int *) xmalloc( X * Y * sizeof(int) );
  forgery_e = (int *) xmalloc( X * Y * sizeof(int) );
  used = (int *) xmalloc( X * Y * sizeof(int) );
  reg_x = (int *) xmalloc( X * Y * sizeof(int) );
  reg_y = (int *) xmalloc( X * Y * sizeof(int) );
  for(i=0; i<X*Y; i++) forgery[i] = 0;
  for(i=0; i<X*Y; i++) forgery_d[i] = 0;
  for(i=0; i<X*Y; i++) forgery_e[i] = 0;
  for(i=0; i<X*Y; i++) used[i] = FALSE;

  /* region growing of zones that voted for other than the main grid */
  for(x=0; x<X; x++)
  for(y=0; y<Y; y++)
    if( used[x+y*X] == FALSE && votes[x+y*X] != main_grid && votes[x+y*X] >= 0 )
      {
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
        for(i=0; i<reg_size; i++)
          for(xx=reg_x[i]-W; xx<=reg_x[i]+W; xx++)
          for(yy=reg_y[i]-W; yy<=reg_y[i]+W; yy++)
            if( xx >=0 && xx < X && yy >= 0 && yy < Y )
              if( used[xx+yy*X] == FALSE && votes[xx+yy*X] == grid )
                {
                  used[xx+yy*X] = TRUE;
                  reg_x[reg_size] = xx;
                  reg_y[reg_size] = yy;
                  ++reg_size;
                  if( xx < x0 ) x0 = xx; /* update region bounding box */
                  if( yy < y0 ) y0 = yy;
                  if( xx > x1 ) x1 = xx;
                  if( yy > y1 ) y1 = yy;
                }

        /* compute NFA for regions with at least the minimal size */
        if( reg_size >= min_size )
          {
            int N = MAX(x1 - x0 + 1, y1 - y0 + 1);
            int n = N * N / 64;
            int k = reg_size / 64;
            double lnfa = log_nfa( n, k, p, logNT );

            if( lnfa < 0.0 ) /* meaningful grid different from the main found */
              {
                printf("forgery found: %d %d - %d %d [%dx%d] ",
                       x0,y0,x1,y1,x1-x0+1,y1-y0+1);
                printf("grid: #%d [%d %d] ", grid, grid % 8, grid / 8 );
                printf("n %d k %d log(nfa) = %g\n",n,k,lnfa);

                /* mark points of the region in the forgery mask */
                for(i=0; i<reg_size; i++)
                  forgery[reg_x[i] + reg_y[i]*X] = 255;
              }
          }
      }

  /* morphologic closing of forgery mask */
  for(x=W; x<X-W; x++)
  for(y=W; y<Y-W; y++)
    if( forgery[x+y*X] != 0 )
      for(xx=x-W; xx<=x+W; xx++)
      for(yy=y-W; yy<=y+W; yy++)
        forgery_d[xx+yy*X] = forgery_e[xx+yy*X] = 255;
  for(x=W; x<X-W; x++)
  for(y=W; y<Y-W; y++)
    if( forgery_d[x+y*X] == 0 )
      for(xx=x-W; xx<=x+W; xx++)
      for(yy=y-W; yy<=y+W; yy++)
        forgery_e[xx+yy*X] = 0;

  /* store forgery detection outputs */
  iio_write_image_int("forgery.png",forgery,X,Y);
  iio_write_image_int("forgery_c.png",forgery_e,X,Y);

  /* free memory */
  free( (void *) forgery );
  free( (void *) forgery_d );
  free( (void *) forgery_e );
  free( (void *) used );
  free( (void *) reg_x );
  free( (void *) reg_y );
}

/*----------------------------------------------------------------------------*/
/*                                    Main                                    */
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
  double * input;
  double * image;
  int * votes;
  int X,Y,C;
  int main_grid;

  if( argc != 2 ) error("use: zero <image>\nfinds JPEG grid and forgeries");

  input = iio_read_image_double_split(argv[1], &X, &Y, &C);

  image = rgb2luminance(input,X,Y,C);

  votes = compute_grid_votes_per_pixel(image,X,Y);

  main_grid = detect_main_grid(votes,X,Y);

  detect_forgery(votes,X,Y,main_grid);

  /* free memory */
  free( (void *) input );
  if( C >= 3 ) free( (void *) image );
  free( (void *) votes );

  return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
