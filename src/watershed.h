/*
 * Copyright 1993-2002 The MathWorks, Inc.
 * $Revision: 1.8 $  $Date: 2002/03/15 15:58:42 $
 */

/*
 * WATERSHED_VS MEX-file
 *
 * L = WATERSHED_VS(A,CONN,IDX) finds the watershed regions of A, producing 
 * a label array L.  IDX is a zero-based sort index for A, which can be 
 * computed in MATLAB as follows:
 *    [junk,idx] = sort(A(:));
 *    idx = idx - 1;
 *
 * Input-output specs
 * ------------------
 * A      N-D full, real, numeric array
 *        empty allowed
 *        +/- Inf allowed
 *        NaN not allowed
 *        logical ignored
 *
 * CONN   See connectivity spec
 *
 * IDX    Double vector containing zero-based sort permutation index for A(:).
 *        This input is *not* checked to make sure that is a valid
 *        sort permutation index.  If it isn't valid, the results will
 *        be unpredictable.
 *
 * L      Full double array of the same size as A.
 */

#include "queue.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mvcd.h"


#define INIT -1
#define MASK -2
#define WSHED 0
#define FICTITIOUS -1

/*
 * Instantiate compute_watershed_vs_TYPE functions.
 */

#undef TYPE
#define DO_NAN_CHECK 0

void compute_watershed(float *I, int **coor, int N, int xdim, int ydim, unsigned long *sort_index, long *L);

void watershed(float **im, int xdim, int ydim, float **water, float sigma);
void watershed_post(float **im, float **water, float **result,int xdim, int ydim);
void gaussian_smooth(float **image, int xdim, int ydim, float sigma, float **smoothedim);

void indexx(long n, float arr[], unsigned long indx[]);
