// imported from matlab
/*
 * Copyright 1993-2002 The MathWorks, Inc.
 * $Revision: 1.8 $  $Date: 2002/04/24 03:00:15 $
 */

/*
 * This file contains the function body for the watershed algorithm.
 * By #defining TYPE and #including this file, a module can instantiate
 * this function for a particular numeric type.  #define DO_NAN_CHECK
 * before #including this file to enable runtime checks for NaNs, which
 * are not allowed.  DO_NAN_CHECK is only recommended if TYPE is float
 * or double.
 *
 * Algorithm reference: L. Vincent and P. Soille, "Watershed in digital 
 * spaces: An efficient algorithm based on immersion simulations," IEEE
 * Transactions on Pattern Analysis and Machine Intelligence, vol. 13,
 * n. 6, 1990, pp. 583-598.  That algorithm has been modified as noted
 * below.
 */

/*
 * void compute_watershed_vs(TYPE *I, int N, NeighborWalker_T walker,
 *                           double *sort_index, double *L)
 */

/*#include "ieeefp.h"*/

#include "watershed.h"

#define int32_T int 
#define true  1
#define false 0 
#define  MAX_int32_T     ((int32_T)(2147483647))    /* 2147483647  */
#define  MIN_int32_T     ((int32_T)(-2147483647-1)) /* -2147483648 */
#define  MAX_uint32_T    ((uint32_T)(0xFFFFFFFFU))  /* 4294967295  */
#define  MIN_uint32_T    ((uint32_T)(0))
 
void compute_watershed(float *I, int **coor, int N, int xdim, int ydim, unsigned long *sort_index, long *L)
{
  int current_label = 0;
  queueADT pixel_queue;
  int32_T *dist;
  int32_T closest_dist;
  int32_T closest_label_value;
  int closest_label_value_is_unique;
  int32_T fictitious = FICTITIOUS;
  int32_T wshed = WSHED;
  int k;
  int num_processed_pixels;
  int k1;
  int k2;
  int mask = MASK;
  int p;
  int q;
  int r;
  int current_distance;
  float current_level;
  
  
#ifdef DO_NAN_CHECK
  for (k = 0; k < N; k++)
    {
      if (isnanf(I[k]))
        {
	  printf("\n Images:watershed:expectedNonnan:%s",
		 "Input image may not contain NaNs.");
        }
    }
#endif /* DO_NAN_CHECK */
  
  /*
   * If the input array is empty, there's nothing to do here.
   */
  if (N == 0)
    {
      return;
    }
  
  /*
   * Initialize output array.
   */
  for (k = 0; k < N; k++)
    {
      L[k] = INIT;
    }
  
  /*
   * Initialize the pixel queue.
   */
  pixel_queue = QueueCreat();
  
  /*
   * Initialize the distance array, filling it with zeros via mxCalloc.
   */
  dist = (int32_T*)malloc(sizeof(int32_T)*N);//new int32_T [N];       
  
  num_processed_pixels = 0;
  while (num_processed_pixels < N)
    {
      /*
       * Find the next set of pixels that all have the same value.
       */
      k1 = num_processed_pixels;
      current_level = I[(int) sort_index[k1]];
      k2 = k1;
      do
        {
	  k2++;
        } 
      while ((k2 < N) && (I[(int) sort_index[k2]] == current_level));
      k2--;
      
      /*
       * Mask all image pixels whose value equals current_level.
       */
      for (k = k1; k <= k2; k++)
        {
	  p = (int) sort_index[k];
	  L[p] = mask;

	  int *nh = (int*) malloc(sizeof(int)*8);
	  int nhsize = 0;
	  int x,y;
	  for(y=-1; y<=1; y++)
	    for(x=-1; x<=1; x++)
	      {
		int nx = coor[p][0]+x;
		int ny = coor[p][1]+y;
		if((nx >= 0) && (nx < xdim)
		   && (ny >= 0) && (ny < ydim)
		   && ((x!=0)||(y!=0)))
		  {
		    nh[nhsize] = ny*xdim+nx;
		    nhsize += 1;
		  }
	      }
	  int i;
	  for(i=0; i<nhsize; i++)
	    {
	      q = nh[i];
	      if ((L[q] > 0) || (L[q] == wshed))
                {
		  /*
		   * Initialize queue with neighbors at current_level
		   * of current basins or watersheds.
		   */
		  dist[p] = 1;
		  QueueEnter(pixel_queue, p);
		  break;
                }
            }
	  num_processed_pixels++;
        }
      
      current_distance = 1;
      QueueEnter(pixel_queue, fictitious);
      
      /*
       * Extend the basins.
       */
      while (true)
        {
	  p = QueueDelete(pixel_queue);
	  if (p == fictitious)
            {
	      if (QueueLength(pixel_queue) == 0)
                {
		  break;
                }
	      else
                {
		  QueueEnter(pixel_queue, fictitious);
		  current_distance++;
		  p = QueueDelete(pixel_queue);
                }
            }
	  
	  /*
	   * NOTE: the code from here down to "detect and process
	   * new minima" is a modified version of the algorithm originally 
	   * published in Vincent and Soille's paper.  That algorithm
	   * could make several changes to L[p] during a single
	   * sweep of its neighbors, which sometimes results in incorrect
	   * labeling.  This seems to be particularly a problem in
	   * higher dimensions with the correspondingly larger number
	   * of neighbors.  Here the algorithm is changed to make a
	   * sweep of the neighborhood, accumulating key information
	   * about it configuration, and then, after the neighborhood
	   * sweep is finished, make one and only one change to L[p].
	   */
	  
	  /*
	   * Find the labeled or watershed neighbors with the closest
	   * distance.  At the same time, put any masked neighbors
	   * whose distance is 0 onto the queue and reset their distance
	   * to 1.
	   */
	  closest_dist = MAX_int32_T;
	  closest_label_value = 0;
	  closest_label_value_is_unique = true;
	  
	  int *nh = (int*)malloc(sizeof(int)*8);//new int[8];
	  int nhsize = 0;
	  int x,y;
	  for(y=-1; y<=1; y++)
	    for(x=-1; x<=1; x++)
	      {
		int nx = coor[p][0]+x;
		int ny = coor[p][1]+y;
		if((nx >= 0) && (nx < xdim)
		   && (ny >= 0) && (ny < ydim)
		   && ((x!=0)||(y!=0)))
		  {
		    nh[nhsize] = ny*xdim+nx;
		    nhsize += 1;
		  }
	      }
	  int i;
	  for(i=0; i<nhsize; i++)
	    {
	      q = nh[i];
            
	      if ((L[q] > 0) || (L[q] == WSHED))
                {
		  if (dist[q] < closest_dist)
                    {
		      closest_dist = dist[q];
		      if (L[q] > 0)
                        {
			  closest_label_value = L[q];
                        }
                    }
		  else if (dist[q] == closest_dist)
                    {
		      if (L[q] > 0)
                        {
			  if ((closest_label_value > 0) &&
			      (L[q] != closest_label_value))
                            {
			      closest_label_value_is_unique = false;
                            }
			  closest_label_value = L[q];
                        }
                    }
                }
	      
	      else if ((L[q] == MASK) && (dist[q] == 0))
                {
		  /*
		   * q is a plateau pixel.
		   */
		  dist[q] = current_distance + 1;
		  QueueEnter(pixel_queue, q);
                }
            }
	  
	  /*
	   * Label p.
	   */
	  if ((closest_dist < current_distance) && (closest_label_value > 0))
            {
	      if (closest_label_value_is_unique && 
		  ((L[p] == MASK) || (L[p] == WSHED)))
                {
		  L[p] = closest_label_value;
                }
	      else if (! closest_label_value_is_unique ||
		       (L[p] != closest_label_value))
                {
		  L[p] = WSHED;
                }
            }
	  else if (L[p] == MASK)
            {
	      L[p] = WSHED;
            }
        }
      
      /*
       * Detect and process new minima at current_level.
       */
      for (k = k1; k <= k2; k++)
        {
	  p = (int) sort_index[k];
	  dist[p] = 0;
	  if (L[p] == mask)
            {
	      /*
	       * p is inside a new minimum.
	       */
	      current_label++;  /* create a new label */
	      QueueEnter(pixel_queue, p);
	      L[p] = current_label;
	      while (QueueLength(pixel_queue) > 0)
                {
		  q = QueueDelete(pixel_queue);
		  
		  /*
		   * Inspect neighbors of q.
		   */
		  int *nh = (int*)malloc(sizeof(int)*8);//new int[8];
		  int nhsize = 0;
		  int x,y;
		  for(y=-1; y<=1; y++)
		    for(x=-1; x<=1; x++)
		      {
			int nx = coor[q][0]+x;
			int ny = coor[q][1]+y;
			if((nx >= 0) && (nx < xdim)
			   && (ny >= 0) && (ny < ydim)
			   && ((x!=0)||(y!=0)))
			  {
			    nh[nhsize] = ny*xdim+nx;
			    nhsize += 1;
			  }
		      }
		  int i;
		  for(i=0; i<nhsize; i++)
		    {
		      r = nh[i];
		      if (L[r] == mask)
                        {
			  QueueEnter(pixel_queue, r);
			  L[r] = current_label;
                        }
                    }
                }
            }
        }
    }
  
  //mxAssert(queue_length(pixel_queue) == 0, "");
  QueueDestroy(pixel_queue);
  free(dist);
}


void watershed(float **im, int xdim, int ydim, float **water, float sigma)
{
  
  long imsize = xdim*ydim;
  float *I = Falloc1d(imsize);//vector(imsize);
  int **coor = Ialloc2d(imsize,2);//imatrix(imsize,2);
  
  float **gx = Falloc2d(xdim,ydim);//matrix(xdim, ydim);
  float **gy = Falloc2d(xdim,ydim);//matrix(xdim, ydim);
  
  int i,j;
  for(j=0; j<ydim; j++)
    for(i=0; i<xdim; i++)
      {
	gx[i][j] = 0;
	gy[i][j] = 0;
      }

  float **sim = Falloc2d(xdim,ydim);//matrix(xdim,ydim);
  
  gaussian_smooth(im, xdim, ydim, sigma, sim);

  
  for(j=1; j<ydim-1; j++)
    for(i=1; i<xdim-1; i++)
      {
	gx[i][j] = sim[i+1][j] - sim[i][j];
	gy[i][j] = sim[i][j+1] - sim[i][j];
      }

  Ffree2d(sim,xdim);
  
  int index = 0;
  for(j=0; j<ydim; j++)
    for(i=0; i<xdim; i++)
      {
	I[index] = sqrt(gx[i][j]*gx[i][j] + gy[i][j]*gy[i][j]);
	coor[index][0] = i;
	coor[index][1] = j;
	index += 1;
      }
  
  Ffree2d(gx,xdim);
  Ffree2d(gy,xdim);

  unsigned long * sort_index = (unsigned long*)malloc(sizeof(unsigned long)*imsize);//lvector(imsize);
  indexx(imsize, I, sort_index);
  
  long *L = (long*)malloc(sizeof(long)*imsize);
  
  FILE *fp;
  fp=fopen("grad.txt","w");
  for(i=0; i<imsize; i++)
    fwrite(&I[i],sizeof(float),1,fp);
  fclose(fp);

  compute_watershed(I, coor, imsize, xdim, ydim, sort_index, L);
  
  for(i=0; i<imsize; i++)
    water[coor[i][0]][coor[i][1]] = L[i];
  
  free(I);
  free(sort_index);
  
  free(L);

  Ifree2d(coor,imsize);//free_imatrix(coor,imsize,2);
}


void watershed_post(float **im, float **water, float **result,int xdim, int ydim)
{
  float tmp;
  int num_region = 0;
  int i,j;
  for(j=0; j<ydim; j++)
    for(i=0; i<xdim; i++)
      {
	if(water[i][j]>num_region)
	  num_region = (int)water[i][j];
      }
  num_region += 1;
  
  float *mean = Falloc1d(num_region);
  float *size = Falloc1d(num_region);
  
  for(i=0; i<num_region; i++)
    mean[i] = size[i] = 0;
  
  for(j=0; j<ydim; j++)
    for(i=0; i<xdim; i++)
      {
	mean[(int)water[i][j]]+=im[i][j];
	size[(int)water[i][j]]+=1;
      }
  
  float max_pixel =0; 
  for(j=0; j<ydim; j++)
    for(i=0; i<xdim; i++)
      if(im[i][j]>max_pixel)
	max_pixel = im[i][j];

  //printf("\n the max pixel is %f",max_pixel);
  for(i=0; i<num_region;i++)
    {
      if(size[i] == 0)
	{
	  mean[i] = 0;
	  printf("\n the size[%d] = 0",i);
	}
      else
	mean[i]/=size[i];
    }
  mean[0] = 0;

  max_pixel =0; 
  for(i=0; i<num_region;i++)
    if(max_pixel<mean[i])
      max_pixel = mean[i];
  
  //  printf("\n the max mean is %f",max_pixel);
  
  for(j=0; j<ydim; j++)
    for(i=0; i<xdim; i++)
      result[i][j]=mean[(int)water[i][j]];
  
  free(mean);
  free(size);
}


void make_gaussian_kernel(float sigma, float **kernel, int *windowsize)
{
   int i, center;
   float x, fx, sum=0.0;

   *windowsize =(int)( 1 + 2 * ceil(2.5 * sigma));
   center = (int)((*windowsize) / 2);

   //printf("  The kernel has %d elements.\n", *windowsize);
   
   if((*kernel = Falloc1d(*windowsize)) == NULL){
     printf("Error callocing the gaussian kernel array.\n");
      exit(1);
   }

   for(i=0;i<(*windowsize);i++){
      x = (float)(i - center);
      fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
      (*kernel)[i] = fx;
      sum += fx;
   }

   for(i=0;i<(*windowsize);i++) (*kernel)[i] /= sum;
   
   
   //  printf("The filter coefficients are:\n");
   //for(i=0;i<(*windowsize);i++)
   //  printf("(*kernel)[%d] = %f\n", i,(*kernel)[i]);   
}

void gaussian_smooth(float **image, int xdim, int ydim, float sigma, float **smoothedim)
{
   int x, y, xx, yy,     /* Counter variables. */
     windowsize,              /* Dimension of the gaussian kernel. */
     center;                  /* Half of the windowsize. */
   float **tmpim,             /* Buffer for separable filter gaussian smoothing. */
     *kernel,              /* A one dimensional gaussian kernel. */
     dot,                  /* Dot product summing variable. */
     sum;                  /* Sum of the kernel weights variable. */

   /****************************************************************************
    * Create a 1-dimensional gaussian smoothing kernel.
    ****************************************************************************/
   // printf("   Computing the gaussian smoothing kernel.\n");
   make_gaussian_kernel(sigma, &kernel, &windowsize);
   center = windowsize / 2;
   //printf("\n center is %d\n",center);
   //printf("The filter coefficients are:\n");
   
   //for(int i=0;i<windowsize;i++)
   //  printf("kernel[%d] = %f\n", i, kernel[i]); 
   /****************************************************************************
    * Allocate a temporary buffer image and the smoothed image.
    ****************************************************************************/
   if((tmpim = Falloc2d(xdim,ydim)) == NULL)
     {
       printf("Error allocating the buffer image.\n");
       exit(1);
     }
 
   /****************************************************************************
   * Blur in the x - direction.
   ****************************************************************************/
   //printf("   Bluring the image in the X-direction.\n");
   
   int indx;
   for(y=0; y<ydim; y++) 
     for(x=0; x<xdim; x++)
       {
	 dot = 0.0;
	 sum = 0.0;
	 for(yy=(-center);yy<=center;yy++)
	   {
	     //  if(((y+yy) >= 0) && ((y+yy) < ydim))
	     indx = y+yy;
	     if(indx<0)
	       indx = - indx;
	     else if(indx>=ydim)
	       indx = 2*ydim - indx - 1;
	     
	     dot += image[x][indx] * kernel[center+yy];
	     sum += kernel[center+yy];
	     
	   }
	 tmpim[x][y] = dot/sum;
       }
   
   /****************************************************************************
    * Blur in the y - direction.
    ****************************************************************************/
   //printf("   Bluring the image in the Y-direction.\n");
 for(y=0; y<ydim; y++)
   for(x=0; x<xdim; x++)
     {
       sum = 0.0;
       dot = 0.0;
       for(xx=(-center);xx<=center;xx++)
	 {
	   //if(((x+xx) >= 0) && ((x+xx) < xdim))
	   indx = x+xx;
	   if(indx<0)
	     indx = -indx;
	   else if(indx>=xdim)
	     indx = 2*xdim- indx-1;
	   
	   dot += tmpim[indx][y] * kernel[center+xx];
	   sum += kernel[center+xx];
	   
	 }
       smoothedim[x][y] = dot/sum;
     }
 
 Ffree2d(tmpim,xdim);
 free(kernel);
}

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx(long n, float arr[], unsigned long indx[])
{
        long i,indxt,ir=n-1,itemp,j,k,l=0;
        int jstack=0,*istack;
        float a;

        istack= Ialloc1d(NSTACK);

        for (j=0;j<n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=0;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[l+1]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l+1];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l+1]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (jstack > NSTACK) printf("NSTACK too small in iindexx.");
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        
        

        
        free(istack);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
