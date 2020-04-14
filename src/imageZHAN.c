#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "general.h"
#include "mvcd.h"
#include "imageZHAN.h"

void Convolution1D(double* input,int inputSize,double* mask,int maskSize,double* output)
{
  int i,ii;
  int halfSize,evenOrNot;
  double sum;

  halfSize=maskSize/2;
  evenOrNot=(maskSize+1)%2;
  
  for (i=0;i<inputSize;i++)
    {
      sum=0;
      if (i-(halfSize-evenOrNot)>=0&&i+halfSize<inputSize)
	{
	  for (ii=-(halfSize-evenOrNot);ii<=halfSize;ii++)
	    sum+=input[i+ii]*mask[ii+halfSize-evenOrNot];
	}
      
      output[i]=sum;	    
 
    }  
}

void Convolution2DSeparate(double** input,Ivector2d inputSize,double* mask,int maskSize,double** output)
{
  int row,col;
  double ** temp;
  double *tempInput,*tempOutput;

  //Convolution along horizontal direction
  tempInput=Dalloc1d(inputSize.y);
  tempOutput=Dalloc1d(inputSize.y);

  temp=Dalloc2d(inputSize.x,inputSize.y);
  
  for (row=0;row<inputSize.x;row++)
    for (col=0;col<inputSize.y;col++)
      temp[row][col]=input[row][col];

  for (row=0;row<inputSize.x;row++)
    {      
      for (col=0;col<inputSize.y;col++)
	tempInput[col]=temp[row][col];
      
      Convolution1D(tempInput,inputSize.y,mask,maskSize,tempOutput);

      for (col=0;col<inputSize.y;col++)
	temp[row][col]=tempOutput[col];
    }

  free(tempInput);
  free(tempOutput);
  
  //Convolution along vertical direction

  tempInput=Dalloc1d(inputSize.x);
  tempOutput=Dalloc1d(inputSize.x);  
  for (col=0;col<inputSize.y;col++)
    {
      for (row=0;row<inputSize.x;row++)
	tempInput[row]=temp[row][col];

      Convolution1D(tempInput,inputSize.x,mask,maskSize,tempOutput);

      for (row=0;row<inputSize.x;row++)
	output[row][col]=tempOutput[row];
    }

  Dfree2d(temp,inputSize.x);
  free(tempInput);
  free(tempOutput);
}

void Convolution2DGeneral(double** input,Ivector2d inputSize,double** mask,Ivector2d maskSize,double** output)
{
  int row,col,ii,jj;
  Ivector2d halfSize,evenOrNot;
  double sum;

  halfSize.x=maskSize.x/2;
  halfSize.y=maskSize.y/2;

  evenOrNot.x=(maskSize.x+1)%2;
  evenOrNot.y=(maskSize.y+1)%2;

  for (row=0;row<inputSize.x;row++)
    for (col=0;col<inputSize.y;col++)
      output[row][col]=0;

  for (row=halfSize.x-evenOrNot.x;row<inputSize.x-halfSize.x;row++)
    for (col=halfSize.y-evenOrNot.y;col<inputSize.y-halfSize.y;col++)
      {
	sum=0;
	for (ii=-(halfSize.x-evenOrNot.x);ii<=halfSize.x;ii++)
	  for (jj=-(halfSize.y-evenOrNot.y);jj<=halfSize.y;jj++)
	    sum+=input[row+ii][col+jj]*mask[ii+halfSize.x-evenOrNot.x][jj+halfSize.y-evenOrNot.y];
	
	output[row][col]=sum;
      }
}

double Gaussian(double x, double sigma)
{
  return(exp((-x*x)/(2*sigma*sigma))/(sigma*sqrt(2*PIE)));
}

#define MAX_WIDTH 1000
#define GAUSSIAN_THRE 0.005

void Smoothing2DImageWithGaussianFilter(double** inputImage,double** outputImage,Ivector2d imageSize,double smoothSigma)
{    
  int i;

  //Determine gaussian filter size
  int filterWidth;
  for(i=0;i<MAX_WIDTH; i++) 	
    {
      if (Gaussian((double)i,smoothSigma)*(smoothSigma*sqrt(2*PIE))<GAUSSIAN_THRE)
	{
	  if (i>0)
	    filterWidth=(2*(i-1))+1;
	  else
	    printf("\n\n\nSomething abnormal in function 'Smoothing2DImageWithGaussianFilter'!\n\n\n");

 	  break;
	}
    }
  //printf("filter_width=%d\n", filterWidth);

  //Contructing Gaussian filter
  double* gaussianFilter1D;
  int halfWidth;

  gaussianFilter1D=Dalloc1d(filterWidth);
  halfWidth=filterWidth/2;

  for (i=-halfWidth;i<=halfWidth;i++)
    gaussianFilter1D[i+halfWidth]=Gaussian((double)i,smoothSigma);
  
  //Convolute Gaussian filter with image  
  Convolution2DSeparate(inputImage,imageSize,gaussianFilter1D,filterWidth,outputImage);
}

void DetectUCImageSize(char* filename,int x_size,int y_size,int* z_size)
{
  FILE* fp;

  fp=myopen(filename,"r");
  fseek(fp,0,SEEK_END);
  *z_size=ftell(fp)/(x_size*y_size);

  fclose(fp);
}


unsigned char ConvertDoubleToChar(double input,double translation,double scale)
{
  double temp;
  unsigned char output;

  temp=(input-translation)/scale;
  
  temp=(temp+1.0)*128;

  if (temp<0)
    output=0;
  else if (temp>255)
    output=255;
  else
    output=(unsigned char)(temp+0.5);

  return output;
}

void Convert3DDoubleImageToCharImage(double *** inputImage,unsigned char *** outputImage,Ivector3d imageSize)
{
  int i,j,k;
  double mean,std;

  mean=0;

  for (k=0;i<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
      for (j=0;j<imageSize.y;j++)
	mean+=inputImage[k][i][j]/(double)(imageSize.z*imageSize.x*imageSize.y);
  
  std=0;
  for (k=0;i<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
      for (j=0;j<imageSize.y;j++)
	std+=((inputImage[k][i][j]-mean)*(inputImage[k][i][j]-mean))/(double)(imageSize.z*imageSize.x*imageSize.y);
  std=sqrt(std);

  for (k=0;i<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
      for (j=0;j<imageSize.y;j++)
	outputImage[k][i][j]=ConvertDoubleToChar(inputImage[k][i][j],mean,3*std);
}

void Convert2DDoubleImageToCharImage(double ** inputImage,unsigned char ** outputImage,Ivector2d imageSize)
{
  int i,j,k;
  double mean,std;

  mean=0;

  for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
      mean+=inputImage[i][j]/(double)(imageSize.x*imageSize.y);
  
  std=0;
  for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
      std+=((inputImage[i][j]-mean)*(inputImage[i][j]-mean))/(double)(imageSize.x*imageSize.y);
  std=sqrt(std);

  for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
      outputImage[i][j]=ConvertDoubleToChar(inputImage[i][j],mean,3*std);
}

#define VALID_THRESHOLD 100

int CalculateGLCMGivenLabel(unsigned char*** origImage,unsigned char*** maskImage,int x_size,int y_size,int z_size,int label,float** GLCM,int GLCMSize,int interPixelDist,char direction)
{
  int i,j,k;
  int countedVolume;
  float** temp;

  countedVolume=0;

  for (i=0;i<GLCMSize;i++)
    for (j=0;j<GLCMSize;j++)
      GLCM[i][j]=0;

  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      for (j=0;j<y_size;j++)
	if (label==maskImage[k][i][j])
	  {
	    //printf("(%d %d %d)\n",i,j,k);
	    switch (direction)
	      {
	      case 'x':
		if (i+interPixelDist>=0&&i+interPixelDist<x_size)
		  if (label==maskImage[k][i+interPixelDist][j])
		    {
		      GLCM[origImage[k][i][j]][origImage[k][i+interPixelDist][j]]++;
		      countedVolume++;
		    }
		break;

	      case 'y':
		if (j+interPixelDist>=0&&j+interPixelDist<y_size)
		  if (label==maskImage[k][i][j+interPixelDist])
		    {
		      GLCM[origImage[k][i][j]][origImage[k][i][j+interPixelDist]]++;
		      countedVolume++;
		    }
		break;

	      case 'z':
		if (k+interPixelDist>=0&&k+interPixelDist<z_size)
		  if (label==maskImage[k+interPixelDist][i][j])
		    {
		      GLCM[origImage[k][i][j]][origImage[k+interPixelDist][i][j]]++;
		      countedVolume++;
		    }
		break;

	      default:
		break;
	      }
	  }
  
  if (countedVolume<VALID_THRESHOLD)
    return NNO;
  else
    {
      temp=Falloc2d(GLCMSize,GLCMSize);

      for (i=0;i<GLCMSize;i++)
	for (j=0;j<GLCMSize;j++)
	  temp[i][j]=GLCM[j][i];
  
      for (i=0;i<GLCMSize;i++)
	for (j=0;j<GLCMSize;j++)
	  GLCM[i][j]=GLCM[i][j]+temp[i][j];
      Ffree2d(temp,GLCMSize);


      for (i=0;i<GLCMSize;i++)
	for (j=0;j<GLCMSize;j++)
	  GLCM[i][j]/=(float)countedVolume*2.;
      
      return YYES;
    }
}

void DrawCircle2D(unsigned char** outputImage,Ivector2d imageSize,Ivector2d circleCenter,int circleSize,int label)
{
  int row,col;
  double r;
  int x,y;

  for (x=-circleSize;x<=circleSize;x++)
    for (y=-circleSize;y<=circleSize;y++)
      {
	r=sqrt((double)x*(double)x+(double)y*(double)y);
	
	if (r<=circleSize&&r>=circleSize-1)
	  {
	    row=circleCenter.x+x;
	    col=circleCenter.y+y;
	    
	    if (row>=0&&row<imageSize.x&&col>=0&&col<imageSize.y)
	      outputImage[row][col]=label;
	  }
      }
}

void DrawCross2D(unsigned char** outputImage,Ivector2d imageSize,Ivector2d crossCenter,int crossSize,int label,int thickness)
{
  int row,col;
  int r,t;
  
  for (r=-crossSize;r<=crossSize;r++)
    for (t=-thickness;t<=thickness;t++)
      {
	row=crossCenter.x+r;
	col=crossCenter.y+t;
	
	if (row>=0&&row<imageSize.x&&col>=0&&col<imageSize.y)	
	  {
	    if (label==-1)
	      if (outputImage[row][col]>128)
		label=0;
	      else
		label=255;
	    
	    outputImage[row][col]=label;      
	  }
      }

      
  for (r=-crossSize;r<=crossSize;r++)
    for (t=-thickness;t<=thickness;t++)
      {
	row=crossCenter.x+t;
	col=crossCenter.y+r;
	
	if (row>=0&&row<imageSize.x&&col>=0&&col<imageSize.y)
	  {
	    if (label==-1)
	      if (outputImage[row][col]>128)
		label=0;
	      else
		label=255;
	    
	    outputImage[row][col]=label;
	  }
      }

}

void DrawDirectionalCross2D(unsigned char** outputImage,Ivector2d imageSize,Ivector2d crossCenter,double crossOrientation,int crossSize,int label)
{
  int row,col;
  int r;

  int row1,col1;
  
  for (r=-2*crossSize;r<=2*crossSize;r++)
    {
      row=crossCenter.x+r*cos(crossOrientation);
      col=crossCenter.y+r*sin(crossOrientation);

      if (row>=0&&row<imageSize.x&&col>=0&&col<imageSize.y)
	outputImage[row][col]=label;

      /*
      row1=(cos(PIE/4.)*(double)row-sin(PIE/4.)*(double)col);
      col1=(sin(PIE/4.)*(double)row+cos(PIE/4.)*(double)col);

      if (row1>=0&&row1<imageSize.x&&col1>=0&&col1<imageSize.y)
	outputImage[row1][col1]=label;
      */
    }
  

  crossCenter.x-=crossSize*cos(crossOrientation);
  crossCenter.y-=crossSize*sin(crossOrientation);

  for (r=-crossSize;r<=crossSize;r++)
    {
      row=crossCenter.x+r*cos(crossOrientation+PIE/2.);
      col=crossCenter.y+r*sin(crossOrientation+PIE/2.);

      if (row>=0&&row<imageSize.x&&col>=0&&col<imageSize.y)
	outputImage[row][col]=label;

      /*
      row1=(cos(PIE/4.)*(double)row-sin(PIE/4.)*(double)col);
      col1=(sin(PIE/4.)*(double)row+cos(PIE/4.)*(double)col);

      if (row1>=0&&row1<imageSize.x&&col1>=0&&col1<imageSize.y)
	outputImage[row1][col1]=label;
      */
    }
}


void EdgeMap2D(double** inputImage,double** outputImage,Ivector2d imageSize,double smoothSigma)
{
  int i,j;
  double** smoothImage,**lx,**ly;
  double** derivativeMask;
  Ivector2d maskSize;
  

  smoothImage=Dalloc2d(imageSize.x,imageSize.y);
  Smoothing2DImageWithGaussianFilter(inputImage,smoothImage,imageSize,smoothSigma);

  lx=Dalloc2d(imageSize.x,imageSize.y);
  ly=Dalloc2d(imageSize.x,imageSize.y);

  // Calculate derivative along x direction
  derivativeMask=Dalloc2d(3,1);
  maskSize.x=3;
  maskSize.y=1;

  derivativeMask[0][0]=-0.5;
  derivativeMask[1][0]=0;
  derivativeMask[2][0]=0.5;
  
  Convolution2DGeneral(inputImage,imageSize,derivativeMask,maskSize,lx);

  Dfree2d(derivativeMask,3);

  // Calculate derivative along y direction
  derivativeMask=Dalloc2d(1,3);
  maskSize.x=1;
  maskSize.y=3;

  derivativeMask[0][0]=-0.5;
  derivativeMask[0][1]=0;
  derivativeMask[0][2]=0.5;
  
  Convolution2DGeneral(inputImage,imageSize,derivativeMask,maskSize,ly);
  Dfree2d(derivativeMask,1);

  for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
      outputImage[i][j]=sqrt(lx[i][j]*lx[i][j]+ly[i][j]*ly[i][j]);
  
  Dfree2d(smoothImage,imageSize.x);
  Dfree2d(lx,imageSize.x);
  Dfree2d(ly,imageSize.y);
}
