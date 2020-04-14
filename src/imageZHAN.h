#ifndef IMAGEZHAN_H
#define IMAGEZHAN_H

#ifdef __cplusplus
 extern "C" {
 #endif


void DetectUCImageSize(char* filename,int x_size,int y_size,int* z_size);
void Convolution1D(double* input,int inputSize,double* mask,int maskSize,double* output);
void Convolution2DSeparate(double** input,Ivector2d inputSize,double* mask,int maskSize,double** output);
void Convolution2DGeneral(double** input,Ivector2d inputSize,double** mask,Ivector2d maskSize,double** output);
void Smoothing2DImageWithGaussianFilter(double** inputImage,double** outputImage,Ivector2d imageSize,double smoothSigma);

unsigned char ConvertDoubleToChar(double input,double translation,double scale);
void Convert3DDoubleImageToCharImage(double *** inputImage,unsigned char *** outputImage,Ivector3d imageSize);
void Convert2DDoubleImageToCharImage(double ** inputImage,unsigned char ** outputImage,Ivector2d imageSize);

void EdgeMap2D(double** inputImage,double** outputImage,Ivector2d imageSize,double smoothSigma);

void DrawCross2D(unsigned char** outputImage,Ivector2d imageSize,Ivector2d crossCenter,int crossSize,int label,int thickness);
void DrawCircle2D(unsigned char** outputImage,Ivector2d imageSize,Ivector2d circleCenter,int circleSize,int label);
void DrawDirectionalCross2D(unsigned char** outputImage,Ivector2d imageSize,Ivector2d crossCenter,double crossOrientation,int crossSize,int label);

#ifdef __cplusplus
 }
 #endif


#endif
