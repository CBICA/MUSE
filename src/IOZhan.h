#ifndef IOZHAN_H
#define IOZHAN_H

#include "mvcd.h"

#define BIN 0
#define ASCII 1

#ifdef __cplusplus
 extern "C" {
 #endif


//Function Prototype
void Read3DImage(char* filename,unsigned char*** Image,int x_size,int y_size,int z_size);
void Write3DImage(char* filename,unsigned char*** Image,int x_size,int y_size,int z_size);

void ReadTrainingDataPara(char* filename,int* sampleNum,int* dimension,int format);
void ReadTrainingData(char* filename,int sampleNum,int dimension,float** sample,int format);
void WriteTrainingData(char* filename,int sampleNum,int dimension,float** sample,int format);

void SaveDeformationField(char* filename,Fvector3d ***deformationField,int image_X,int image_Y,int image_Z);
void LoadDeformationField(char* filename,Fvector3d ***deformationField,int image_X,int image_Y,int image_Z);

void Read3DFloatImage(char* filename,float*** Image,int x_size,int y_size,int z_size);
void Write3DFloatImage(char* filename,float*** Image,int x_size,int y_size,int z_size);

void Read3DShortImage(char* filename,short*** Image,int x_size,int y_size,int z_size);
void Write3DShortImage(char* filename,short*** Image,int x_size,int y_size,int z_size);

void Read2DFloatImage(char* filename,float** Image,int x_size,int y_size);
void Write2DFloatImage(char* filename,float** Image,int x_size,int y_size);

void ReadVertexNum(char* filename,int* vertexNum);
void ReadVertex(char* filename,Fvector3d* vertex,int vertexNum);
void WriteVertex(char* filename,Fvector3d* vertex,int vertexNum, int asciiOutput);



#ifdef __cplusplus
 }
 #endif

#endif
