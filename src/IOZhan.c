#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mvcd.h"
#include "IOZhan.h"

void Read3DImage(char* filename,unsigned char*** Image,int x_size,int y_size,int z_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
      
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      fread(Image[k][i],1,y_size,fp);
  fclose(fp);
}


void Write3DImage(char* filename,unsigned char*** Image,int x_size,int y_size,int z_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"w")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
  
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      fwrite(Image[k][i],1,y_size,fp);  
  fclose(fp);
}

void ReadTrainingDataPara(char* filename,int* sampleNum,int* dimension,int format)
{
  FILE* fp;
  int s,n;

  printf("\nReading Data Parameters...\n");
  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }

  if (format==BIN)
    {
      fread(&s,sizeof(int),1,fp);
      fread(&n,sizeof(int),1,fp);
      
    }
  else
    {
      fscanf(fp, "%d", &s);
      fscanf(fp, "%d", &n);
    }


  if (format==BIN)
    *dimension=n;
  else
    *dimension=n-1;

  *sampleNum=s;

  fclose(fp);
}

void ReadTrainingData(char* filename,int sampleNum,int dimension,float** sample,int format)
{
  FILE* fp;
  int s,n,i,j;
  float* symbol;

  printf("\nReading Data...\n");

  symbol=Falloc1d(sampleNum);
  
  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
      
  if (format==BIN)
    {
      fseek(fp,sizeof(int)*2,SEEK_SET);
      
      for (i=0;i<sampleNum;i++)
	fread(sample[i],sizeof(float),dimension,fp);

      fread(symbol,sizeof(float),sampleNum,fp);

      for (i=0;i<sampleNum;i++)
	sample[i][dimension]=symbol[i];

    }
  else
    {
      fscanf(fp, "%d", &s);
      fscanf(fp, "%d", &n);
      for (i=0;i<sampleNum;i++)
        for (j=0;j<dimension+1;j++)          
          fscanf(fp, "%f", &sample[i][j]);            
    }
      
  free(symbol);
  fclose(fp);
  printf("Reading data finished!\n");
}

void WriteTrainingData(char* filename,int sampleNum,int dimension,float** sample,int format)
{
  FILE* fp;
  int s,n,i,j;
  float* symbol;

  printf("\nWriting Data...\n");

  symbol=Falloc1d(sampleNum);

  for (i=0;i<sampleNum;i++)
    symbol[i]=sample[i][dimension];

  if (NULL==(fp=fopen(filename,"w")))
      {
	printf("Unable to open %s!\n",filename);
      }

  if (format==BIN)
    {
      fwrite(&sampleNum,sizeof(int),1,fp);
      fwrite(&dimension,sizeof(int),1,fp);
      
      for (i=0;i<sampleNum;i++)
	fwrite(sample[i],sizeof(float),dimension,fp);

      fwrite(symbol,sizeof(float),sampleNum,fp);
      
    }
  else
    {
      fprintf(fp,"%d %d\n",sampleNum,dimension+1);
      
      for (i=0;i<sampleNum;i++)
	{
	  for (j=0;j<dimension+1;j++)
	    fprintf(fp,"%f ",sample[i][j]);
	  
	  fprintf(fp,"\n");
	}
    }
  
  fclose(fp);

  free(symbol);
}

void LoadDeformationField(char* filename,Fvector3d ***deformationField,int image_X,int image_Y,int image_Z)
{
  int i,j,k;
  FILE* fp;
  
  fp=myopen(filename,"r");

  for(k=0;k<image_Z;k++)
    for(i=0;i<image_X;i++)
      fread(deformationField[k][i],sizeof(Fvector3d),image_Y,fp);
  fclose(fp);
}

void SaveDeformationField(char* filename,Fvector3d ***deformationField,int image_X,int image_Y,int image_Z)
{
  int i,j,k;
  FILE* fp;
  
  fp=myopen(filename,"w");

  for(k=0;k<image_Z;k++)
    for(i=0;i<image_X;i++)
      fwrite(deformationField[k][i],sizeof(Fvector3d),image_Y,fp);
  fclose(fp);
}

void Read3DFloatImage(char* filename,float*** Image,int x_size,int y_size,int z_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
      
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      fread(Image[k][i],sizeof(float),y_size,fp);
  fclose(fp);
}

void Write3DFloatImage(char* filename,float*** Image,int x_size,int y_size,int z_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"w")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
  
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      fwrite(Image[k][i],sizeof(float),y_size,fp);  
  fclose(fp);
}

void Read2DFloatImage(char* filename,float** Image,int x_size,int y_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
      
  for (i=0;i<x_size;i++)
      fread(Image[i],sizeof(float),y_size,fp);
  fclose(fp);
}

void Write2DFloatImage(char* filename,float** Image,int x_size,int y_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"w")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
      
  for (i=0;i<x_size;i++)
      fwrite(Image[i],sizeof(float),y_size,fp);
  fclose(fp);
}

void Read3DShortImage(char* filename,short*** Image,int x_size,int y_size,int z_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
      
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      fread(Image[k][i],sizeof(short),y_size,fp);
  fclose(fp);
}

void Write3DShortImage(char* filename,short*** Image,int x_size,int y_size,int z_size)
{
  int i,j,k;
  FILE * fp;

  if (NULL==(fp=fopen(filename,"w")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
  
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      fwrite(Image[k][i],sizeof(short),y_size,fp);  
  fclose(fp);
}

void ReadVertexNum(char* filename,int* vertexNum)
{
  FILE * fp;
  
  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
  
  fscanf(fp,"%d",vertexNum);

  fclose(fp);
}

void ReadVertex(char* filename,Fvector3d* vertex,int vertexNum)
{
  FILE * fp;
  int i,temp;

  if (NULL==(fp=fopen(filename,"r")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
  
  fscanf(fp,"%d",&temp);

  for (i=0;i<vertexNum;i++)
    {
      fscanf(fp,"%f %f %f",&vertex[i].x,&vertex[i].y,&vertex[i].z);
    }
  
  fclose(fp);
}

void WriteVertex(char* filename,Fvector3d* vertex,int vertexNum, int asciiOutput)
{
  FILE * fp;
  int i;

  if (NULL==(fp=fopen(filename,"w")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }


  if (!asciiOutput){
      fwrite(vertex, sizeof(Fvector3d), vertexNum, fp);
      printf("Writing vertex file in BINARY format...\n");
     }
  else
     {
       fprintf(fp,"%d\n",vertexNum);
       printf("Writing vertex file in ASCII format...\n");

       for (i=0;i<vertexNum;i++)
          {
            fprintf(fp,"%f %f %f\n",vertex[i].x,vertex[i].y,vertex[i].z);
          }
     } 
  fclose(fp);
}
