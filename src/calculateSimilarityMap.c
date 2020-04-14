#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <general.h>
#include <mvcd.h>
#include <matrixZHAN.h>  
#include <cres.h>
#include <IOZhan.h>
#include <imageZHAN.h>
#include <watershed.h>
#include <unistd.h>

#define YYES 1
#define NNO 0
#define pi 3.1415926

char inputImgNameA[200], inputImgNameB[200], outputImgName[200];
char featureListNameA[200], featureListNameB[200], featureListPrefix[200];
char fileName[200];

int singlePointOrNot;
Ivector3d point;

int i,j,k;
int ii,jj,kk;
int s,t;
int x,y,z;
int featureIndex;
Ivector3d imageSize;
int numFeatures, numFeaturesA, numFeaturesB;

unsigned char ***inputImgA, ***inputImgB;
unsigned char ***mask;
float ***outputImg;
float ***similarityMap;
unsigned char ****featureMapA, ****featureMapB;
//float ****featureMapAFloat, ****featureMapBFloat;
float ****featureMapANormalized, ****featureMapBNormalized;



void ReadFeaturesFromFeatureList(char* featureImageListFile,unsigned char ****featureMap, int x_size,int y_size, int z_size);
//void calculateSaliencyMap(float ****featureMap, unsigned char ***mask, int numFeatures, Ivector3d imageSize, float ratio, int radiusXY, int radiusZ, float ***salienceMap);
float calculateEuclideanDistanceBetweenTwoFeatureVectors(unsigned char ****featureMapA, unsigned char ****featureMapB, int x, int y, int z, int numFeatures);
void generateMask(unsigned char ***inputImgA, unsigned char ***inputImgB, Ivector3d imageSize, unsigned char ***mask);
void featureVectorNormalization(unsigned char ****featureMap, float ****featureMapNormalized, int numFeatures, Ivector3d imageSize);
void calculateSimilarityMap(unsigned char ****featureMapA, unsigned char ****featureMapB, unsigned char ***mask, int numFeatures, Ivector3d imageSize, int similaritydefinition, float ***similarityMap);
void convertUCFeatureMaptoFloat(unsigned char ****featureMap, float ****featureMapFloat, Ivector3d imageSize, int numFeatures);
void calculateSimilarityMapForAPoint(unsigned char ****featureMapA, unsigned char ****featureMapB, unsigned char ***mask, int numFeatures, Ivector3d imageSize, Ivector3d point, int similaritydefinition, float ***similarityMap);
void DrawCross2DFloat(float** outputImage,Ivector2d imageSize,Ivector2d crossCenter,int crossSize,float label,int thickness);

void show_usage()
{
printf("\n USAGE:  calculateSimilarityMap <inImgA><inImgB><featureListPrefix><outImg> \n\
\t This program calculates voxel-wise similarity between two images based on their feature profiles, or between a point in image A to all points in image B.\n\
\t -d  <int>,<int>         : image dimension (default 256,256)\n\
\t -p  <int>,<int>,<int> : calculate the similarity between a point in A and all points in B, the input is the x,y,z coordinates of point p\n\
\t -w  <int>               : similarity definition (default: 1 -- sim=1/(1+dist);  2 -- sim=exp(-dist/50) ) \n\
  example1: calculateSimilarityMap A.img B.img Gabor_ similarityMap.img -d256,256 -w2\n\
  example2: calculateSimilarityMap A.img B.img Gabor_ similarityMap.img -d256,256 -p100,100,1\n\n\
  note: the point under consideration is from image A\n\
");

  exit(1);
}


int main(int argc,char *argv[])
{
  int c, num;
  FILE *fpA, *fpB;
  float dist;
  
  
  num=4;
  if(argc<num)
    show_usage();

  imageSize.x = 256;
  imageSize.y = 256;
  singlePointOrNot = NNO;
  int similaritydefinition = 1;

    
  while((c=getopt(argc-num,argv+num,"d:p:w:")) != -1)
    {
      switch(c)
		{
			case 'd':
				sscanf(optarg,"%d,%d",&imageSize.x,&imageSize.y);
				break;
				
			case 'p':
				sscanf(optarg,"%d,%d,%d",&point.x,&point.y,&point.z);
				singlePointOrNot = YYES;
				break;

			case 'w':
				sscanf(optarg, "%d", &similaritydefinition);
				break;
				
			default:
				break;
		}
    }

   sprintf(inputImgNameA, "%s", argv[1]);
   sprintf(inputImgNameB, "%s", argv[2]);
   sprintf(featureListPrefix, "%s", argv[3]);
   sprintf(outputImgName, "%s", argv[4]);
   
   printf("\ninput image A= %s, input image B = %s\nfeatureListPrefix = %s\n", inputImgNameA, inputImgNameB, featureListPrefix);
   printf("output image = %s\n", outputImgName);

   DetectUCImageSize(inputImgNameA, imageSize.x, imageSize.y, &imageSize.z);
   printf("\nimage size = (%d, %d, %d)\n", imageSize.x, imageSize.y, imageSize.z);
   printf("\nsimilarity definition = %d\n", similaritydefinition);
   
   inputImgA = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);
   inputImgB = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);
   outputImg = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
   mask = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);
   similarityMap = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
   printf("Read input images...\n");
   Read3DImage(inputImgNameA, inputImgA, imageSize.x, imageSize.y, imageSize.z);
   Read3DImage(inputImgNameB, inputImgB, imageSize.x, imageSize.y, imageSize.z);
   
   // generate a mask
   printf("generate mask...");
   generateMask(inputImgA, inputImgB, imageSize, mask);
   printf("done!\n");
   
   // read features for both image under registration
   sprintf(featureListNameA, "%sA.lst", featureListPrefix);
   if (NULL==(fpA=fopen(featureListNameA,"rb"))){
       printf("File %s doesn't exist!\n",featureListNameA);
	   exit(1);
	   }
   fscanf(fpA,"%d",&numFeaturesA);
   
   sprintf(featureListNameB, "%sB.lst", featureListPrefix);
   if (NULL==(fpB=fopen(featureListNameB,"rb"))){
       printf("File %s doesn't exist!\n",featureListNameB);
	   exit(1);
	   }
   fscanf(fpB,"%d",&numFeaturesB);
   
   if (numFeaturesA!=numFeaturesB)
     {
	 printf("Error: the numbers of features must match!\n");
	 exit(1);
	 }
   else
     numFeatures = numFeaturesA;
   printf("total number of features = %d\n", numFeatures);

   // allocate space for feature map
   printf("Allocate memory for feature images...\n");
   featureMapA = (unsigned char****)malloc(sizeof(unsigned char***)*numFeatures);	   
   featureMapB = (unsigned char****)malloc(sizeof(unsigned char***)*numFeatures);
   //featureMapANormalized = (float****)malloc(sizeof(float ***)*numFeatures);
   //featureMapBNormalized = (float****)malloc(sizeof(float ***)*numFeatures);
   //featureMapAFloat = (float****)malloc(sizeof(float ***)*numFeatures);
   //featureMapBFloat = (float****)malloc(sizeof(float ***)*numFeatures);
   
   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
     {
	 featureMapA[featureIndex] = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);	 
	 featureMapB[featureIndex] = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);	 
	 //featureMapAFloat[featureIndex] = Falloc3d(imageSize.x, imageSize.y, imageSize.z);	 
	 //featureMapBFloat[featureIndex] = Falloc3d(imageSize.x, imageSize.y, imageSize.z);	 
	 //featureMapANormalized[featureIndex] = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
	 //featureMapBNormalized[featureIndex] = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
	 }
   
   // read features 
   printf("\nReading feature images A...\n");
   ReadFeaturesFromFeatureList(featureListNameA, featureMapA, imageSize.x, imageSize.y, imageSize.z);
   printf("\nReading feature images B...\n");
   ReadFeaturesFromFeatureList(featureListNameB, featureMapB, imageSize.x, imageSize.y, imageSize.z);
	
   //normalize feature vector
   //featureVectorNormalization(featureMapA, featureMapANormalized, numFeatures, imageSize);
   //featureVectorNormalization(featureMapB, featureMapBNormalized, numFeatures, imageSize);
   //convertUCFeatureMaptoFloat(featureMapA, featureMapAFloat, imageSize, numFeatures);
   //convertUCFeatureMaptoFloat(featureMapB, featureMapBFloat, imageSize, numFeatures);
   
   if (singlePointOrNot==NNO)   // similarity between two images
	{
	printf("\ncalculate similarity map between two images...\n");
    //calculateSimilarityMap(featureMapANormalized, featureMapBNormalized, mask, numFeatures, imageSize, similarityMap);
    calculateSimilarityMap(featureMapA, featureMapB, mask, numFeatures, imageSize, similaritydefinition, similarityMap);
	}
   else
    {
	printf("\ncalculate similarity between point (%d, %d, %d) in image A and all points in image B...\n", point.x, point.y, point.z);
	calculateSimilarityMapForAPoint(featureMapA, featureMapB, mask, numFeatures, imageSize, point, similaritydefinition, similarityMap);
	}
   
   // save similarity map
   printf("\nSaving similarity map into %s\n\n", outputImgName);
   Write3DFloatImage(outputImgName, similarityMap, imageSize.x, imageSize.y, imageSize.z);
   
   UCfree3d(inputImgA, imageSize.z, imageSize.x);
   UCfree3d(inputImgB, imageSize.z, imageSize.x);
   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
     {
	 UCfree3d(featureMapA[featureIndex], imageSize.z, imageSize.x);
	 UCfree3d(featureMapB[featureIndex], imageSize.z, imageSize.x);
	 //Ffree3d(featureMapAFloat[featureIndex], imageSize.z, imageSize.x);
	 //Ffree3d(featureMapBFloat[featureIndex], imageSize.z, imageSize.x);
	 }
   free(featureMapA);
   free(featureMapB);
   //free(featureMapAFloat);
   //free(featureMapBFloat);
}



void ReadFeaturesFromFeatureList(char* featureImageListFile,unsigned char ****featureMap, int x_size,int y_size, int z_size)
{
  unsigned char ***featureImage;
  FILE *fp;
  int featureIndex, featureNum;
  int i,j,k;
  char filename[200];

  fp=fopen(featureImageListFile,"rb");
  fscanf(fp,"%d",&featureNum);
  
  featureImage=UCalloc3d(x_size,y_size, z_size);
    
  for (featureIndex=0;featureIndex<featureNum;featureIndex++)
    {
      fscanf(fp,"%s",filename);

      printf("feature image #%d out of %d: %s\n", (featureIndex+1), featureNum, filename);
      Read3DImage(filename, featureImage, x_size, y_size, z_size);

	  for (k=0;k<z_size;k++)
       for (i=0;i<x_size;i++)
		for (j=0;j<y_size;j++)
		  featureMap[featureIndex][k][i][j]=featureImage[k][i][j];	      
    }//for (featureIndex=0;featureIndex<featureNum;featureIndex++)
  
  //UCfree3d(featureImage,z_size,x_size);
  UCfree3d(featureImage,z_size, x_size);
  fclose(fp);
}


/*
void calculateSaliencyMap(float ****featureMap, unsigned char ***mask, int numFeatures, Ivector3d imageSize, float ratio, int radiusXY, int radiusZ, float ***salienceMap)
{
  int i,j,k,ii,jj,kk;
  int x,y,z;
  float maxDist, dist;
  float meanDisSimilarity, stdDisSimilarity, disSimilarity;
  int numVoxelsInNeighborhood;
  float salienceMin=1.0;
  float salienceMax=0.0;
  
 
  for (k=0;k<imageSize.z;k++)
     {
	 printf("\nslice #%d\n",k);
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 salienceMap[k][i][j]=0.0;
		 
		 if (mask[k][i][j]>0)
		   {
			salienceMap[k][i][j] = 0.00; // initilize
			numVoxelsInNeighborhood = 0;
			maxDist = 0.0;
			meanDisSimilarity = 0.0;
		 
			// calculate mean dis-similarity in the neighborhood of (i,j,k)
			for (x=-radiusXY;x<=radiusXY;x++)
			  for (y=-radiusXY;y<=radiusXY;y++)
				for (z=-radiusZ;z<=radiusZ;z++)
				  {
					ii = i+x;
					jj = j+y;
					kk = k+z;
			   
					dist = 0.0;
					if (ii>=0 && ii<imageSize.x && jj>=0 && jj<imageSize.y && kk>=0 && kk<imageSize.z && !(ii==i && jj==j && kk==k))
					 {
						numVoxelsInNeighborhood++;
						for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
							dist += pow((featureMap[featureIndex][k][i][j]-featureMap[featureIndex][kk][ii][jj]), 2.0);
						dist = sqrt(dist);
						meanDisSimilarity += dist;
				 
						if (dist>maxDist)
							maxDist = dist;
					 } //if 
				  } // for x, for y, for z
		 
			meanDisSimilarity /= (maxDist*(float)numVoxelsInNeighborhood);
		 
			// calculate std of dis-similarities in the neighborhood of (i,j,k)
			stdDisSimilarity = 0.0;
			for (x=-radiusXY;x<=radiusXY;x++)
			  for (y=-radiusXY;y<=radiusXY;y++)
				for (z=-radiusZ;z<=radiusZ;z++)
				  {
					ii = i+x;
					jj = j+y;
					kk = k+z;
			   
					dist = 0.0;
					if (ii>=0 && ii<imageSize.x && jj>=0 && jj<imageSize.y && kk>=0 && kk<imageSize.z && !(ii==i && jj==j && kk==k))
					  {
						for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
							dist += pow((featureMap[featureIndex][k][i][j]-featureMap[featureIndex][kk][ii][jj]), 2.0);
						dist = sqrt(dist);
						disSimilarity = dist/maxDist;
				 
						stdDisSimilarity += pow((disSimilarity-meanDisSimilarity), 2.0);
					  } //if 
				  } // for x, for y, for z
			stdDisSimilarity = sqrt(stdDisSimilarity/(float)numVoxelsInNeighborhood);
		 
			// combine mean and std dis-similarity into distinctiveness
			salienceMap[k][i][j] = meanDisSimilarity - ratio*stdDisSimilarity;
		 
			if (salienceMap[k][i][j]>salienceMax)
			  salienceMax = salienceMap[k][i][j];
			if (salienceMap[k][i][j]<salienceMin)
			  salienceMin = salienceMap[k][i][j];
		   
		   } // if mask>0
			//printf("At (%d, %d, %d), mean and std dis-similarity = %f and %f, distinctiveness = %f\n", i,j,k,meanDisSimilarity, stdDisSimilarity, outputImg[k][i][j]);
		} // for i,j
     } // for k
   
   
   // linearly normalize the salience value at each voxel into (0,1)
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 if (mask[k][i][j]==0)
		   salienceMap[k][i][j] = 0.0;
		 else
	 	   salienceMap[k][i][j] = (salienceMap[k][i][j]-salienceMin)/(salienceMax-salienceMin);
		 }
	
   //Write3DFloatImage(outputImgName, outputImg, imageSize.x, imageSize.y, imageSize.z);
}

*/

float calculateEuclideanDistanceBetweenTwoFeatureVectors(unsigned char ****featureMapA, unsigned char ****featureMapB, int x, int y, int z, int numFeatures)
{
   float dist=0.0;
   int featureIndex;
   
   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
     {
     if ( (featureMapB[featureIndex][z][x][y]!=0)&&(featureMapA[featureIndex][z][x][y]!=0) )
      dist += pow((featureMapA[featureIndex][z][x][y]-featureMapB[featureIndex][z][x][y]), 2.0);
	 }
   
   dist = sqrt(dist/(float)numFeatures);
   
   return dist;
}

void calculateSimilarityMap(unsigned char ****featureMapA, unsigned char ****featureMapB, unsigned char ***mask, int numFeatures, Ivector3d imageSize, int similaritydefinition, float ***similarityMap)
{
   int i,j,k;
   float dist;
   float maxSim=0;
   
   for (k=0;k<imageSize.z;k++)
     {
	 printf("slice #%d\n",k);
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 if (mask[k][i][j]>0)
		   {
		   dist = calculateEuclideanDistanceBetweenTwoFeatureVectors(featureMapA, featureMapB, i, j, k, numFeatures);
		   
		   if (similaritydefinition==1) {
			   similarityMap[k][i][j] = 1.0/(dist+1.0);
			   //similarityMap[k][i][j] = 20/(dist+20);
			   //printf("at (%d, %d, %d), dist = %f, sim = %f\n", i,j,k,dist, similarityMap[k][i][j]);
			}
		   else
			similarityMap[k][i][j] = exp(-dist/50.0);
		   }
		 else
		   similarityMap[k][i][j] = 0;
		   
		 if (similarityMap[k][i][j]>maxSim)
		   maxSim = similarityMap[k][i][j];
		 } // for...for
	 } //for
	
  /*	
  // normalize the similarity into [0,1] by dividing the maximum similarity
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
	  for (j=0;j<imageSize.y;j++)
	    similarityMap[k][i][j] = similarityMap[k][i][j]/maxSim;
     */
}


void generateMask(unsigned char ***inputImgA, unsigned char ***inputImgB, Ivector3d imageSize, unsigned char ***mask)
{
   int i,j,k;
   
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 //if ( inputImgA[k][i][j]>0 || inputImgB[k][i][j]>0 )
		 if (inputImgB[k][i][j]>0)
		    mask[k][i][j] = 255;
	     else
		    mask[k][i][j] = 0;
		 }
}



void featureVectorNormalization(unsigned char ****featureMap, float ****featureMapNormalized, int numFeatures, Ivector3d imageSize)
{
   int featureIndex, i,j,k;
   double mag;
   double maxMag=0.0;
   
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 mag = 0;
		 for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
		   mag += pow((double)featureMap[featureIndex][k][i][j], 2.0);
		 
		 mag = sqrt(mag);
		 
		 if (mag>maxMag)
		   maxMag = mag;
		 }
		 
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
		   {
		   if (mag==0)
		     featureMapNormalized[featureIndex][k][i][j] = 0.0;
		   else
		     featureMapNormalized[featureIndex][k][i][j] = (double)featureMap[featureIndex][k][i][j]/maxMag;
		   }
		 }
}



void convertUCFeatureMaptoFloat(unsigned char ****featureMap, float ****featureMapFloat, Ivector3d imageSize, int numFeatures)
{
   int i,j,k, featureIndex;
   
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
		   featureMapFloat[featureIndex][k][i][j] = (float)featureMap[featureIndex][k][i][j];
		 }
}



void calculateSimilarityMapForAPoint(unsigned char ****featureMapA, unsigned char ****featureMapB, unsigned char ***mask, int numFeatures, Ivector3d imageSize, Ivector3d point, int similaritydefinition, float ***similarityMap)
{
   int i,j,k, featureIndex;
   float dist;
   float maxSim=0;
   Ivector2d imageSize2D, crossCenter;
   imageSize2D.x = imageSize.x;
   imageSize2D.y = imageSize.y;
   crossCenter.x = point.x;
   crossCenter.y = point.y;
   
   
   // for calculate similarity on feature vectors (begin)
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 dist = 0.0;
		 //if ( (mask[k][i][j]!=0)&(sqrt(pow((i-point.x),2.0)+pow((j-point.y),2.0))<=25) )
		 if (mask[k][i][j]==0)
		   similarityMap[k][i][j]=0.0;
		 else
		   {
		   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
			 dist += pow((featureMapA[featureIndex][point.z][point.x][point.y]-featureMapB[featureIndex][k][i][j]), 2.0);
		   dist = sqrt(dist/(float)numFeatures);

		   if (similaritydefinition==1)
			   similarityMap[k][i][j] = 20.0/(20.0+dist);
		   else
			similarityMap[k][i][j]=exp(-dist/75.0);
		   }
        
		   
		 if (similarityMap[k][i][j]>maxSim)
		   maxSim = similarityMap[k][i][j];
		 }
		 printf("max similarity = %f\n", maxSim);
    
   /*	
   // normalize similarity into [0,1] by dividing similarity at each pixel by the largest similarity 
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 similarityMap[k][i][j] = similarityMap[k][i][j]/maxSim;
		 }
	*/
	
   // for calculate similarity on feature vectors (end)
   
   
   /*
   // for calculate similarity on each individual feature (begin)
   char filename[200];
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 // get the maximum similarity
		 //if ( (mask[k][i][j]!=0)&(sqrt(pow((i-point.x),2.0)+pow((j-point.y),2.0))<=10) )
		 if (mask[k][i][j]!=0)
		   {
		 
		   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
		   {
		   dist = fabs((float)featureMapA[featureIndex][point.z][point.x][point.y]-(float)featureMapB[featureIndex][k][i][j]);
		   
		   if (dist!=0)
		     similarityMap[k][i][j] = 1/dist;
		   else
		     similarityMap[k][i][j] = 0;
		   
		   if (similarityMap[k][i][j]>maxSim)
		     maxSim = similarityMap[k][i][j];
		   }
		   }
		 }
   printf("max sim=%f\n", maxSim);
   
   // normalize similarity into [0,1] by dividing similarity at each pixel by the largest similarity 
   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
   {
   for (k=0;k<imageSize.z;k++)
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 //if ( (mask[k][i][j]!=0)&(sqrt(pow((i-point.x),2.0)+pow((j-point.y),2.0))<=10) )
		 if (mask[k][i][j]!=0)
		   {
		   dist = fabs((float)featureMapA[featureIndex][point.z][point.x][point.y]-(float)featureMapB[featureIndex][k][i][j]);
		   
		   if (dist!=0)
		     similarityMap[k][i][j]=(1/dist)/maxSim;
		   else
		     similarityMap[k][i][j]=0.0;
		   }
		 else
		   similarityMap[k][i][j]=0;
		 }
   sprintf(filename, "sim_x%dy%d_Gabor%d.img", point.x, point.y, (featureIndex+1));
   printf("save %s\n", filename);
   Write3DFloatImage(filename, similarityMap, imageSize.x, imageSize.y, imageSize.z);
   } // for featureIndex
   // for calculate similarity on each individual feature (end)
     */
   
   
   //printf("test1\n");	 
   //DrawCross2DFloat(similarityMap[point.z],imageSize2D,crossCenter,5,0.5,1);
   //printf("test2\n");
}



void DrawCross2DFloat(float** outputImage,Ivector2d imageSize,Ivector2d crossCenter,int crossSize,float label,int thickness)
{
  int row,col;
  int r,t;
  
  //printf("crossCenter = (%d,%d), cross size = %d, lable = %f, thickness = %d\n", crossCenter.x, crossCenter.y, crossSize, label, thickness);
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



