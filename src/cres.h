#ifndef CRES_H
#define CRES_H

#define TINY 1.0e-20;

#ifdef __cplusplus
 extern "C" {
 #endif


typedef struct
{
  float H;
  float K;
  float k1;
  float k2;
  float error;
} Skappa;

void nrerror(char []);
double *vector(int,int);
int *ivector(int,int);
double *devecotr(int,int);
double **matrix(int,int,int,int);
double **dmatrix(int,int,int,int);
int **imatrix(int,int,int,int);
double **submatrix(double **,int,int,int,int,int,int);
void free_vector(double *,int,int);
void free_ivector(int *,int,int);
void free_dvector(double *,int,int);
void free_matrix(double **,int,int,int,int);
void free_dmatrix(double **,int,int,int,int);
void free_imatrix(int **,int,int,int,int);
void ludcmp(double **,int,int *,double *);
void lubksb(double **,int,int *,double *);
void matrixInv(double **,double **,int);
void matrixDet(double *,double **,int);
Skappa estimate_kappa(Fvector3d *,int,Fvector3d);
void Fvector3dcurvature2d(int,int,Fvector3d **,float **,int,float);


 #ifdef __cplusplus
 }
 #endif


#endif
