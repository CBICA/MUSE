#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mvcd.h"
#include "cres.h"

void nrerror(char error_text[])
{
  
  printf("Numerical Recipes run-time error...\n");
  printf("%s\n",error_text);
  printf("...now exiting to system...\n");
  exit(1);
}

double *vector(int nl,int nh)
{
  double *v;
  
  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl;
}

int *ivector(int nl,int nh)
{
  int *v;
  
  v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl;
}

double *dvector(int nl,int nh)
{
  double *v;
  
  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl;
}

double **matrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  int i;
  double **m;
  
  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  int i;
  double **m;
  
  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

int **imatrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  int i,**m;
  
  m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
  if (!m) nrerror("allocation failure 1 in imatrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
    if (!m[i]) nrerror("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
  return m;
}

double **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
     double **a;
     int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
  int i,j;
  double **m;
  
  m=(double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(double*));
  if (!m) nrerror("allocation failure in submatrix()");
  m -= newrl;
  
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;
  
  return m;
}

void free_vector(double *v,int nl,int nh)
{
  free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
     int *v,nl,nh;
{
  free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
     double *v;
     int nl,nh;
{
  free((char*) (v+nl));
}

void free_matrix(m,nrl,nrh,ncl,nch)
     double **m;
     int nrl,nrh,ncl,nch;
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
     double **m;
     int nrl,nrh,ncl,nch;
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
     int **m;
     int nrl,nrh,ncl,nch;
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void ludcmp(double **a,int n,int *indx,double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;
  
  
  vv=vector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,1,n);
}

void lubksb(double **a,int n,int *indx,double *b)
{
  int i,ii=0,ip,j;
  double sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


Skappa estimate_kappa(nb,n,normal)
     Fvector3d *nb,normal; /* an array of relative position of neighbors w.r.t. central point */
     int n;                   /* number of neighbors */

{
  int i,indx[7],locN;
  double A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T;
  double a,b,c,d,f,g;
  Skappa curvature;
  double **m,vec[7],error;
  double xx,yy,xxx,yyy,x,y,z,denom,*dd,*px,*py,*pz;
  Fvector3d u,v,w;  /* The orthonormal basis based on the surface tangent ***/
  double temp,temp1;
  A=B=C=D=E=F=G=H=I=J=K=L=M=N=O=P=Q=R=S=T=0.0;

  px=(double *)calloc(n,sizeof(double));
  py=(double *)calloc(n,sizeof(double));
  pz=(double *)calloc(n,sizeof(double));
  m=Dalloc2d(7,7);

  for(i=0;i<7;i++)
    {
      indx[i]=0;
      vec[i]=0;
    }
  /*** Compute the new basis ***/

  temp=sqrt((double) (nb[1].x*nb[1].x+nb[1].y*nb[1].y+nb[1].z*nb[1].z));
  u.x=nb[1].x/temp;
  u.y=nb[1].y/temp;
  u.z=nb[1].z/temp;
  
  temp1=nb[2].x*u.x+nb[2].y*u.y+nb[2].z*u.z;
  v.x=nb[2].x-temp1*u.x;
  v.y=nb[2].y-temp1*u.y;
  v.z=nb[2].z-temp1*u.z;
  temp=sqrt((double) (v.x*v.x+v.y*v.y+v.z*v.z));
  v.x/=temp;
  v.y/=temp;
  v.z/=temp;

  
  w.x=u.y*v.z-v.y*u.z;
  w.y=v.x*u.z-u.x*v.z;
  w.z=u.x*v.y-v.x*u.y;

  for(i=0;i<n;i++)
    {
      /*** Find the coordinates wrt the new orthonormal coordinate system ***/
      x=nb[i].x*u.x+nb[i].y*u.y+nb[i].z*u.z;
      y=nb[i].x*v.x+nb[i].y*v.y+nb[i].z*v.z;
      z=nb[i].x*w.x+nb[i].y*w.y+nb[i].z*w.z;
      px[i]=x;py[i]=y;pz[i]=z;
/*      x=nb[i].x;y=nb[i].y;z=nb[i].z;*/
      xx=x*x;
      yy=y*y;
      xxx=xx*x;
      yyy=yy*y;
      A+=pow(x,4.0);
      B+=xx*yy;
      C+=xxx*y;
      D+=xxx;
      E+=xx*y;
      F+=xx;
      G+=pow(y,4.0);
      H+=x*yy;
      I+=yyy;
      J+=yy;
      K+=x*y;
      L+=x;
      M+=y;
      N+=xx*z;
      O+=yy*z;
      P+=x*y*z;
      Q+=x*z;
      R+=y*z;
      S+=z;
      T+=yyy*x;
    }  
  m[1][1]=A;
  m[1][2]=B;
  m[1][3]=C;
  m[1][4]=D;
  m[1][5]=E;
  m[1][6]=F;
  m[2][1]=B;
  m[2][2]=G;
  m[2][3]=T;
  m[2][4]=H;
  m[2][5]=I;
  m[2][6]=J;
  m[3][1]=C;
  m[3][2]=T;
  m[3][3]=B;
  m[3][4]=E;
  m[3][5]=H;
  m[3][6]=K;
  m[4][1]=D;
  m[4][2]=H;
  m[4][3]=E;
  m[4][4]=F;
  m[4][5]=K;
  m[4][6]=L;
  m[5][1]=E;
  m[5][2]=I;
  m[5][3]=H;
  m[5][4]=K;
  m[5][5]=J;
  m[5][6]=M;
  m[6][1]=F;
  m[6][2]=J;
  m[6][3]=K;
  m[6][4]=L;
  m[6][5]=M;
  m[6][6]=n;
  vec[1]=N;
  vec[2]=O;
  vec[3]=P;
  vec[4]=Q;
  vec[5]=R;
  vec[6]=S;

  locN=6;
  dd=(double *)calloc(1,sizeof(double));
  ludcmp(m,locN,indx,dd);
  lubksb(m,locN,indx,vec);
  
  a=vec[1];
  b=vec[2];
  c=vec[3];
  d=vec[4];
  f=vec[5];
  g=vec[6];
  denom=1+f*f+d*d;
  H=(b*(d*d+1)+a*(1+f*f)-c*d*f)/pow((double) denom,1.5);
  K=(4*a*b-c*c)/(denom*denom);
  
  curvature.H=H;
  curvature.K=K;
  /*  printf("Exiting est.kappa:  curvatures are: (%f, %f, %f, %f)\n",(float) H,curvature.K,curvature.k1,curvature.k2);*/
  if(H*H-K<0) printf("OOOOPS!!!\n");
 temp=sqrt((double) (H*H-K));
  curvature.k1=H-temp;
  curvature.k2=H+temp;

  /*** Now Calculate the Error of the Fit (Goodness of fit measure) *****/

  error=0;
  for(i=0;i<n;i++)
    error+=(pz[i]-a*px[i]*px[i]-b*py[i]*py[i]-c*px[i]*py[i]-d*px[i]-f*py[i]-g)*(pz[i]-a*px[i]*px[i]-b*py[i]*py[i]-c*px[i]*py[i]-d*px[i]-f*py[i]-g);
  curvature.error=error;
  free(px);
  free(py);
  free(pz);
  free(dd);
  Dfree2d(m,7);
  return(curvature);
  
}      

/*Find Matrix Inverse*/

void matrixInv(double **Ainv,double **A,int N)
{
  int i,j;
  int *indx;
  double **a,d,*col;
  
  a=Dalloc2d(N+1,N+1);
  indx=Ialloc1d(N+1);
  col=Dalloc1d(N+1);

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      a[i+1][j+1]=A[i][j];
  
  ludcmp(a,N,indx,&d);
  for(j=1;j<=N;j++)
    {
      for(i=1;i<=N;i++)
	col[i]=0.0;
      col[j]=1.0;
      lubksb(a,N,indx,col);
      for(i=1;i<=N;i++)
	Ainv[i-1][j-1]=col[i];
    }
  Dfree2d(a,N+1);
  free(indx);
  free(col);
}

/*Find Matrix Determinant*/

void matrixDet(double *d,double **A,int N)
{
  int i,j;
  int *indx;
  double **a;
  
  a=Dalloc2d(N+1,N+1);
  indx=Ialloc1d(N+1);

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      a[i+1][j+1]=A[i][j];

  ludcmp(a,N,indx,d);
  for(i=1;i<=N;i++)
    *d*=a[i][i];

  Dfree2d(a,N+1);
  free(indx);
}

/*Calculate Curvature*/

void Fvector3dcurvature2d(int NNN,int MMM,Fvector3d **p,float **curvature,int flag,float max_kap_error)
{
  int a,b,i,j,k,l,m,n,temp_size,count,start_k,start_l,end_k,end_l,size,max;
  Fvector3d *nb,dummy;
  Skappa temp;
  float maga,magb,kap_error;

  if(flag<4&&flag>=0)
    {
      max=15;
      nb=Fvector3dalloc1d((max*2+1)*(max*2+1));
      
      for(i=max;i<NNN-max-1;i++)
	for(j=max;j<MMM-max-1;j++)
	  {
	    size=4;
	    kap_error=0;
	    while(kap_error<max_kap_error&&size<=max)
	      {
		nb[0].x=0;
		nb[0].y=0;
		nb[0].z=0;
		start_k=i-size;
		start_l=j-size;
		m=1;
		n=1;
		a=1;
		b=2;
		
		nb[a].x=p[i][n+j].x-p[i][j].x;
		nb[a].y=p[i][n+j].y-p[i][j].y;
		nb[a].z=p[i][n+j].z-p[i][j].z;
		nb[b].x=p[m+i][j].x-p[i][j].x;
		nb[b].y=p[m+i][j].y-p[i][j].y;
		nb[b].z=p[m+i][j].z-p[i][j].z;
		
		end_k=start_k+2*size+1;
		end_l=start_l+2*size+1;
		count=3;
		for(k=start_k;k<end_k;k++)
		  for(l=start_l;l<end_l;l++)
		    if(k!=i||l!=j)
		      if(k!=(m+i)||l!=j)
			if(k!=i||l!=(n+j))
			  {
			    nb[count].x=p[k][l].x-p[i][j].x;
			    nb[count].y=p[k][l].y-p[i][j].y;
			    nb[count].z=p[k][l].z-p[i][j].z;
			    count++;
			  }
		/*
		  printf("count=%d",count);
		  */
		fflush(stdout);
		if(flag!=4)
		  temp=estimate_kappa(nb,count,dummy);
		kap_error=temp.error;
		/*
		  printf("kap_error=%f\n",kap_error);
		  */          
		switch (flag)
		  {
		  case 0:
		    temp.k1*=-1;
		    curvature[i][j]=temp.k1;
		    break;
		  case 1:
		    curvature[i][j]=temp.k2;
		    break;
		  case 2:
		    temp.k1*=-1;
		    curvature[i][j]=sqrt((double) (temp.k1*temp.k2)*(temp.k1*temp.k2));
		    break;
		  case 3:
		    temp.k1*=-1;
		    curvature[i][j]=temp.k1+temp.k2;
		    break;
		  default:
		    break;
		  }
		size++;
	      }
	  }
      free(nb);
    }
  else
    {
      for(i=0;i<NNN;i++)
	for(j=0;j<MMM;j++)
	  {
	    if(i%5==0||j%10==0)
	      curvature[i][j]=1;
	  }
    }
}

static float at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void svdcmp(a,m,n,w,v)
float **a,*w,**v;
int m,n;
{
	int flag,i,its,j,jj,k,l,nm;
	float c,f,h,s,x,y,z;
	float anorm=0.0,g=0.0,scale=0.0;
	double *rv1;

	if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
	rv1=vector(1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=PYTHAG(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) nrerror("No convergence in 30 SVDCMP iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}

#undef SIGN
#undef MAX
#undef PYTHAG
