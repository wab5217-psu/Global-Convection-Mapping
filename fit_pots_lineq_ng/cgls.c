#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>

#define DIFF_THRESH 0.
#define MIN(q,p) (((q)<(p))?(q):(p))
#define MAX(q,p) (((q)>(p))?(q):(p))

int cgls(int m,int n,double* a[],double* b, double* x){

  int stat;
  int jr,jc;
  int sixfour=64;
  double afac_m=-1;
  double afac_p=1;
  double bfac=1;
  double bfac_z=0;
  int step=1;
  double delta = 1;
  double delta_x=1;
  double res_norm;
  double res_norm_last;
  double x_norm;
  double x_norm_last;

  
  double *res;
  res=(double*)mkl_calloc(m,sizeof(double),sixfour);

  double *d;
  d=(double*)mkl_calloc(m,sizeof(double),sixfour);


  double *A;
  A=(double*)mkl_calloc(n*m,sizeof(double),sixfour);
  for( jc=0; jc<n; jc++ )for( jr=0; jr<m; jr++ )*(A+jr+jc*m)=a[jr][jc];
  
  
  double *ATR;
  ATR=(double*)mkl_calloc(n,sizeof(double),sixfour);
  double ATR_norm_l;
  double ATR_norm;
  
  double *AD;
  AD=(double*)mkl_calloc(m,sizeof(double),sixfour);
  double AD_norm;
  
  double alpha=1;
  double beta;
  
  fprintf(stderr,"in cgls %d rows %d columns\n",m,n);
  
  /*  res = b - A * x  */
  
  for( jr=0; jr<m; jr++ )res[jr]=b[jr];
  cblas_dgemv(CblasRowMajor,CblasNoTrans,m,n,afac_m,A,n,x,step,bfac,res,step);
  res_norm=cblas_dnrm2(m,res,step);

  x_norm=cblas_dnrm2(n,x,step);
  
  /*  d = AT * res  */
  for( jr=0; jr<n; jr++ )d[jr]=0;
  cblas_dgemv(CblasRowMajor,CblasTrans,m,n,afac_p,A,n,res,step,bfac_z,d,step);

  fprintf(stderr,"in cgls dgemv 2\n");

  int itmax=500;
  int itr=0;
  int imax;
  /* double dtol=.2; */
  double dtol=1e-6;
  double dmax=100;
  double xmax;
  double davg;
  double rmax=10;
  double drmax=100;
  double ravg;
  int nan_flag=0;
  
  FILE *fp;

  fp=fopen("deltas","w");

  /* while( alpha*drmax/rmax > dtol && itr<itmax){ */
  while( alpha*drmax/x_norm > dtol && itr<itmax){
    itr++;
    
    /*  alpha = ||AT Res ||^2 / ||A d ||^2  */
    
    cblas_dgemv(CblasRowMajor,CblasTrans,m,n,afac_p,A,n,res,step,bfac_z,ATR,step);
    ATR_norm_l=cblas_dnrm2(n,ATR,step);
    
    cblas_dgemv(CblasRowMajor,CblasNoTrans,m,n,afac_p,A,n,d,step,bfac_z,AD,step);
    AD_norm=cblas_dnrm2(m,AD,step);

    alpha=(ATR_norm_l*ATR_norm_l)/(AD_norm*AD_norm);
    
    fprintf(stderr,"cgls iter:%d %lf %lf\n",itr,AD_norm,ATR_norm);

    nan_flag=0;
    for( jc=0; jc<n; jc++ )if(!isfinite(d[jc]) || !isfinite(alpha)){
	fprintf(stderr,"****************************************had a nan: %lf %lf %lf\n",d[jc],alpha,AD_norm);
	nan_flag=1;
	break;
      }

    if( nan_flag )break;
    
    /*  x = x + alpha * d  */
    cblas_daxpy(n,alpha,d,1,x,1);
    
    /* res = res - alpha * A * d  */
    cblas_daxpy(m,(-1)*alpha,AD,1,res,1);

  /*  beta = ||AT res ||^2 /|| AT res_l ||^2   */    
    cblas_dgemv(CblasRowMajor,CblasTrans,m,n,afac_p,A,n,res,step,bfac_z,ATR,step);
    ATR_norm=cblas_dnrm2(n,ATR,step);
    
    beta=(ATR_norm*ATR_norm)/(ATR_norm_l*ATR_norm_l);
    
  /*  d = AT * res  + beta * d  */
    for( jc=0; jc<n; jc++ )d[jc] = beta*d[jc]+ATR[jc];
    
    res_norm_last=res_norm;
    
    res_norm=cblas_dnrm2(m,res,step);

    x_norm_last=x_norm;
    x_norm=cblas_dnrm2(n,x,step);

    imax=cblas_idamax(m,res,1);
    rmax=fabs(res[imax]);
    imax=cblas_idamax(m,AD,1);
    drmax=fabs(AD[imax]);
    
    delta = x_norm/res_norm;
    delta_x=(x_norm-x_norm_last)/x_norm;
    fprintf(fp,"%d  %le %le %le %le %le\n",itr,alpha*drmax,rmax,alpha*drmax/rmax,res_norm,alpha*drmax/x_norm);
    fflush(fp);

  }
  

  mkl_free(res);
  fprintf(stderr,"res free\n");
  mkl_free(d);
  fprintf(stderr,"d free\n");
  mkl_free(A);
  fprintf(stderr,"A free\n");
  mkl_free(ATR);
  fprintf(stderr,"ATR free\n");
  mkl_free(AD);
  fprintf(stderr,"AD free\n");

  fclose(fp);
  stat=0;
  return(stat);
}
