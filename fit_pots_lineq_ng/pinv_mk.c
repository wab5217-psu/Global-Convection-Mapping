#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>

#define SV_THRESH 1e-8
#define DIFF_THRESH 0.
#define MIN(q,p) (((q)<(p))?(q):(p))
#define MAX(q,p) (((q)>(p))?(q):(p))

/*

ordering of input matrix is [r1c1, r1c2, r1c3,...,r2c1, r2c2...]

output matrix "inv" must be declared to be double [n*m]

*/


int pinv(MKL_INT m, MKL_INT n, double* a, double* inv){
  int info;
  MKL_INT lda=n;
  MKL_INT ldu=m;
  MKL_INT ldvt=n;
  FILE *svfl;
  FILE *ufl,*vfl;
  
  double t1,t2;
  double alpha;

  int jr,jc;
  int sixfour=64;

  /*  mkl_set_dynamic(0); */
  /* mkl_set_num_threads(6); */

  double *aa;
  aa=(double*)mkl_calloc(m*n,sizeof(double),sixfour);
  cblas_dcopy(m*n,a,(MKL_INT)1,aa,(MKL_INT)1);
  
  double *sv;
  sv=(double*)mkl_calloc(MIN(m,n),sizeof(double),sixfour);
  double *u;
  u=(double*)mkl_calloc(ldu*m,sizeof(double),sixfour);
  double *vt;
  vt=(double*)mkl_calloc(ldvt*n,sizeof(double),sixfour);

  t1=dsecnd();

  info=LAPACKE_dgesdd( LAPACK_ROW_MAJOR,'A',m,n,aa,lda,sv,u,ldu,vt,ldvt);  

  if( info > 0 ) {
    printf( "The algorithm computing SVD failed to converge.\n" );
    return( -1 );
  }
  
  t2=dsecnd();
  fprintf(stderr,"\n\ntime for dgesdd: %f\n",(float)((t2-t1)));

  svfl=fopen("singVals","w");
  fprintf(svfl,"%d\n",MIN(m,n));
  for(jr=0; jr<MIN(m,n); jr++)fprintf(svfl," %10.4lf ",sv[jr]);
  fclose(svfl);
    
  /* Check for convergence */
  if( info > 0 ) {
    fprintf(stderr,"The algorithm computing SVD failed to converge. info=%d;\n",info );
    fprintf(stderr,"the least squares solution could not be computed.\n" );

    mkl_free((double*)aa);
    mkl_free((double*)sv);
    mkl_free((double*)u);
    mkl_free((double*)vt);
    return( -1 );
  }


  /* /\* Determine SV threshold *\/ */
  /* int dif_sp=2; /\* spacing between points for smoothed derivative *\/ */
  /* int npts=MIN(m,n); */

  /* double avg_diff=0; */
  /* double *diffs; */
  /* int nsig=0; */
  /* diffs=calloc(npts,sizeof(double)); */
  /* for( jc=0; jc<npts; jc++){ */
  /*   diffs[jc]= */
  /*   avg_diff+=diffs[jc]; */
  /* } */

  
  /* avg_diff/=(double)jc; */
  /* for( jc=npts/2; jc<npts; jc++ ){ */
  /*   if(log(diffs[jc]/avg_diff) > (double)DIFF_THRESH) break; */
  /*   nsig=jc; */
  /*   /\*    nsig=jc*.98; *\/ /\*sometimes needed...like for the potential mapping *\/ */
  /* } */
  /* free(diffs); */

  int nsig=1;
  for( jc=1; jc<MIN(m,n); jc++){
    if(sv[jc]/sv[0]<SV_THRESH) break;
    nsig++;
  }
  
  fprintf(stderr,"NSIG: %d\n",nsig);
  if( nsig <= 0 ){
    mkl_free((double*)aa);
    mkl_free((double*)sv);
    mkl_free((double*)u);
    mkl_free((double*)vt);
    return( -1);
  }

  double temp;
  alpha = SV_THRESH;
  MKL_INT i1=1;
  for( jr=0; jr<n; jr++)if(sv[jr]!=0){
      temp=(double)(sv[jr]/(sv[jr]*sv[jr]+alpha*alpha));
    cblas_dscal(ldvt,temp,(vt+jr*ldvt),i1);
  }
  /* for( jr=nsig; jr<n; jr++)if(sv[jr]!=0){ */
  /*   temp=(double)(0.); */
  /*   cblas_dscal(ldvt,temp,(vt+jr*ldvt),i1); */
  /* } */

  t1=dsecnd();

  double one=1;
  double zero=0;

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,n,m,n,one,vt,ldvt,u,ldu,zero,inv,m);

  t2=dsecnd();
  fprintf(stderr,"time for cblas_dgemm: %f\n\n\n",(float)((t2-t1)));

  mkl_free((double*)aa);
  mkl_free((double*)sv);
  mkl_free((double*)u);
  mkl_free((double*)vt);

  return(0);
}

