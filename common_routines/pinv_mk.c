#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>

#define SV_THRESH 2e-4
#define DIFF_THRESH 0.
#define MIN(q,p) (((q)<(p))?(q):(p))
#define MAX(q,p) (((q)>(p))?(q):(p))

/*

ordering of input matrix is [r1c1, r1c2, r1c3,...,r2c1, r2c2...]

output matrix "inv" must be declared to be double [n*m]

*/


int pinv(int m_in, int n_in, double* a, double* inv){
  int info;
  MKL_INT m=(MKL_INT)m_in;
  MKL_INT n=(MKL_INT)n_in;
  
  MKL_INT lda=n;
  MKL_INT ldu=m;
  MKL_INT ldvt=n;
  FILE *svfl;
  
  double t1,t2;

  int jr,jc;
  int sixfour=64;

  /*  mkl_set_dynamic(0); */
  /* mkl_set_num_threads(6); */

  double *aa;
  double alpha=(double)SV_THRESH;
  aa=(double*)mkl_calloc(m*n,sizeof(double),sixfour);
  cblas_dcopy(m*n,a,(MKL_INT)1,aa,(MKL_INT)1);
  
  double *sv;
  sv=(double*)mkl_calloc((size_t)MIN(m,n),sizeof(double),sixfour);
  double *superb;
  superb=(double*)mkl_calloc((size_t)MIN(m,n),sizeof(double),sixfour);
  double *u;
  u=(double*)mkl_calloc((size_t)ldu*(size_t)m,sizeof(double),sixfour);
  double *vt;
  vt=(double*)mkl_calloc((size_t)ldvt*(size_t)n,sizeof(double),sixfour);

  t1=dsecnd();

  /* info=LAPACKE_dgesdd( LAPACK_ROW_MAJOR,'A',(lapack_int)m,(lapack_int)n,aa,(lapack_int)lda,sv,u,(lapack_int)ldu,vt,(lapack_int)ldvt);   */
  info=LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'A','A',(lapack_int)m,(lapack_int)n,aa,(lapack_int)lda,sv,u,(lapack_int)ldu,vt,(lapack_int)ldvt,superb);
  if( info != 0 ) {
    printf( "The algorithm computing SVD failed to converge. info = %d\n",info );

    mkl_free((double*)aa);
    mkl_free((double*)sv);
    mkl_free((double*)u);
    mkl_free((double*)vt);
    return( -1 );
  }
  
  t2=dsecnd();
  /* fprintf(stderr,"\n\ntime for dgesdd: %f\n",(float)((t2-t1))); */

  /* svfl=fopen("singVals","w"); */
  /* fprintf(svfl,"%lld\n",MIN(m,n)); */
  /* for(jr=0; jr<(int)MIN(m,n); jr++)fprintf(svfl," %10.6le \n",sv[jr]); */
  /* fclose(svfl); */
    
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

  int nsig=1;
  for( jc=1; jc<MIN(m,n); jc++){
    if(sv[jc]/sv[0]<SV_THRESH) break;
    nsig++;
  }
  
  /* fprintf(stderr,"NSIG: %d\n",nsig); */
  if( nsig <= 0 ){
    mkl_free((double*)aa);
    mkl_free((double*)sv);
    mkl_free((double*)u);
    mkl_free((double*)vt);
    return( -1);
  }

  double temp;  
  /* MKL_INT i1=1; */
/* #pragma omp parallel for */
/*   for( jr=0; jr<n; jr++){ */
/*     if( jr<m){ */
/*       temp=(double)(sv[jr]/(sv[jr]*sv[jr]+alpha*alpha)); */
/*     }else{ */
/*       temp=0; */
/*     } */
/*     cblas_dscal(ldvt,temp,(vt+jr*ldvt),i1); */
/*   } */

  double *vsig;
  vsig=(double*)mkl_calloc((size_t)n*(size_t)m,sizeof(double),sixfour);

  /* if( m<=n){ */
  alpha=alpha*sv[0];
#pragma omp parallel for
  for( jr=0; jr<n; jr++)for(jc=0; jc<MIN(n,m); jc++ ){
      temp=(double)(sv[jc]/(sv[jc]*sv[jc]+alpha*alpha));
      vsig[jr*m+jc]=temp*vt[jr+jc*n];
    }
  
  t1=dsecnd();
  /* fprintf(stderr,"pinv going to dgemm %lld %lld %lld %lld\n",n,m,ldvt,ldu); */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, m, 1.0, vsig, m, u, m, 0.0, inv, m);

  t2=dsecnd();
  /* fprintf(stderr,"time for cblas_dgemm: %f\n\n\n",(float)((t2-t1))); */

  mkl_free((double*)aa);
  mkl_free((double*)sv);
  mkl_free((double*)u);
  mkl_free((double*)vt);
  mkl_free((double*)vsig);
  return(0);
}

