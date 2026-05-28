#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>

double secant(double (*)( double, int, double),int, double, double, double, double, int*); 
double hyper_geo(double l,int m, double x);
double **find_zeros_1(int mx_l, double xc);

double hyper_geo(double l,int m, double x){

  double z=(1-x)/2;
    
  double epsilon=1.e-12;
  double ak=1;

  double F21_l=1;
  int k=1;
  
  while( fabs(ak)>epsilon ){        
    ak=ak*(m-l+k-1)*(m+l+k)*z/(k*(m+k));
    F21_l+=ak;
    k+=1;
  }
  return(F21_l);
}

double hyper_deriv(double l,int m, double x){
  return(l*x*hyper_geo(l,m,x)-(l-m)*hyper_geo(l-1,m,x));
}

double **find_zeros_1(int mx_l, double xc){

  double delta;
  double tol=1.e-12;
  int m=1;
  int flag;
  double l0=3;
  double fl,fl1,sign,csign;

  double l,l1;

  int jn;
  
  double **zs;
  zs=malloc(mx_l*sizeof(double*));
  for( jn=0; jn<mx_l; jn++){
    zs[jn]=(double *)calloc(mx_l,sizeof(double));
  }

  delta=1/(20*(1-xc));
  l=l0;
  l1=l+delta;

  for( m=0; m<=mx_l; m++){
    l=0;
    l1=l+delta;

    for( jn=0; jn<mx_l; jn++){
      if( m>jn )continue;
      zs[jn][m]=secant((&hyper_geo),m,xc, l,l1,tol, &flag);
      l=zs[jn][m]+delta;
      fl=hyper_geo(l,m,xc);
    
      l1=l+delta;
      fl1=hyper_geo(l1,m,xc);
    
      sign=(fl1-fl)/fabs(fl1-fl);
      csign=sign;
    
      while( csign==sign ){
	fl=fl1;
    
	l1+=delta;
	fl1=hyper_geo(l1,m,xc);

	csign=(fl1-fl)/fabs(fl1-fl);
      }
      l=l1;/* +delta;*/
      l1=l+delta;
    }
  }
  return(zs);
}


double **find_zeros_2(int mx_l, double xc){

  double delta;
  double tol=1.e-12;
  int m=1;
  int flag;
  double fl,fl1,sign,csign;

  double l,l1;
  int jn;
  
  double **zs;
  zs=malloc(mx_l*sizeof(double*));
  for( jn=0; jn<mx_l; jn++){
    zs[jn]=(double *)calloc(mx_l,sizeof(double));
  }

  delta=fabs(1/(20*(1-xc)));
  for( m=0; m<=mx_l; m++){

    l=0;
    l1=l+delta;
    
    for( jn=0; jn<mx_l; jn++){      
      if( m>jn ) continue;
      zs[jn][m]=secant((&hyper_deriv),m,xc, l,l1,tol, &flag);

      l=zs[jn][m]+delta;
      fl=hyper_deriv(l,m,xc);
    
      l1=l+2*delta;
      fl1=hyper_deriv(l1,m,xc);
    
      sign=(fl1-fl)/fabs(fl1-fl);
      csign=sign;
    
      while( csign==sign ){
	fl=fl1;
    
	l1+=delta;
	fl1=hyper_deriv(l1,m,xc);

	csign=(fl1-fl)/fabs(fl1-fl);
      }
      l=l1; /*+delta;*/
      l1=l+delta;
    }
  }
  return(zs);
}


double secant( double (*f)(double l, int m, double x), int m, double x, double x1,double x2, double eps, int *flag) {
  double x3;
  int i, iter=1000;
  double alpha=0.25;
  *flag = 1;
  i = 0;
  while (fabs(x2 - x1) >= eps)
    {
      i = i + 1;
      x3 = x2 - alpha*(f(x2,m,x)*(x2-x1))/(f(x2,m,x)-f(x1,m,x));
      x1 = x2;
      x2 = x3;
      if(i >= iter) break;
    }
  if (i == iter) *flag = 0;
  return x3;
} 


double G_func(double l, int m, double x){
  
  double z;
  double epsilon=1.e-10;
  double ak=1;
  double b,beta,G;
  double test=1;
  int k,n;
  
  z=(1-x)/2;    
    
  if( m==0 ){
    b=0;
  }else{
    b=0;
    for( n=0; n<2*m; n++){
      b+=(l-m+1+n);
    }   
  }
  
  G=b;
  k=1;
  while( fabs(test)>epsilon ){
    
    beta=0;
    if( k==0 ){
      beta=0;
    }else{
      for( n=0; n<k; n++ ){
	beta+=1/((m+l+1+n)*(m-l+n));
      }
      beta*= -(2*l+1);
    }
    ak=ak*(m-l+k-1)*(m+l+k)*z/(k*(m+k));
    G+=ak*(b+beta);
    k+=1;
    test=ak; /* *(b+beta); */
  }             
  return(G);
}


double calc_plm(double l, int m, double x, double xc,int p_typ){

  double F21;
  double plm;

  double delt=1.e-6;

  if( (m!=0) & ((1-x)<delt) ){
    return((double)0);
  }

  F21=hyper_geo(l,m,x);
  if( m==0 ){
    return(F21);
  }
  
  plm=pow((1-x*x),(double)m/2)*F21;

  return(plm);
}


double plm_norm(double l, int m, double xc, int p_type){

  int jp;
  double x, dx;
  double norm;
  int npts = 1000;
  double plm;

  dx=(1.0-xc)/(double)npts;
  x=xc;
  norm=0;
  for( jp=0; jp<npts; jp++ ){
    plm=calc_plm(l, m, x, xc, p_type);
    norm+=dx*plm*plm;
    x+=dx;
  }
  return(sqrt(norm));
}
  
	       
