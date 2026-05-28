#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern int sub_sphazm(double Alon,double Alat,double Clon,double Clat,double *azm,double *range);

int sub_sphcal(double Alon, double Alat, double azm, double range, double *Clon, double *Clat){

  double A;
  double B;
  double gcd;
  double Re=6362.;
  double aside,aside_r;
  double bside_r;
  double cside,cside_r;
  double A_r;
  double arg;
  double dtor;
  double ck_azm,ck_range;
  int stat=1;

  dtor=acos(-1.)/180.;
  A=azm;
  gcd=range;
  cside=90.-Alat;

  A_r=A*dtor;
  cside_r=cside*dtor;
  bside_r=gcd/Re;

  arg=cos(bside_r)*cos(cside_r)+sin(bside_r)*sin(cside_r)*cos(A_r);
  aside_r=acos(arg);

  aside=aside_r/dtor;
  *Clat=90.-aside;

  arg=(sin(A_r)*sin(bside_r))/sin(aside_r);
  B =asin(arg)/dtor;
  *Clon=Alon+B;

  stat=sub_sphazm(Alon,Alat,*Clon,*Clat,&ck_azm,&ck_range);

  /* fprintf(stderr,"sub_sphcal: %f %f %f %f %f %f\n",Alon,Alat,*Clon,*Clat,range,ck_range); */
  
  if (fabs(fabs(ck_range)-fabs(range))>1.){*Clon=Alon+(180.-B);}

  while(*Clon > 360.){*Clon -=360.;}
  while(*Clon < 0.){*Clon += 360.;}

  return(stat);
}
