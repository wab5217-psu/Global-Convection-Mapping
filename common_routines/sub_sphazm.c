#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int sub_sphazm(double Alon,double Alat,double Clon,double Clat,double *azm,double *range){
  /*
; subroutine to determine azm and range of point (Clon,Clat) from
; source point (Alon,Alat)
;
; last written: Dec 8/94, Apr 11/95 - azm range
  */

  double A;
  double B;
  double Re=6362.;
  double aside,aside_r;
  double bside_r;
  double cside,cside_r;
  double B_r;
  double arg;
  double dtor;
  int stat=1;


  dtor=acos(-1.)/180.;


  cside=90.-Alat;
  aside=90.-Clat;

  B=Clon-Alon;
  cside_r=cside*dtor;
  aside_r=aside*dtor;
  B_r=B*dtor;

  arg=cos(aside_r)*cos(cside_r)+sin(aside_r)*sin(cside_r)*cos(B_r);

  /*
;if(arg >  1.){arg=1.;}
;if(arg < -1.){arg= -1.;}
  */

  bside_r=acos(arg);
  *range=Re*bside_r;

  arg=sin(B_r)*sin(aside_r)/sin(bside_r);

  /*
;if(arg >  1.){arg=1.;}
;if(arg < -1.){arg=-1.;}
  */
  A=asin(arg)/dtor;

  if (Clat < Alat){A=180.-A;}

  *azm=A;
  return(stat);
    }
