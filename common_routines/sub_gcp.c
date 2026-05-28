#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define RE 6371.
#define PI acos(-1)
#define sind(x) (sin(fmod((x),360)*PI/180))
#define cosd(x) (cos(fmod((x),360)*PI/180))
#define tand(x) (tan(fmod((x),360)*PI/180))

/* 

Determines points along the great circle path from pt 1 to pt 2.
Calculates lat and lon of a point at a distance "range" form point 1 along the path
Returns 0 if successful -1 if there is an error 

Routine based upon Wikipedia page on Great-circle Naviation
 
 */
int sub_gcp(double lon1, double lat1, double lon2, double lat2, double range, double *lon, double *lat){

  double lon12=lon2-lon1;
  double alpha1=atan2(sind(lon12),cosd(lat1)*tand(lat2)-sind(lat1)*cosd(lon12));
  double alpha0=atan2(sin(alpha1)*cosd(lat1),
		      sqrt(cos(alpha1)*cos(alpha1)+sin(alpha1)*sin(alpha1)*sind(lat1)*sind(lat1)));
  double sig01=atan2(tand(lat1),cos(alpha1));
  double lon0=lon1-atan2(sin(alpha0)*sin(sig01),cos(sig01))*180/PI;
  double sig=sig01+range/RE;
  *lon=lon0+atan2(sin(alpha0)*sin(sig),cos(sig))*180/PI;
  *lat=atan2(cos(alpha0)*sin(sig),sqrt(cos(sig)*cos(sig)+sin(alpha0)*sin(alpha0)*sin(sig)*sin(sig)))*180/PI;

  return 0;
}
