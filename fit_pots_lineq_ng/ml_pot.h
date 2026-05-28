typedef struct VDSTR{
  int yr;
  int mo;
  int dy;
  int hr;
  int mt;
  int sc;
  float Bx,By,Bz,Au,Al,v_sw;
  int npts;
  double* lats;
  double* lons;
  double* vn;
  double* ve;
  double* vn_cov;
  double* ve_cov;
  double* vmag;
  double* vaz;
  int* dcount;
}VDstr;
