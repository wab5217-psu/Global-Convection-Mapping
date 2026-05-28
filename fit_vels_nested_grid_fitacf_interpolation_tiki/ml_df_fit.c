/* #include <stdlib.h> */
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <zlib.h>
#include "rtypes.h"
#include "dmap.h"
#include "option.h"
#include "rtime.h"
#include "radar.h"
#include "rprm.h"
#include "rpos.h"
#include "fitdata.h"
#include "cfitdata.h"
#include "scandata.h"
#include "fitread.h"
#include "fitscan.h"
#include "fitindex.h"
#include "fitseek.h"
#include "rtypes.h"
#include "dmap.h"
#include "invmag.h"
#include "griddata.h"
#include "gridread.h"

#include "ml_df.h"

#include "cnvgrid.h"
#include "cnvmap.h"
#include "cnvmapindex.h"
#include "cnvmapseek.h"
#include "cnvmapread.h"
#include "cnvmapsolve.h"
#include "aacgmlib_v2.h"
#include "aacgm.h"
#include "mlt_v2.h"
#include "igrflib.h"
#include "grid.h"

#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>


extern int sub_sphazm(double Alon,double Alat,double Clon,double Clat,double *azm,double *range);
extern int sub_sphcal(double Alon, double Alat, double azm, double range, double *Clon, double *Clat);

double **find_zeros_1(int mx_l, double xc);
double **find_zeros_2(int mx_l, double xc);
double calc_plm(double l, int m, double x, double xc,int p_typ);
double plm_norm(double l, int m, double xc,int p_typ);

#define LENGTH(x,y) sqrt(x*x+y*y)

#define INERTIAL 0

#define C 299792458.0 
#define PI 3.14159265359
#define LRE 6371.
#define MIN_RANGE 600
#define MAX_RANGE 3000
#define SECS_P_DAY 86400 /* 24*60*60 */
#define MAX(q,p) (((q)>(p))?(q):(p))
#define sind(x) (sin(fmod((x),360)*PI/180))
#define cosd(x) (cos(fmod((x),360)*PI/180))
#define tand(x) (tan(fmod((x),360)*PI/180))
/* #define MIN_ERR 100.0 */
#define MIN_ERR 50.0
#define ERR_SCALE 1.0
#define MIN_ML_ERR 300.0
/* #define MAX_ML_COV 10000.0 */
#define MAX_ML_COV 5000.0
/* #define DIV_ERR 0.5 */
#define DIV_ERR 0.2
/* #define DIV_ERR 1 */
#define MAX_V 5000.0
#define MAX_V_ERR 200.0
/* #define MAX_V_ERR 0.2 */
#define MIN_V 30.
#define MIN_COUNT 5
#define MAX_BEAMS 24
#define MIN_LAT 55
#define MAX_LAT 89

#define ML_LAT_0 55
#define ML_DLAT 2
#define ML_DLON 15
#define ML_NLON 24

#define BAD_VALUE -99999.9
#define BAD_INT -9999

#define WRITE_COEF 0

#define MAX_DATA_LAT 89

#define TRUE  1
#define FALSE 0

/* Externals from grid lib */

double min_lat=MIN_LAT;
double max_lat=MAX_LAT;

double dlat;
double dlon;
double dlat_ng;
double dlon_ng;
double min_lon_ng;
double max_lon_ng;
double *start_lon;
int ngrid;
int nedge;

int nbp;
B_POINT* bp=NULL;

CELL* grid;
NEIGHBOR* neighbors;

/***************************/

double max_range=MAX_RANGE;
double min_range=MIN_RANGE;
long start_time;
long end_time;
int avg_ival;
char *radar_list[30];
int nrad;
M_ARRAY* kazm_array;
int *edge; 
M_ARRAY* edge_array;
double* edge_data;
double* los_data;
double* los_kazm;
double* los_lats;
double* los_lons;
double* los_err;
M_ARRAY* div_array;
M_ARRAY* smooth_array;
int mdays[]={31,59,90,120,151,181,212,243,273,304,334};
FILE *grd_file;
char ml_file[128];
char hemisphere[16]="north";

#define MX_L 10
double **ls_1,**ls_2;
double **norms_1,**norms_2;


int dayofweek(int d, int m, int y)
{
    static int t[] = { 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 };
    y -= m < 3;
    return ( y + y/4 - y/100 + y/400 + t[m-1] + d) % 7;
}
 
char *choppy( char *s )
{
  char *n = malloc( strlen( s ? s : "\n" ) );
  if( s )
    strcpy( n, s );
  if( s[strlen(s)-1] =='\n')
    n[strlen(n)-1]='\0';
  return n;
}


void parse_instructions(FILE *fp)
{
  char *line=NULL;
  size_t len=0;
  char *param=NULL;
  char *token;

  nrad=0;
  nbp=0;
  while( getline(&line,&len,fp)!= EOF ){
    if(line[0]=='#' || line[0]==' ')continue;
    param=strtok(line," ");
    if(strcmp(param,"boundary_point")==0)
      {
	nbp++;
	bp=realloc(bp,nbp*sizeof(B_POINT));
	sscanf(strtok(NULL," "),"%lf",&bp[nbp-1].lat);
	sscanf(strtok(NULL," "),"%lf",&bp[nbp-1].lon);
	if( bp[nbp-1].lon<0 )bp[nbp-1].lon+=360;
      }
   else if(strcmp(param,"min_lat")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&min_lat);
     }
   else if(strcmp(param,"max_lat")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&max_lat);
     }
   else if(strcmp(param,"lat_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlat);
     }
   else if(strcmp(param,"lon_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlon);
     }
   else if(strcmp(param,"ng_lat_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlat_ng);
     }
   else if(strcmp(param,"ng_lon_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlon_ng);
     }
   else if(strcmp(param,"max_range")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&max_range);
     }
   else if(strcmp(param,"min_range")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&min_range);
     }
   else if(strcmp(param,"start_time")==0)
     {
	sscanf(strtok(NULL," "),"%ld",&start_time);
     }
   else if(strcmp(param,"end_time")==0)
     {
	sscanf(strtok(NULL," "),"%ld",&end_time);
     }
   else if(strcmp(param,"avg_interval")==0)
     {
	sscanf(strtok(NULL," "),"%d",&avg_ival);
     }
   else if(strcmp(param,"ml_file")==0)
     {
	sscanf(strtok(NULL," "),"%s",ml_file);
     }
   else if(strcmp(param,"hemisphere")==0)
     {
	sscanf(strtok(NULL," "),"%s",hemisphere);
     }
   else if(strcmp(param,"radar_list")==0)
     {
       while( (token=strtok(NULL," ")) != NULL)
	 {
	   radar_list[nrad]=choppy(token);
	   nrad++;
	 }
     }
 }
  if(line)
    free(line);
}



void make_div_ar(){
  int igr,inr0,inr1,inr2;
  double fr0,fr1,fr2,chk;
  int idv;
  int ind;
  double dx,dy,dlon_l;
  double dlon0;
  FILE *divfile;

  divfile=fopen("divergence_info.dat","w");
  

  idv=nedge;
  div_array=calloc(ngrid-nedge,sizeof(struct mod_array));
  for( igr=0; igr<ngrid; igr++)div_array[igr].coef=calloc(2*ngrid,sizeof(double));

  /* first nedge points are edge points, so start calculating divergence after those */
  for( igr=nedge; igr<ngrid; igr++){
    /* calculate d/dx term as vx_cell_right - vx_cell
       if edge of nested grid area, use half from lower half from upper */    
    
    inr0=neighbors[igr].right[0];
    inr1=neighbors[igr].right[1];
    
    dlon_l=grid[igr].lon[1]-grid[igr].lon[0];
      
    dx=fabs(dlon_l*DTOR*LRE*cosd(grid[igr].center_lat));
    div_array[idv].coef[igr]=1/dx;
    if( inr1==BAD_INT){
      div_array[idv].coef[inr0]=-1/dx;
    }else{
      div_array[idv].coef[inr0]=-1/(2*dx);
      div_array[idv].coef[inr1]=-1/(2*dx);
    }
    fprintf(divfile,"%d %d %d %f ",igr,inr0,inr1,dx);
    
    
    /* calculate d/dy term as vy_cell - vy_cell_lower */

    dy=fabs((grid[igr].lat[2]-grid[igr].lat[1])*DTOR*LRE);    
    div_array[idv].coef[igr+ngrid]=1/dy;

    inr0=neighbors[igr].lower[0];
    inr1=neighbors[igr].lower[1];
    inr2=neighbors[igr].lower[2];
    fr0=0;
    fr1=0;
    fr2=0;
    if( inr0 == BAD_INT ){ 
      div_array[idv].coef[inr1+ngrid]=-1/dy;
    }else if(inr1 == BAD_INT){
      fr0=1;
      div_array[idv].coef[inr0+ngrid]=-1/dy;
    } else if( inr2==BAD_INT) { 
      dlon0=(grid[inr0].lon[2]-grid[igr].lon[0]);
      if( dlon0>360 )dlon0-=360;
      fr0=dlon0/(grid[igr].lon[1]-grid[igr].lon[0]);
      fr1=1-fr0;
      div_array[idv].coef[inr0+ngrid]=-fr0/dy;
      div_array[idv].coef[inr1+ngrid]=-fr1/dy;
    }else{
      dlon0=(grid[inr0].lon[2]-grid[igr].lon[0]);
      if( dlon0>360 )dlon0-=360;
      fr0=dlon0/(grid[igr].lon[1]-grid[igr].lon[0]);
      
      dlon0=(grid[inr1].lon[2]-grid[inr1].lon[3]);
      if( dlon0>360 )dlon0-=360;
      fr1=dlon0/(grid[igr].lon[1]-grid[igr].lon[0]);
      
      fr2=1-(fr0+fr1);
      div_array[idv].coef[inr0+ngrid]=-fr0/dy;
      div_array[idv].coef[inr1+ngrid]=-fr1/dy;
      div_array[idv].coef[inr2+ngrid]=-fr2/dy;      
    }
    fprintf(divfile,"%d %d %d %f %f %f %f\n",inr0,inr1,inr2,dy,fr0,fr1,fr2);
    if( fabs(fr0+fr1+fr2-1)>.001 )fprintf(stderr,"make_div_ar: %d %d %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d\n",igr,idv,dx,dy,fr0,fr1,fr2,inr0,inr1,inr2);

    chk=0;
    for(ind=0; ind<2*ngrid; ind++){chk+=div_array[idv].coef[ind];}
    if( fabs(chk)>.01 )fprintf(stderr,"make_div_ar: %d check: %lf\n",idv,chk);
	
    idv++;
  }
  fclose(divfile);
}


struct tm *parse_date_str( long t_i)
{
  struct tm *t_o;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();

  t_o=malloc(sizeof(struct tm));
  t_o->tm_year=(int)(t_i/1e8);
  t_o->tm_mon=(int)((t_i-1e8*t_o->tm_year)/1e6);
  t_o->tm_mday=(int)((t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon)/1e4);
  t_o->tm_hour=(int)((t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon-t_o->tm_mday*1e4)/1e2);
  t_o->tm_min=(int)(t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon-t_o->tm_mday*1e4-t_o->tm_hour*1e2);
  t_o->tm_yday=mdays[t_o->tm_mon-1]+t_o->tm_mday;
  t_o->tm_wday=dayofweek(t_o->tm_mday,t_o->tm_mon,t_o->tm_year);
  if( IS_LEAPYEAR(t_o->tm_year) && t_o->tm_mon>2) t_o->tm_yday++;
  t_o->tm_sec=0;
  t_o->tm_year-=1900;
  t_o->tm_mon-=1;
  t_o->tm_isdst=0;
  return t_o;
}

time_t fname_to_time(char *fname)
{
  int datev;
  char datestr[10]="0";
  char hr_str[3]="0";
  char mn_str[3]="0";
  struct tm f_tm;
  char *tz;
  
  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();
  
  memset(datestr, '\0', sizeof datestr);
  memset(hr_str, '\0', sizeof hr_str);
  memset(mn_str, '\0', sizeof mn_str);

  strncpy(datestr,fname,8);
  strncpy(hr_str,fname+9,2);
  strncpy(mn_str,fname+11,2);
  datev=atoi(datestr);
  f_tm.tm_year=datev/10000;
  f_tm.tm_mon=(datev-10000*f_tm.tm_year)/100;
  f_tm.tm_mday=datev-10000*f_tm.tm_year-100*f_tm.tm_mon;
  f_tm.tm_yday=mdays[f_tm.tm_mon-1]+f_tm.tm_mday;
  if( IS_LEAPYEAR(f_tm.tm_year) && f_tm.tm_mon>2) f_tm.tm_yday++;
  f_tm.tm_hour=atoi(hr_str);
  f_tm.tm_min=atoi(mn_str);
  f_tm.tm_sec=0;
  f_tm.tm_isdst=0;
  f_tm.tm_wday=dayofweek(f_tm.tm_mday,f_tm.tm_mon,f_tm.tm_year);
  f_tm.tm_year-=1900; /* unix epoch year correction */
  f_tm.tm_mon-=1;     /* unix epoch month 0 to 11 */

  return mktime(&f_tm);
}

int split_dateline(char* date_line, MLDstr *mlDstr){

  char * token0 = strtok(date_line," ");
  char * token1 = strtok(NULL," ");

  sscanf(token0,"%d-%d-%d",&mlDstr->yr, &mlDstr->mo, &mlDstr->dy);
  sscanf(token1,"%d:%d:%d",&mlDstr->hr, &mlDstr->mt, &mlDstr->sc);

  return(0);
}


int read_ml_record(FILE *fp, MLDstr *mlDstr, char *hemisphere){

  double n_rmse,e_rmse;
  char date_line[256];
  int stat;
  int npts;
  int i;

  if( fgets(date_line, sizeof(date_line), fp)==NULL ){
    return(-1);
  }

  
    stat=split_dateline(date_line,mlDstr);
    fprintf(stderr,"date: %d %d %d\n  time: %d %d %d\n",mlDstr->yr,mlDstr->mo,mlDstr->dy,mlDstr->hr,mlDstr->mt,mlDstr->sc);

    fscanf(fp,"%f %f %f %f %f %f\n",&mlDstr->Bx,&mlDstr->By,&mlDstr->Bz,&mlDstr->Au,&mlDstr->Al,&mlDstr->v_sw);

    fscanf(fp,"%d",&npts);
    fprintf(stderr,"npts %d \n",npts);


    if( mlDstr->lats!=NULL )free(mlDstr->lats);
    mlDstr->lats=(double *)calloc(npts,sizeof(double));
    if( mlDstr->lons!=NULL) free(mlDstr->lons);
    mlDstr->lons=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vn!=NULL) free(mlDstr->vn);
    mlDstr->vn=(double *)calloc(npts,sizeof(double));
    if( mlDstr->ve!=NULL) free(mlDstr->ve);
    mlDstr->ve=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vn_cov!=NULL) free(mlDstr->vn_cov);
    mlDstr->vn_cov=(double *)calloc(npts,sizeof(double));
    if( mlDstr->ve_cov!=NULL) free(mlDstr->ve_cov);
    mlDstr->ve_cov=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vmag!=NULL) free(mlDstr->vmag);
    mlDstr->vmag=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vaz!=NULL) free(mlDstr->vaz);
    mlDstr->vaz=(double *)calloc(npts,sizeof(double));

   
    for( i=0; i<npts; i++ ){
      if((stat=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
		      (mlDstr->lats+i),(mlDstr->lons+i),(mlDstr->vn+i),
		      (mlDstr->ve+i),&n_rmse,&e_rmse,
		      (mlDstr->vmag+i),(mlDstr->vaz+i)))==EOF)
	{return(-1);}
      mlDstr->vn_cov[i]=MIN(n_rmse*n_rmse,MAX_ML_COV);
      mlDstr->ve_cov[i]=MIN(e_rmse*e_rmse,MAX_ML_COV);

      if( strcmp(hemisphere,"south")==0 ){
	*(mlDstr->vaz+i)+=180;
	*(mlDstr->ve+i)*=-1.0;
	*(mlDstr->vn+i)*=-1.0;
      }
      

      mlDstr->lons[i]=(float)inv_MLTConvertYMDHMS_v2(mlDstr->yr,
						     mlDstr->mo,mlDstr->dy,
						     mlDstr->hr,mlDstr->mt,
						     mlDstr->sc,mlDstr->lons[i]);

      
      if( mlDstr->lons[i]<0 )mlDstr->lons[i]+=360;
    }
    mlDstr->npts=npts;
    return(ftell(fp));
}


int fit_ml(MLDstr *mlDstr, double *coefs){

  MKL_INT n_rows, n_cols, nrhs = 1, lda, ldb, info;

  int mx_l=MX_L;
  int jp,jl,m,col,ncoef;
  int sixfour=64;

  double lat,lon,x,plm;

  double theta_c=(90.0-mlDstr->lats[0])*DTOR;
  double xc=cos(theta_c);

  ncoef=0;
  for( jl=0; jl<mx_l; jl++ )for( m=0; m<=jl; m++)ncoef+=1;  

  fprintf(stderr,"fit_ml npts: %d\n",mlDstr->npts);
  
  n_rows=mlDstr->npts-24; // 24 is number of hours in model
  n_cols=2*ncoef-mx_l;
  lda=n_cols;
  ldb=1;

  
  double *v1;
  v1=(double*)mkl_calloc(n_rows,sizeof(double),sixfour);

  double *v2;
  v2=(double*)mkl_calloc(n_rows,sizeof(double),sixfour);
  
  double *A;
  A=(double*)mkl_calloc(n_rows*n_cols,sizeof(double),sixfour);

  int aindx=0;

  for( jp=0; jp<n_rows; jp++ ) *(v1+jp)=mlDstr->vn[jp+24];
  for( jp=0; jp<n_rows-24; jp++ ) *(v2+jp)=mlDstr->ve[jp+24];
  
  for( jp=24; jp<mlDstr->npts; jp++ ){
    col=0;
    for( jl=0; jl<mx_l; jl++ )for( m=0; m<=jl; m++){
	if( m==0 ){
	  
	  aindx=(jp-24)*n_cols+col;
	  col+=1;
	  lat=mlDstr->lats[jp];
	  lon=mlDstr->lons[jp];
	  
	  x=cos((90.0-lat)*DTOR);
	  
	  plm=calc_plm(ls_1[jl][m],m,x,xc,1)/norms_1[jl][m];
	  
	  A[aindx]=cos((double)m*lon*DTOR)*plm;
	}else{
	  aindx=(jp-24)*n_cols+col;
	  col+=2;
	  lat=mlDstr->lats[jp];
	  lon=mlDstr->lons[jp];
	  
	  x=cos((90.0-lat)*DTOR);
	  
	  plm=calc_plm(ls_1[jl][m],m,x,xc,1)/norms_1[jl][m];
	  
	  A[aindx]=cos((double)m*lon*DTOR)*plm;	    
	  A[aindx+1]=sin((double)m*lon*DTOR)*plm;	    
	}
      }
  }
  fprintf(stderr,"fit_ml %d %d %d\n",n_rows,n_cols,nrhs);
  info = LAPACKE_dgels( LAPACK_ROW_MAJOR, 'N', n_rows, n_cols, nrhs, A, lda, v1, ldb );
  
  /* Check for the full rank */
  if( info > 0 ) {
    fprintf(stderr,"The diagonal element %lld of the triangular factor ", info );
    fprintf(stderr,"of A is zero, so that A does not have full rank;\n" );
    fprintf(stderr,"the least squares solution could not be computed.\n" );
    exit( 1 );
  }
  
  for( jp=0; jp<n_cols; jp++ )coefs[jp]=v1[jp];
				
  aindx=0;    
  for( jp=24; jp<mlDstr->npts; jp++ ){
    col=0;
    for( jl=0; jl<mx_l; jl++ )for( m=0; m<=jl; m++){
	if( m==0 ){
	  
	  aindx=(jp-24)*n_cols+col;
	  col+=1;
	  lat=mlDstr->lats[jp];
	  lon=mlDstr->lons[jp];
	  
	  x=cos((90.0-lat)*DTOR);
	  
	  plm=calc_plm(ls_2[jl][m],m,x,xc,1)/norms_2[jl][m];
	  
	  A[aindx]=cos((double)m*lon*DTOR)*plm;
	}else{
	  aindx=(jp-24)*n_cols+col;
	  col+=2;
	  lat=mlDstr->lats[jp];
	  lon=mlDstr->lons[jp];
	  
	  x=cos((90.0-lat)*DTOR);
	  
	  plm=calc_plm(ls_2[jl][m],m,x,xc,1)/norms_2[jl][m];
	  
	  A[aindx]=cos((double)m*lon*DTOR)*plm;	    
	  A[aindx+1]=sin((double)m*lon*DTOR)*plm;
	  
	}
      }
  }
  fprintf(stderr,"fit_ml 2 %d %d %d\n",n_rows,n_cols,nrhs);
  info = LAPACKE_dgels( LAPACK_ROW_MAJOR, 'N', n_rows, n_cols, nrhs, A, lda, v2, ldb );
  
    /* Check for the full rank */
  if( info > 0 ) {
    fprintf(stderr,"The diagonal element %lld of the triangular factor ", info );
    fprintf(stderr,"of A is zero, so that A does not have full rank;\n" );
    fprintf(stderr,"the least squares solution could not be computed.\n" );
    exit( 1 );
  }
  for( jp=n_cols; jp<2*n_cols; jp++ )coefs[jp]=v2[jp-n_cols];
  
  mkl_free(A);
  mkl_free(v1);
  mkl_free(v2);

  return(0);
}
  
int interpolate_ml(MLDstr *mlDstr, GridMLDstr *rgMLDstr, double *coefs){
  
  int j,jl,m,col,ll,indx;
  int n_cols,ncoef;
  double lat,lon;
  double vn,ve;
  double plm;
  double x;
  int mx_l=MX_L;

  double theta_c=(90.0-fabs(mlDstr->lats[0]))*DTOR;
  double xc=cos(theta_c);

  int ml_lat_indx; 
  double min_lon=1000;

  ncoef=0;
  for( jl=0; jl<mx_l; jl++ )for( m=0; m<=jl; m++)ncoef+=1;  
    
  n_cols=2*ncoef-mx_l;
  
  for( j=0; j<2*ML_NLON; j++ ){
    if( mlDstr->lons[j]<min_lon )min_lon=mlDstr->lons[j];
  }
  
  for( j=0; j<ngrid; j++ ){
    lat=fabs(grid[j].center_lat);
    lon=grid[j].center_lon;    

    if( fabs(lat) < fabs(mlDstr->lats[0])){
      rgMLDstr->vn[j]=0;
      rgMLDstr->ve[j]=0;      
      rgMLDstr->vn_cov[j]=MAX_ML_COV;
      rgMLDstr->ve_cov[j]=MAX_ML_COV;
      continue;
    }
    
    x=cos((90.0-lat)*DTOR);
	
      vn=0;
      ve=0;
      col=0;
      for( jl=0; jl<mx_l; jl++ ) for( m=0; m<=jl; m++){
	  
	  if( m==0 ){
	    
	    plm=calc_plm(ls_1[jl][m],m,x,xc,1)/norms_1[jl][m];
	  
	    if(!isfinite(plm)){
	      fprintf(stderr,"interpolate_ml plm 1 nan %d %d %lf %lf %lf\n",m,jl,x,ls_1[jl][m],norms_1[jl][m]);
	    }
	    vn+=coefs[col]*cos((double)m*lon*DTOR)*plm ;
	    
	    if(isnan(vn)){
	      fprintf(stderr,"interpolate_ml vn nan %d %d %lf %lf %lf\n",m,jl,x,coefs[col],coefs[jl*(jl+1)+2*m+1]);
	    }
	    
	    
	    plm=calc_plm(ls_2[jl][m],m,x,xc,2)/norms_2[jl][m];
	    if(!isfinite(plm)){
	      fprintf(stderr,"interpolate_ml plm2 nan %d %d %lf %lf %lf\n",m,jl,x,ls_2[jl][m],norms_2[jl][m]);
	    }
	    ve+=coefs[col+n_cols]*cos((double)m*lon*DTOR)*plm;
	    
	    col+=1;
	    
	    if(isnan(ve)){
	      fprintf(stderr,"interpolate_ml ve nan %d %d %lf %lf %lf\n",m,jl,x,coefs[jl*(jl+1)+2*m],coefs[jl*(jl+1)+2*m+1]);
	    }
	    
	  }else{	    
	    
	    plm=calc_plm(ls_1[jl][m],m,x,xc,1)/norms_1[jl][m];
	    
	    if(!isfinite(plm)){
	      fprintf(stderr,"interpolate_ml plm 1 nan %d %d %lf %lf %lf\n",m,jl,x,ls_1[jl][m],norms_1[jl][m]);
	    }
	    vn+=coefs[col]*cos((double)m*lon*DTOR)*plm  + coefs[col+1]*sin((double)m*lon*DTOR)*plm;
	    
	    if(isnan(vn)){
	      fprintf(stderr,"interpolate_ml vn nan %d %d %lf %lf %lf\n",m,jl,x,coefs[jl*(jl+1)+2*m],coefs[jl*(jl+1)+2*m+1]);
	    }
	    
	    
	    plm=calc_plm(ls_2[jl][m],m,x,xc,2)/norms_2[jl][m];
	    if(!isfinite(plm)){
	      fprintf(stderr,"interpolate_ml plm2 nan %d %d %lf %lf %lf\n",m,jl,x,ls_2[jl][m],norms_2[jl][m]);
	    }
	    ve+=coefs[col+n_cols]*cos((double)m*lon*DTOR)*plm  + coefs[col+n_cols+1]*sin((double)m*lon*DTOR)*plm;
	    
	    col+=2;
	    
	    if(isnan(ve)){
	      fprintf(stderr,"interpolate_ml ve nan %d %d %lf %lf %lf\n",m,jl,x,coefs[jl*(jl+1)+2*m],coefs[jl*(jl+1)+2*m+1]);
	    }
	  }	
	}
      rgMLDstr->vn[j]=vn;
      rgMLDstr->ve[j]=ve;


      
      ml_lat_indx=(int)((fabs(lat)-ML_LAT_0)/ML_DLAT);
      indx=ml_lat_indx*ML_NLON;
      
      if( lon<mlDstr->lons[indx] )while( mlDstr->lons[indx]>mlDstr->lons[indx-1] )indx++;
      
      if( lon < min_lon ){
	indx=MAX(0,indx-ML_NLON);      
	while((mlDstr->lats[indx]<lat) && (mlDstr->lons[indx+1]>mlDstr->lons[indx]) && indx<mlDstr->npts){
	  indx++;
	}
      }
      
      ll=indx;
      while(( mlDstr->lons[indx] < lon )&&( indx< mlDstr->npts )){
	ll=indx;
	indx++;
      }
      rgMLDstr->vn_cov[j]=mlDstr->vn_cov[ll];
      rgMLDstr->ve_cov[j]=mlDstr->ve_cov[ll];
  }    
  return(0);
}

   
FILE_INFO *file=NULL;
FILE_INFO *select_file(char *radar, time_t time){
  
  char yr_str[5], mo_str[3], dy_str[3];
  char *raid_path;
  char dir_path[PATH_LEN];
  time_t fntime;
  time_t difntime;
  time_t mindif=100000;
  DIR *dp;
  struct dirent *ep;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();

  raid_path=getenv("RAID_PATH");
  struct tm *in_time;
  in_time=gmtime(&time);
  int yr=in_time->tm_year+1900;
  int mo=in_time->tm_mon+1;
  int dy=in_time->tm_mday;

  CNV_TO_STR(yr,yr_str);
  CNV_TO_STR(mo,mo_str);
  CNV_TO_STR(dy,dy_str);
  
  sprintf(file->dir_path,"%s","");
  sprintf(file->fname,"%s","");

  sprintf(dir_path,"%s%s/%s.%s/",raid_path,yr_str,mo_str,dy_str);
  fprintf(stderr,"\n %s%s/%s.%s\n",raid_path,yr_str,mo_str,dy_str);
  fprintf(stderr,"\n dirpath: %s\n",dir_path);
  if((dp=opendir(dir_path))==NULL) {
    fprintf(stderr,"---- COULDN'T OPEN DATA DIRECTORY ----\n");
    closedir(dp);
    return(file);
  }

  while( (ep=readdir(dp)) )
    {
      if( strstr(ep->d_name,".gz") != NULL) continue;
      if( strstr(ep->d_name,".bz2") != NULL) continue;
      if( strstr(ep->d_name,radar) != NULL)
	{
	  fntime=fname_to_time(ep->d_name);
	  difntime=time-fntime;
	  if( difntime>=0 && difntime<mindif )
	  {
	    sprintf(file->dir_path,"%s",dir_path);
	    sprintf(file->fname,"%s",ep->d_name);
	    mindif=difntime;
	  }
	}
    }
  closedir(dp);
  return file;
}



time_t ftime(struct RadarParm *prm){
  int yr,mo,dy,hr,mt,sc;
  yr=prm->time.yr;
  mo=prm->time.mo;
  dy=prm->time.dy;
  hr=prm->time.hr;
  mt=prm->time.mt;
  sc=prm->time.sc;
  return (time_t)TimeYMDHMSToEpoch(yr,mo,dy,hr,mt,sc);
}

void get_pos_ar(struct RadarSite *site, struct RadarParm *prm, struct RadarPos *rdrpos){
  
  double lat,lon,alt=300;
  double magazm,srnge;
  int rn,bm,rsep,frang;
  int s;
  int rxrise,yr;
  int chisham=1;
  rsep=prm->rsep;
  frang=prm->frang;
  rxrise=prm->rxrise;
  yr=prm->time.yr;
  
  for( bm=0; bm<site->maxbeam; bm++) for( rn=0; rn<=prm->nrang; rn++){

      if((double)frang+(double)rsep*(double)rn > (double)max_range)continue;

      /* Calculate magnetic latitude, longitude, and azimuth of range/beam
       * position */
      s=RPosInvMag(bm,rn,yr,site,frang,rsep,rxrise,alt,&lat,&lon,&magazm,&srnge,chisham,0);      

      
      rdrpos->lat[bm][rn]=lat;
      rdrpos->lon[bm][rn]=lon;      
      rdrpos->kazm[bm][rn]=magazm;
    }
}


void map_pos_to_grid(struct RadarPos rdrpos,int frang,int rsep,double bmsep, struct RadarMap *rmap){
  int ib,ir,ig,ig_found,icnt;
  int bsteps,ibstep;
  int count;
  double mlat,mlon,kazm,mlat_c,mlon_c;
  double bm_delta;

  bsteps=MAX_COUNT/2;

  for( ib=0; ib<MAX_BEAM; ib++)for( ir=0; ir<MAX_GATES; ir++){
      mlat_c=rdrpos.lat[ib][ir];
      mlon_c=rdrpos.lon[ib][ir];
      kazm=rdrpos.kazm[ib][ir];

      if((double)frang+(double)rsep*(double)ir > (double)max_range)continue;
      
      count=0;
      bm_delta=((double)frang+(double)rsep*(double)ir)*bmsep/(double)bsteps;
      /* fprintf(stderr,"\n map_pos_pt_grid: %d %d %d %f %f %f\n",ir,frang,rsep,(double)frang+(double)rsep*(double)ir,bmsep,bm_delta); */
      for(ibstep=0; ibstep<bsteps; ibstep++){
	sub_sphcal(mlon_c, mlat_c, kazm-90, (double)ibstep*bm_delta, &mlon, &mlat);
	for(ig=0; ig<ngrid; ig++)if(cell_in_out(mlat,mlon,grid[ig])){
	    for( icnt=0; icnt<=count; icnt++){
	      ig_found=0;
	      if( ig == rmap->cell[ib][ir][icnt] ) break;
	      ig_found=1;
	    }
	    if( ig_found == 1 ){
	      rmap->cell[ib][ir][count]=ig;
	      count++;
	    }
	    
      }
    }
      
      for(ibstep=0; ibstep<bsteps; ibstep++){
	sub_sphcal(mlon_c, mlat_c, kazm+90, (double)ibstep*bm_delta, &mlon, &mlat);
	for(ig=0; ig<ngrid; ig++)if(cell_in_out(mlat,mlon,grid[ig])){
	    for( icnt=0; icnt<=count; icnt++){
	      ig_found=0;
	      if( ig == rmap->cell[ib][ir][icnt] ) break;
	      ig_found=1;
	    }
	    if( ig_found == 1 ){
	      rmap->cell[ib][ir][count]=ig;
	      count++;
	    }
	  }
      }
      rmap->cell_count[ib][ir]=count;      
    }
}
  
double grid_lat(int ig){
  return((grid[ig].lat[0]+grid[ig].lat[1]+grid[ig].lat[2]+grid[ig].lat[3])/4);
}

double grid_lon(int ig){
  return((grid[ig].lon[0]+grid[ig].lon[1]+grid[ig].lon[2]+grid[ig].lon[3])/4);
}


struct RadarParm *prm;
struct FitData *fit;
struct FitIndex *inx;

struct RadarNetwork *network;  
struct Radar *radar;

struct RadarPos *rpos;


extern int pinv(int m,int n,double* a,double* inv);  
MLDstr mlDstr;
GridMLDstr rgMLDstr;

double *coefs;


int main(int argc, char *argv[]){

  /* read instruction file
     file should have date, time, lat and lon of area corners, radars to contribute
     filtering instructions
  */

  FILE *fp;
  FILE *vf;
  FILE *mlf;
  FILE *rptrs[NRAD];
  FILE *dataout;
  FILE *losout;
  FILE *divfile;
  char file_name[128];
  int nbms;
  struct tm *t_start;
  struct tm *t_end; 
  time_t fit_time, esec;
  char *envstr;
  struct RadarSite *site=NULL;
  struct RadarMap *rmap;
  int *grid_data_count;

  double *coef;
  double *coefTrans;
  double *data;
  double *A;
  double *INV_ar;
  double *CDI;
  double *GD;
  double alpha=1;
  double beta=0;
  double dot,row_mag0,row_mag1;
  
  int stat;
  double *solution_l;
  double ml_time;

  double min_lon=0;
  double max_lon=360;

  int fpos;
  int yr,mo,dy,hr,mt;
  int j,jr,jc,jeqn;
  int neqn,ndata;
  double t1,t2;
  double sc;
  
  int yr_st,mo_st,dy_st,hr_st,mt_st;
  int fit_yr,fit_mo,fit_dy,fit_hr,fit_mt;
  double sc_st,fit_sc;
  int jj,jjc,jjr,count;
  int jcmn,jcmx,jrmn,jrmx;
  int frang,rsep;
  int *rsepl;
  int *frangl;
  int *nrangl;
  double bmsep;
  double kazm;
  double avg_val;
  double cor_mag,cor_azm;
  int grid_cell;
  int try_count;
  int count_lim=20;
  int s,jrmin,jrmax;

  double hold;
  double **filt_ar;
  double **filt_ar1;
  double **filt_err;
  double **filt_cnt;
  double med_filt_ar[9],median,tmp,tp;
  int **good_ar;

  int ncoef, n_cols;
  int jl,m;
  int mx_l=MX_L;
  ncoef=0;
  for( jl=0; jl<mx_l; jl++ )for( m=0; m<=jl; m++)ncoef+=1;  
    
  n_cols=2*ncoef-mx_l;

  coefs=(double *)calloc(2*n_cols,sizeof(double));
  
  mlDstr.lats=NULL;
  mlDstr.lons=NULL;
  mlDstr.ve=NULL;
  mlDstr.vn=NULL;
  mlDstr.ve_cov=NULL;
  mlDstr.vn_cov=NULL;
  mlDstr.vmag=NULL;
  mlDstr.vaz=NULL;
  
  if (argc <= 1){
    fprintf(stderr,"********NO INSTRUCTION FILE GIVEN********\n");
    exit(-1);
  }
  strcpy(file_name,argv[1]);
  if ((fp=fopen(file_name,"r")) == 0) {
    fprintf(stderr,"******FIT INSTRUCTION FILE %s NOT FOUND******\n",argv[1]);
    exit(-1);
  }
  /* parse file */
  parse_instructions(fp);
  fclose(fp);


  t_start=parse_date_str(start_time);
  t_end=parse_date_str(end_time);

  fit_time=mktime(t_start);
  esec=mktime(t_end);  

  
  TimeEpochToYMDHMS(fit_time,&yr,&mo,&dy,&hr,&mt,&sc);
  sc=0;
  fprintf(stderr,"%d %d %d %d %d %f\n",yr,mo,dy,hr,mt,sc);
  IGRF_SetDateTime(yr,mo,dy,hr,mt,(int)sc);
  AACGM_v2_SetDateTime(yr,mo,dy,hr,mt,(int)sc);  
  
  min_lon_ng=99999.9;
  max_lon_ng=-99999.9;
  for( j=0; j<nbp; j++ ){
    if( bp[j].lon<min_lon_ng )min_lon_ng=bp[j].lon;
    if( bp[j].lon>max_lon_ng )max_lon_ng=bp[j].lon;
  }
  

  /*   create grid   */

  int nlat=(max_lat-min_lat)/dlat;
  start_lon=(double *)calloc(nlat,sizeof(double));

  bool random_start=FALSE;
  set_start_lon(random_start, nlat, min_lat, dlat, start_lon);  

  ngrid=get_grid_size(min_lat, max_lat, nbp, bp);

  grid=(CELL *)calloc(ngrid,sizeof(CELL));
  make_grid(min_lat, max_lat, nbp, bp, hemisphere, grid);
  
  if( strcmp(hemisphere,"south")==0 ){
    hold=max_lat;
    max_lat=-min_lat;
    min_lat=-hold;
  }

  /* write grid to file and determine cell neighbors */
  
  grd_file=fopen("grid.dat","w");
  
  fprintf(grd_file,"min_lat: %f\n",min_lat);
  fprintf(grd_file,"max_lat: %f\n",max_lat);
  fprintf(grd_file,"min_lon: %f\n",min_lon);
  fprintf(grd_file,"max_lon: %f\n",max_lon);
  fprintf(grd_file,"ngrid: %d\n",ngrid);
  
  for( j=0; j<ngrid; j++ ){
    fprintf(grd_file,"%d %f  %f\n",j,grid[j].center_lat,grid[j].center_lon);
    fprintf(grd_file,"%f  %f  %f  %f\n",grid[j].lat[0],grid[j].lat[1],grid[j].lat[2],grid[j].lat[3]);
    fprintf(grd_file,"%f  %f  %f  %f\n",grid[j].lon[0],grid[j].lon[1],grid[j].lon[2],grid[j].lon[3]);    
  }
  fclose(grd_file);

  grid_data_count=malloc(ngrid*sizeof(int));
  
  rgMLDstr.vn=(double *)calloc(ngrid,sizeof(double));
  rgMLDstr.ve=(double *)calloc(ngrid,sizeof(double));
  rgMLDstr.vn_cov=(double *)calloc(ngrid,sizeof(double));
  rgMLDstr.ve_cov=(double *)calloc(ngrid,sizeof(double));
  
  neighbors=(NEIGHBOR *)calloc(ngrid,sizeof(struct Neighbor));
  find_neighbors(hemisphere, grid, ngrid, neighbors);
  
  /* create array of coefficients for calculating divergence */
  divfile=fopen("divergence.dat","w");
  make_div_ar();

  for( jr=0; jr<ngrid-nedge; jr++){
    fprintf(divfile,"%d ",jr);
    for( jc=0; jc<2*ngrid; jc++)if(div_array[jr].coef[jc] !=0 )fprintf(divfile,"   %d %f\n",jc,div_array[jr].coef[jc]);
    fprintf(divfile,"\n");
  }
  fclose(divfile);
  
  envstr=getenv("SD_RADAR");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_RADAR' must be defined.\n");
    exit(-1);
  }

  fp=fopen(envstr,"r");

  if (fp==NULL) {
    fprintf(stderr,"Could not locate radar information file.\n");
    exit(-1);
  }

  network=RadarLoad(fp);
  fclose(fp); 
  if (network==NULL) {
    fprintf(stderr,"Failed to read radar information.\n");
    exit(-1);
  }

  envstr=getenv("SD_HDWPATH");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_HDWPATH' must be defined.\n");
    exit(-1);
  }

  RadarLoadHardware(envstr,network);

  /* find fit files and read into model-array and data vector
     requires determining k-vectors
   */

  if(file)free(file);
  file=calloc(1,sizeof(FILE_INFO));

  int res;
  nrad--;
  for( jj=0; jj<nrad; jj++ )
    {
      rptrs[jj]=NULL;
      file=select_file(radar_list[jj],fit_time);
      if ( (res=strcmp(file->fname,"")) !=0){
	fprintf(stderr,"%s%s\n",file->dir_path,file->fname);
	sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	if(rptrs[jj])fclose(rptrs[jj]);
	if((rptrs[jj]=fopen(file_name,"r"))== NULL)fprintf(stderr,"could not open file: %s",file_name);
      }
    }
  fprintf(stderr,"first files open\n");
  
  if ((mlf=fopen(ml_file,"r")) == 0) {
    fprintf(stderr,"******ML FILE %s NOT FOUND******\n",ml_file);
    exit(-1);
  }

  fpos=read_ml_record(mlf, &mlDstr, hemisphere);  

  fprintf(stderr,"back from ml read npts %d\n",mlDstr.npts);
  double theta_c=(90.0-mlDstr.lats[0])*DTOR;
  double xc=cos(theta_c);
  
  ls_1=find_zeros_1(mx_l,xc);
  ls_2=find_zeros_2(mx_l,xc);

  fprintf(stderr,"have ls\n");

  norms_1=malloc(mx_l*sizeof(double*));
  for( jl=0; jl<mx_l; jl++){
    norms_1[jl]=(double *)calloc(mx_l,sizeof(double));
  }

  for( m=0; m<=mx_l; m++){
    for( jl=0; jl<mx_l; jl++ ){
      if( jl<m )continue;
      norms_1[jl][m]=plm_norm(ls_1[jl][m],m,xc,(int) 1);
      /* fprintf(stderr,"norms 1: %d %d %lf %le\n",m,jl,ls_1[jl][m],norms_1[jl][m]); */
    }
  }    
  fprintf(stderr,"have norms 1\n");

  norms_2=malloc(mx_l*sizeof(double*));
  for( jl=0; jl<mx_l; jl++){
    norms_2[jl]=(double *)calloc(mx_l,sizeof(double));
  }

  for( m=0; m<=mx_l; m++){
    for( jl=0; jl<mx_l; jl++ ){
      if( jl<m )continue;
      norms_2[jl][m]=plm_norm(ls_2[jl][m],m,xc,(int) 2);	
      /* fprintf(stderr,"norms 2: %d %d %lf %le\n",m,jl,ls_2[jl][m],norms_2[jl][m]); */
    }
  }

  fprintf(stderr,"have norms 2\n");
  
  
  s=fit_ml(&mlDstr, coefs);
  
  fprintf(stderr,"read ml_record\n");
  
  vf=fopen("vel_out","w");
  losout=fopen("los_out","w");
  /* fp=fopen("coef_array","w"); */

  nrangl=calloc(nrad, sizeof(int));
  rsepl=calloc(nrad, sizeof(int));
  frangl=calloc(nrad, sizeof(int));
  
  prm=RadarParmMake();
  fit=FitMake();
  rpos=calloc(nrad, sizeof(struct RadarPos));
  rmap=calloc(nrad, sizeof(struct RadarMap));
  
  
  while( fit_time<esec )
    {
      TimeEpochToYMDHMS(fit_time,&fit_yr,&fit_mo,&fit_dy,&fit_hr,&fit_mt,&fit_sc);      
      TimeEpochToYMDHMS(fit_time-(time_t)avg_ival/2,&yr_st,&mo_st,&dy_st,&hr_st,&mt_st,&sc_st);
      printf("%d %d %d %d %d %f\n",yr_st,mo_st,dy_st,hr_st,mt_st,sc_st);
      fprintf(stderr,"fit_time time: %d %d %d %d %d %f\n",fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,fit_sc);
      fprintf(stderr,"Interval start time: %d %d %d %d %d %f\n",yr_st,mo_st,dy_st,hr_st,mt_st,sc_st);
      ndata=0;
      for(jj=0; jj<ngrid; jj++)grid_data_count[jj]=0;
      for( jj=0; jj<nrad; jj++){


	if( rptrs[jj]!=NULL ){
	  s=fseek(rptrs[jj],0L,SEEK_SET);
	  if( s==0 )s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	  if( s==0 )s=FitFread(rptrs[jj],prm,fit);
	  if( s==0 )fprintf(stderr,"%s file time: %d %d %d %d %d %d %d %ld\n",radar_list[jj],prm->time.yr,
			    prm->time.mo,prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,prm->stid,
			    ftime(prm)-fit_time);
	  if( s!=0 ){
	    fprintf(stderr,"1. File does not contain the requested interval. %d:%d\n",hr_st,mt_st);
	    if(rptrs[jj]){
	      fprintf(stderr,"*******Closing %s file\n",radar_list[jj]);
	      fclose(rptrs[jj]);
	      rptrs[jj]=NULL;
	    }
	  }
	}

	if( rptrs[jj]==NULL ){
	  file=select_file(radar_list[jj],fit_time-(time_t)avg_ival/2);
	  if ( (res=strcmp(file->fname,"")) ==0)continue;
	  /* if( file->fname==NULL ) continue; */
	  sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	  fprintf(stderr,"*******Opening file %s\n",file_name); 
	  rptrs[jj]=fopen(file_name,"r");
	  if( rptrs[jj]!=NULL ){
	    s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	    if( s==0 )s=FitFread(rptrs[jj],prm,fit);
	    /* TimeEpochToYMDHMS(fit_time+(time_t)avg_ival,&yr,&mo,&dy,&hr,&mt,&sc); */
	  }
	}
	  
	radar=RadarGetRadar(network,prm->stid);
	site=RadarYMDHMSGetSite(radar,yr_st,mo_st,dy_st,hr_st,mt_st,(int) 0);

	fprintf(stderr,"past radar and site\n");

	if (site==NULL) {fprintf(stderr,"NULL site\n"); continue;}

	nbms=(int)site->maxbeam;
	filt_ar=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_ar[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_ar1=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_ar1[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_err=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_err[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_cnt=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_cnt[jr]=(double*)calloc(prm->nrang,sizeof(double));
	good_ar=malloc(nbms*sizeof(int*));
	for( jr=0; jr<nbms; jr++)good_ar[jr]=(int*)calloc(prm->nrang,sizeof(int));
	nrangl[jj]=prm->nrang;
	
	fprintf(stderr,"allocated filter arrays\n");
	
	try_count=0;
	if( rptrs[jj]!=NULL ){
	  s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	  while( ((ftime(prm)-fit_time)<=(long)avg_ival/2) && (try_count<count_lim)){
	    if( s!=0 ){
	      try_count++;
	      fprintf(stderr,"2. File does not contain the requested interval. %d:%d\n",hr_st,mt_st);
	      fprintf(stderr,"%s  %ld\n",radar_list[jj],fit_time+(time_t)avg_ival);
	      file=select_file(radar_list[jj],fit_time+(time_t)avg_ival);
	      sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	      if(rptrs[jj]){
		fprintf(stderr,"*******Closing %s file\n",radar_list[jj]);
		fclose(rptrs[jj]);
	      }
	      rptrs[jj]=NULL;
	      /* if( file->fname == NULL )continue; */
	      if ( (res=strcmp(file->fname,"")) ==0)continue;
	      fprintf(stderr,"*******Opening file %s\n",file_name); 
	      if( (rptrs[jj]=fopen(file_name,"r")) ==NULL )continue;
	      /*	      rewind(rptrs[jj]);*/
	      /* TimeEpochToYMDHMS(fit_time+(time_t)avg_ival,&yr,&mo,&dy,&hr,&mt,&sc); */
	      if((s=FitFseek(rptrs[jj],fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,fit_sc,NULL,inx))!=0 ||
		 (s=FitFread(rptrs[jj],prm,fit))!=0){
		fprintf(stderr,"problem reading %s\n",file_name);
		fclose(rptrs[jj]);
		rptrs[jj]=NULL;
		break;
	      }
	      fprintf(stderr,"main: %s\n",file->fname);
	    }


	    if( s==0 ){

	      rsep=prm->rsep;
	      frang=prm->frang;
	      if((rsep != rsepl[jj]) || (frang != frangl[jj])){
		fprintf(stderr,"recalculating position array %s %d %d  %d  %d %d\n",
			radar_list[jj],rsep,rsepl[jj],frang,frangl[jj],prm->nrang);
		radar=RadarGetRadar(network,prm->stid);
		site=RadarYMDHMSGetSite(radar,fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,(int) fit_sc);
		if (site==NULL){fprintf(stderr,"site was null\n"); continue;}
		bmsep=(double)site->bmsep*(double)PI/180.0;
		get_pos_ar(site,prm,&rpos[jj]);
		map_pos_to_grid(rpos[jj],frang,rsep,bmsep,&rmap[jj]);
		rsepl[jj]=rsep;
		frangl[jj]=frang;

	      }
	      if(prm->nrang != nrangl[jj]){
		fprintf(stderr,"%d station: %s  nrang mismatch old: %d   new: %d...reallocating\n",jj,radar_list[jj],nrangl[jj],prm->nrang);

		for( jr=0; jr<nbms; jr++)if(filt_ar[jr] != NULL)free(filt_ar[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_ar1[jr] != NULL)free(filt_ar1[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_err[jr] != NULL)free(filt_err[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_cnt[jr] != NULL)free(filt_cnt[jr]); 
		for( jr=0; jr<nbms; jr++)if(good_ar[jr] != NULL)free(good_ar[jr]);


		for( jr=0; jr<nbms; jr++)filt_ar[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_ar1[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_err[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_cnt[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)good_ar[jr]=(int*)calloc(prm->nrang,sizeof(int));
		
		if(prm->nrang == 0)for( jr=0; jr<nbms; jr++){
		    filt_ar[jr]=NULL;
		    filt_ar1[jr]=NULL;
		    filt_err[jr]=NULL;
		    filt_cnt[jr]=NULL;
		    good_ar[jr]=NULL;
		  }		  
		nrangl[jj]=prm->nrang;
	      }
	      jrmin=(min_range-frang)/rsep;
	      jrmax=(max_range-frang)/rsep;
	      if( labs(ftime(prm)-fit_time) < (double)avg_ival/2){
		yr=prm->time.yr;
		mo=prm->time.mo;
		dy=prm->time.dy;
		hr=prm->time.hr;
		mt=prm->time.mt;
		sc=(double)prm->time.sc;
		/* fprintf(stderr,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc); */
		/* fprintf(stderr,"fit %d %d %d %d %d %d\n",fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,(int)fit_sc); */
		for( jr=jrmin; jr<jrmax; jr++ ){
		  if(fit->rng[jr].qflg == 1 && fit->rng[jr].gsct == 0 && 
		     (fit->rng[jr-1].gsct == 0 && fit->rng[jr+1].gsct == 0) &&
		     fabs(fit->rng[jr].v_err) < MAX_V_ERR &&
		     fabs(fit->rng[jr].v) < MAX_V && fabs(fit->rng[jr].v) > MIN_V){
		    /* fprintf(stderr,"v %f\n",fit->rng[jr].v); */
		    filt_ar[prm->bmnum][jr]+=fit->rng[jr].v;
		    filt_err[prm->bmnum][jr]+=ERR_SCALE*fit->rng[jr].v_err*fit->rng[jr].v_err;
		    filt_cnt[prm->bmnum][jr]+=1.;		    
		  }
		}
	      }
	    }
	    s=FitFread(rptrs[jj],prm,fit);
	  }
	}

	for( jc=0; jc<site->maxbeam; jc++ )for( jr=1; jr<prm->nrang; jr++){
	    jcmn=MAX(jc-1,0);
	    jcmx=MIN(jc+1,site->maxbeam-1);
	    jrmn=MAX(jr-1,0);
	    jrmx=MIN(jr+1,prm->nrang-1);
	    count=0;
	    avg_val=0;
	    filt_ar1[jc][jr]=0;
	    good_ar[jc][jr]=0;
	    if( filt_cnt[jc][jr]==0 ) continue;
	    
	    for(jjc=jcmn; jjc<=jcmx; jjc++)for(jjr=jrmn; jjr<=jrmx; jjr++){
		if( filt_cnt[jjc][jjr]!=0 ){
		  med_filt_ar[count]=filt_ar[jjc][jjr]/filt_cnt[jjc][jjr];
		  count++;
		  avg_val+=(filt_ar[jjc][jjr]/filt_cnt[jjc][jjr]);
		}
	      }

	    if( count<MIN_COUNT)continue;
	    
	    for( jjc=0; jjc<count; jjc++ )for( jjr=jjc+1; jjr<count; jjr++ ){
		tp=med_filt_ar[jjc];
		if( tp>med_filt_ar[jjr] ){
		  tmp=med_filt_ar[jjr];
		  med_filt_ar[jjr]=tp;
		  med_filt_ar[jjc]=tmp;
		  tp=tmp;
		}
	      }
	    
	    if( count%2 ){
	      median=med_filt_ar[(int)(count/2)];
	    }else{
	      median=(med_filt_ar[(int)(count/2)]+med_filt_ar[(int)(count/2)+1])/2;
	    }
		
	    /* for(jjc=0; jjc<count; jjc++) */
	    /*   fprintf(stderr,"%d med_filt_ar %lf \n",jjc,med_filt_ar[jjc]); */
	    
	    /* fprintf(stderr,"The median was: %lf\n",median); */

	    
	    avg_val/=(double)count;
	    
	    good_ar[jc][jr]=1;

	    /* filt_ar1[jc][jr]=avg_val; */
	    filt_ar1[jc][jr]=median;
	    if(count<MIN_COUNT){
	      good_ar[jc][jr]=0;
	      filt_ar1[jc][jr]=0;
	      continue;
	    }
	  }
	for( jc=0; jc<nbms; jc++ )for( jr=1; jr<prm->nrang; jr++){
	    for( jjc=0; jjc<rmap[jj].cell_count[jc][jr]; jjc++){
	      if ((grid_cell=rmap[jj].cell[jc][jr][jjc]) != -1 && good_ar[jc][jr] == 1 &&
		  filt_cnt[jc][jr]>0 && fabs(rpos[jj].lat[jc][jr]) < MAX_DATA_LAT){	
		grid_data_count[grid_cell]++;
		    
		kazm=rpos[jj].kazm[jc][jr];

		if( strcmp(hemisphere,"south")==0 ){
		  kazm+=180;
		}
	      
		/* printf("%d %d %d %f %f %f %f %f %f %f\n",prm->stid,jc,jr,rpos[jj].lat[jc][jr], */
		/* 	     rpos[jj].lon[jc][jr],kazm,filt_ar1[jc][jr]/filt_cnt[jc][jr], */
		/* 	     grid_lat(grid_cell),grid_lon(grid_cell)); */
		kazm_array=realloc(kazm_array,(ndata+1)*sizeof(struct mod_array));
		kazm_array[ndata].coef=calloc(2*ngrid,sizeof(double));
		los_data=realloc(los_data,(ndata+1)*sizeof(double));
		los_kazm=realloc(los_kazm,(ndata+1)*sizeof(double));
		los_lats=realloc(los_lats,(ndata+1)*sizeof(double));
		los_lons=realloc(los_lons,(ndata+1)*sizeof(double));
		los_err=realloc(los_err,(ndata+1)*sizeof(double));
		kazm_array[ndata].coef[grid_cell]=sind(kazm);
		kazm_array[ndata].coef[grid_cell+ngrid]=cosd(kazm);

		if( strcmp(hemisphere,"south")==0 ){
		  kazm_array[ndata].coef[grid_cell]=-sind(kazm);
		}

		
		los_lats[ndata]=grid_lat(grid_cell);
		los_lons[ndata]=grid_lon(grid_cell);
		los_kazm[ndata]=kazm;
		los_data[ndata]=(-1.)*filt_ar1[jc][jr]; /* +cor_mag*cosd(kazm-cor_azm); */
		los_err[ndata]=filt_err[jc][jr]/filt_cnt[jc][jr];
		/* fprintf(stderr,"%d %d %f %f %f %f %f\n",ndata,prm->stid,los_lats[ndata],los_lons[ndata],los_data[ndata],los_err[ndata],cor_mag*cosd(kazm-cor_azm)); */
		
		ndata++;
	      }
	    }
	  }
	for( jr=0; jr<nbms; jr++)if(filt_ar[jr] != NULL)free(filt_ar[jr]); 
	free(filt_ar);
	for( jr=0; jr<nbms; jr++)if(filt_ar1[jr] != NULL)free(filt_ar1[jr]); 
	free(filt_ar1);
	for( jr=0; jr<nbms; jr++)if(filt_err[jr] != NULL)free(filt_err[jr]); 
	free(filt_err);
	for( jr=0; jr<nbms; jr++)if(filt_cnt[jr] != NULL)free(filt_cnt[jr]); 
	free(filt_cnt);
	for( jr=0; jr<nbms; jr++)if(good_ar[jr] != NULL)free(good_ar[jr]); 
	free(good_ar);
      }
    
    /* synchronize ML model with grid data */
    
    ml_time=TimeYMDHMSToEpoch(mlDstr.yr,mlDstr.mo,mlDstr.dy,mlDstr.hr,
			      mlDstr.mt,(double)mlDstr.sc);
    while( ml_time<fit_time ){
      if((fpos=read_ml_record(mlf, &mlDstr, hemisphere))==-1)break;
      s=fit_ml(&mlDstr, coefs);
      
      ml_time=TimeYMDHMSToEpoch(mlDstr.yr,mlDstr.mo,mlDstr.dy,
				mlDstr.hr,mlDstr.mt,(double)mlDstr.sc);
      
    }

    fprintf(stderr,"grid: %d %d %d %d %d %lf\n",yr,mo,dy,hr,mt,sc);
    fprintf(stderr,"ML: %d %d %d %d %d %d\n",mlDstr.yr,mlDstr.mo,mlDstr.dy,mlDstr.hr,mlDstr.mt,mlDstr.sc);
    
    /* Interpolate ML model to grid */
    /* if( regrid_ml(&mlDstr,&rgMLDstr)!=-1){ */
    if( interpolate_ml(&mlDstr,&rgMLDstr,coefs)!=-1){

      fprintf(stdout,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc);
      fprintf(stdout,"%d\n",ngrid);
      /* for( j=0; j<ngrid; j++ ){ */
	
      /* 	fprintf(stdout,"%d %5.2f %5.2f %8.3f %8.3f %8.3f %8.3f\n",j,grid[j].center_lat,grid[j].center_lon,rgMLDstr.ve[j],rgMLDstr.vn[j],rgMLDstr.ve_cov[j],rgMLDstr.vn_cov[j]); */
      /* 	fflush(stdout); */
      /* } */
      
      neqn=3*ngrid-nedge+ndata;

      fprintf(stderr,"NGRID=%d ",ngrid);
      fprintf(stderr,"NEDGE=%d ",nedge);
      fprintf(stderr,"NDATA=%d ",ndata);
      fprintf(stderr,"NEQN=%d\n",neqn);
      coef=(double *)mkl_calloc((size_t)neqn*(size_t)ngrid*2,sizeof(double), 64);
      coefTrans=(double *)mkl_calloc((size_t)ngrid*2*(size_t)neqn,sizeof(double), 64);
      CDI=(double *)mkl_calloc((size_t)neqn*(size_t)neqn,sizeof(double), 64);
      data=(double *)mkl_calloc((size_t)neqn,sizeof(double), 64);
      
      fprintf(stderr,"***ALLOCATED***");
      
      
      jeqn=0;
      /* Edges are low-latitude boundary points, which are the first nedge points of the grid */
      for( jr=0; jr<ngrid; jr++){
	if(isnan(rgMLDstr.ve[jr])||isnan(rgMLDstr.vn[jr])){
	  fprintf(stderr,"NAN: %d %f %f\n",jr,rgMLDstr.ve[jr],rgMLDstr.vn[jr]);
	}
	data[jeqn]=rgMLDstr.ve[jr];
	coef[jeqn*ngrid*2+jr]=1;
	CDI[jeqn*neqn+jeqn]=1/MAX(rgMLDstr.ve_cov[jr],MIN_ML_ERR);
	
	jeqn++;

	data[jeqn]=rgMLDstr.vn[jr];
	coef[jeqn*ngrid*2+jr+ngrid]=1;
	CDI[jeqn*neqn+jeqn]=1/(MAX(rgMLDstr.vn_cov[jr],MIN_ML_ERR));
	jeqn++;
      }
      fprintf(stderr,"MODEL SET  %d***",jeqn);
      
      
      for( jr=nedge; jr<ngrid; jr++){
	for( jc=0; jc<2*ngrid; jc++)coef[jeqn*2*ngrid+jc]=div_array[jr].coef[jc];
	data[jeqn]=0;
	CDI[jeqn*neqn+jeqn]=1/(double)DIV_ERR;;
	
	jeqn++;
      }
      fprintf(stderr,"DIVERGENCE SET  %d***",jeqn);
      
      for( jr=0; jr<ndata; jr++){
	
	data[jeqn]=los_data[jr];
	
	for( jc=0; jc<2*ngrid; jc++)coef[jeqn*2*ngrid+jc]=kazm_array[jr].coef[jc];
	CDI[jeqn*neqn+jeqn]=1/(MAX(los_err[jr],MIN_ERR));
	jeqn++;
      }
      fprintf(stderr,"DATA SET***jeqn %d****\n",jeqn);

      t1=dsecnd();
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,(MKL_INT)(2*ngrid),(MKL_INT)neqn,(MKL_INT)neqn,(double)1.0,coef,(MKL_INT)2*ngrid,CDI,(MKL_INT)neqn,(double)0.,coefTrans,(MKL_INT)neqn);

      t2=dsecnd();
      
      fprintf(stderr,"time for %d by %d dgemm: %f\n",neqn,ngrid*2,(t2-t1));
      fprintf(stderr,"back from gt=g^T*CD^-1\n");
      
      GD=(double *)mkl_calloc(2*(size_t)ngrid, sizeof(double), 64);
      cblas_dgemv(CblasRowMajor,CblasNoTrans,(MKL_INT)2*ngrid,(MKL_INT)neqn,(double)1.,coefTrans,neqn,data,1,0.,GD,(MKL_INT)1);
      
      A=(double *)mkl_calloc(2*(size_t)ngrid*2*(size_t)ngrid,sizeof(double),64);
      
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(MKL_INT)(2*ngrid),(MKL_INT)(2*ngrid),(MKL_INT)neqn,(double)1.0,coefTrans,(MKL_INT)neqn,coef,(MKL_INT)(2*ngrid),0.0,A,(MKL_INT)(2*ngrid));
    

      if( WRITE_COEF ){
    
	dataout=fopen("data_out","w");
	for( jr=0; jr<2*ngrid; jr++ ){
	  fprintf(dataout,"%lf\n",GD[jr]);
	}
	for( jr=0; jr<2*ngrid; jr++ )for(jc=0; jc<2*ngrid; jc++ )if( A[jr*2*ngrid+jc]!=0 )fprintf(dataout,"coefficient %d %d %lf\n",jr,jc,A[jr*2*ngrid+jc]); 
	/* for( jr=0; jr<2*ngrid; jr++ )for(jc=0; jc<ngrid; jc++ )if( aa[jr][jc]!=0 )fprintf(dataout,"coefficient %d %d %lf\n",jr,jc,aa[jr][jc]);  */
	fclose(dataout);
      }
      
      for( jr=0; jr<2*ngrid; jr++)for( jc=0; jc<2*ngrid; jc++)if( !isfinite(A[jr*2*ngrid+jc]) )fprintf(stderr,"coefficient %d %d %lf\n",jr,jc,A[jr*2*ngrid+jc]); 
      
      solution_l=calloc(2*ngrid,sizeof(double));
      
      t1=dsecnd();
      
      INV_ar=(double *)mkl_calloc(2*(size_t)ngrid*2*(size_t)ngrid,sizeof(double),64);
      stat=pinv((int)(2*ngrid),(int)(2*ngrid),A,INV_ar);
      if( stat<0 ){
	fprintf(stderr,"pinv_returned %d\n",stat);

	dataout=fopen("data_out","w");
	for( jr=0; jr<2*ngrid; jr++ ){
	  fprintf(dataout,"%lf\n",GD[jr]);
	}
	for( jr=0; jr<2*ngrid-1; jr++ ){
	  /* fprintf(stderr,"row %d\n",jr); */
	  for( jjr=jr+1; jjr<2*ngrid-1; jjr++ ){
	    dot=0;
	    row_mag0=0;
	    row_mag1=0;
	    for(jc=0; jc<2*ngrid; jc++ ){
	      /* if( A[jr*2*ngrid+jc]!=0 )fprintf(dataout,"coefficient %d %d %lf\n",jr,jc,A[jr*2*ngrid+jc]); */
	      row_mag0+=A[jr*2*ngrid+jc]*A[jr*2*ngrid+jc];
	      row_mag1+=A[jjr*2*ngrid+jc]*A[jjr*2*ngrid+jc];
	      dot +=A[jr*2*ngrid+jc]*A[jjr*2*ngrid+jc];
	    }
	    if( dot/(sqrt(row_mag0)*sqrt(row_mag1)) > 0.9 ){fprintf(stderr,"nearly parallel ");
	      fprintf(stderr,"rows %d %d magnitude: %lf %le %le\n",jr,jjr,dot/(sqrt(row_mag0)*sqrt(row_mag1)),row_mag0,row_mag1);
	    }
	  }
	}
	/* for( jr=0; jr<2*ngrid; jr++ )for(jc=0; jc<ngrid; jc++ )if( aa[jr][jc]!=0 )fprintf(dataout,"coefficient %d %d %lf\n",jr,jc,aa[jr][jc]);  */
	fclose(dataout);

	exit(-1);
      }
      cblas_dgemv(CblasRowMajor,CblasNoTrans,2*ngrid,2*ngrid,alpha,INV_ar,2*ngrid,GD,(const MKL_INT)1,beta,solution_l,(const MKL_INT)1);
      
      t2=dsecnd();
      
      fprintf(stderr,"time for %d by %d pinv solution: %f\n",neqn,ngrid*2,(t2-t1));
      
      fprintf(stderr,"grid: %d %d %d %d %d %lf\n",yr,mo,dy,hr,mt,sc);
      fprintf(stderr,"ML: %d %d %d %d %d %d\n",mlDstr.yr,mlDstr.mo,mlDstr.dy,mlDstr.hr,mlDstr.mt,mlDstr.sc);    
      TimeEpochToYMDHMS(fit_time+(time_t)avg_ival/2,&yr,&mo,&dy,&hr,&mt,&sc);
      fprintf(stderr,"fit_time: %d %d %d %d %d %lf\n",yr,mo,dy,hr,mt,sc);

      fprintf(vf,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc);
      fprintf(vf,"%f %f %f %f %f %f\n",mlDstr.Bx,mlDstr.By,mlDstr.Bz,mlDstr.v_sw,mlDstr.Au,mlDstr.Al);
      fprintf(vf,"%d\n",ngrid);
      for( jc=0; jc<ngrid; jc++)
	fprintf(vf,"%f %f %f %f %f %f %d\n",grid[jc].center_lat,grid[jc].center_lon,solution_l[jc],solution_l[jc+ngrid],rgMLDstr.ve[jc],rgMLDstr.vn[jc],grid_data_count[jc]);
    
      fflush(vf);

      fprintf(losout,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc);
      fprintf(losout,"%f %f %f %f %f %f\n",mlDstr.Bx,mlDstr.By,mlDstr.Bz,mlDstr.v_sw,mlDstr.Au,mlDstr.Al);
      fprintf(losout,"%d\n",ndata);
      for( jc=0; jc<ndata; jc++)
	if( strcmp(hemisphere,"south")==0 ){
	  los_kazm[jc]-=180;
	}

	fprintf(losout,"%f %f %f %f\n",los_lats[jc],los_lons[jc],los_data[jc],los_kazm[jc]);
      
      
      free(solution_l);
          
      mkl_free(coef);
      mkl_free(coefTrans);
      mkl_free(CDI);
      mkl_free(GD);
      mkl_free(data);
      mkl_free(A);
      mkl_free(INV_ar);

      for( jr=0; jr<ndata; jr++)free(kazm_array[jr].coef);
    }
    
    fit_time+=(time_t)avg_ival;
    
  }
  free(div_array);
  free(grid);
  free(grid_data_count);
  free(neighbors);
  free(rgMLDstr.vn);
  free(rgMLDstr.ve);
  free(rgMLDstr.vn_cov);
  free(rgMLDstr.ve_cov);

  fclose(vf);
}
