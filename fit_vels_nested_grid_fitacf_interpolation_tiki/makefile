CFLAGS= -D _GNU_SOURCE -D __USE_XOPEN  -I$(RSTPATH)/include -O3 -ipo -mcmodel medium -fPIC -Wall -D_GNU_SOURCE -D_LINUX -Wno-unused-but-set-variable -DMKL_ILP64  -qmkl-ilp64=parallel 

INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn -I$(IPATH)/analysis -I$(HOME)/src/lib -I$(HOME)/src/include -I$(HOME)/src/global_convection_maps/include

LIBS=-L/opt/intel/oneapi/compiler/latest/lib/ -L$(RSTPATH)/lib -lfit.1 -lrscan.1 -lradar.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 -laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lrpos.1 -lrfile.1 -lcfit.1 -ligrf.1 -lmlt_v2.1 -L$(HOME)/lib -lgrid -lgcc -laacgm.1


COMMON_FUNCTIONS=$(HOME)/src/global_convection_maps/common_routines


.c.o:
	icx $(CFLAGS) $(INCLUDE) -c -o $@ $<

OBJS=ml_df_fit.o $(COMMON_FUNCTIONS)/sub_sphazm.o $(COMMON_FUNCTIONS)/sub_sphcal.o $(COMMON_FUNCTIONS)/sub_gcp.o $(COMMON_FUNCTIONS)/pinv_mk.o $(COMMON_FUNCTIONS)/plms.o

ml_df_fit:	$(OBJS)
	icx -o ~/bin/ml_df_fit $(CFLAGS) $(INCLUDE) $(OBJS) -lm  -lz -qmkl-ilp64=parallel -qmkl=parallel -qopenmp -v $(LIBS) 

