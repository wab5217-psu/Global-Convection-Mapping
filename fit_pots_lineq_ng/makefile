CFLAGS= -D _GNU_SOURCE -D __USE_XOPEN  -I$(RSTPATH)/include -L /usr/local/lib -L/usr/intel -L$(RSTPATH)/lib -O3 -ipo -mcmodel medium -shared-intel -fPIC -Wall -D_GNU_SOURCE -D_LINUX


INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn -I$(IPATH)/analysis -I$(HOME)/src/lib -I$(HOME)/src/include -I$(HOME)/src/global_convection_maps/include

# LIBS=-lfit.1 -lrscan.1 -lradar.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 -laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lrpos.1  -lcnvmap.1 -loldcnvmap.1 -lshf.1 -lgrd.1 -loldgrd.1 -laacgm.1 -ldmap.1 -lrfile.1 -lopt.1 -lcfit.1 -ligrf.1

LIBS= -lcnvmap.1 -lshf.1 -loldcnvmap.1 -lgrd.1 -loldgrd.1 -lradar.1 \
      -ldmap.1 -lrfile.1 -lrtime.1 -lopt.1 -lrcnv.1 -laacgm.1 \
	  -laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lmlt_v2.1 -lrtime.1 -L$(HOME)/lib -lgrid



.c.o:
	icx $(CFLAGS) $(INCLUDE) -g -c -o $@ $<

OBJS=fit_ml_pot_lineq_ng.o pinv_mk.o

fit_ml_pot_lineq:	$(OBJS)
	icx -g -o ~/bin/fit_ml_pot_lineq $(CFLAGS) $(INCLUDE) $(OBJS) -lm  -lz -qmkl=parallel -qopenmp $(LIBS)


