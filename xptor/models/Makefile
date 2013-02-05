include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB = ../driver/xptor_lib

OBJECTS = callglf2d.o \
          etaw17a.o \
          forcebal.o \
          glf2d_dv.o \
	  glf2d.o \
	  cgg_solver.o \
	  glf2d_max.o \
          ifsppl_dv.o \
          ip_chi1.o \
          ip_chi2.o \
          mmm95.o \
          model.o \
          nclass_mod.o \
          neoclass.o \
          neoclass_one.o \
	  xptor_neo_map.o \
          weiland14.o

.SUFFIXES : .o .f .F


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(FC) $(FFLAGS) -c $<

.F.o :
	$(FC) $(FFLAGS)  -c $<

clean:
	rm -f *.o *.lst
