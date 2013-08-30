include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB = ../driver/xptor_lib

OBJECTS = glf2d_dv.o \
          glf2d.o    \
	  cgg_solver.o \
          neoclass.o \
          neoclass_one.o \
	  xptor_neo_map.o \
          ip_chi1.o \
          ip_chi2.o 


.SUFFIXES : .o .f .f90 .F


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(F77) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FFLAGS) -c $<

.F.o :
	$(F77) $(FFLAGS)  -c $<

clean:
	rm -f *.o *.lst
