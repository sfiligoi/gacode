include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
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
	$(FC) $(FMATH) $(FFLAGS) -c $<

.F.o :
	$(FC) $(FMATH) $(FFLAGS)  -c $<

clean:
	rm -f *.o 
