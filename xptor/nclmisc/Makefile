include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB = ../driver/xptor_lib

OBJECTS = fluxav_efit.o \
          fluxav_vmec.o \
          geom_circ.o \
          geom_efit.o \
          geom_vmec.o \
          get_pro_xptor.o \
	  read_eqdsk.o \
          read_efit_eqdsk.o \
          read_vmec_wout.o \
	  store_efit.o \
          track_mod.o

.SUFFIXES : .o .f


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst
