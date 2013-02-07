include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB=../driver/xptor_lib

OBJECTS = advancefields.o \
          norm_dv.o \
          ptor_dv.o \
          ptorinit_dv.o \
          trcoef_dv.o \
          vtrans_dv.o \
          tglf_startup.o \
	  neo_flows.o

.SUFFIXES : .o .f 


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(F77) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst
