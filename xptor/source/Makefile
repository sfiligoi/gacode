include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB=../driver/xptor_lib

OBJECTS = neucg.o \
          neuden.o \
          neut.o \
	  palpha_dv.o \
	  prad_dv.o \
	  pohmic_dv.o \
          sigvth.o  

.SUFFIXES : .o .f 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst
