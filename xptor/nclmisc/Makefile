include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
include ../install/make.ext.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} ${XPTORFLAGS}
else
   FFLAGS=${FOPT} ${XPTORFLAGS}
endif

LLIB = ../driver/xptor_lib

OBJECTS = 

.SUFFIXES : .o .f


lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst
