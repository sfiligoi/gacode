include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

.SUFFIXES : .o .f90

all: vpro_lib.a

vpro_lib.a: gapy_vpro.o gapy_vpro_util.o
	$(ARCH) vpro_lib.a gapy_vpro.o gapy_vpro_util.o

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o *.a *.mod
