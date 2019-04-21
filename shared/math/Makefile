include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT} -fPIC
endif

LLIB=math_lib

OBJECTS = math_constants.o \
          energy_integral.o \
          gauss_legendre.o \
          i0.o \
          parallel_dim.o

.SUFFIXES : .o .f90 .F90 .f

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a *.mod *.lst
	rm -f ${GACODE_ROOT}/modules/math_constants.mod
	rm -f ${GACODE_ROOT}/modules/timer_lib.mod

-include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
