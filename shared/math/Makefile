include ../install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif
      
LLIB=math_lib

OBJECTS = math_constants.o \
          bound_deriv.o \
          bound_extrap.o \
          cub_spline.o \
          cub_spline_deriv.o \
          cub_spline_deriv2.o \
          curtis_rad.o \
          energy_integral.o \
          erf.o \
          gacode_system.o \
          gamma.o \
          gauss_legendre.o \
          gauss_integ.o \
          ifs_pppl.o \
          invert_p32.o \
          i0.o \
          j0y0.o \
          linint.o \
          logint.o \
          parallel_dim.o \
          pascal.o \
          polydiff.o \
          poly2diff.o \
          proc_time.o \
          p32.o \
          parallel_lib.o \
          pseudospec_lib.o \
          remap_grid.o \
          rjbesl.o \
          sigv.o \
          sivukhin.o \
          timer_lib.o

.SUFFIXES : .o .f90 .f

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<
 
clean:
	rm -f *.o $(LLIB).a *.mod
	rm -f ${GACODE_ROOT}/modules/math_constants.mod
	rm -f ${GACODE_ROOT}/modules/parallel_lib.mod
	rm -f ${GACODE_ROOT}/modules/timer_lib.mod

-include install/make.ext.${GACODE_PLATFORM}
