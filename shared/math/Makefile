include ../install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif
      
LLIB=math_lib

OBJECTS = math_constants.o \
          j0y0.o \
          i0.o \
          erf.o \
          proc_time.o \
	  parallel_dim.o \
          p32.o \
          invert_p32.o \
          gauss_legendre.o \
          gauss_integ.o \
          polydiff.o \
          poly2diff.o \
          pascal.o \
          bound_extrap.o \
          bound_deriv.o \
          remap_grid.o \
          rjbesl.o \
          ribesl.o \
          cub_spline.o \
          cub_spline_deriv.o \
          cub_spline_deriv2.o \
          gamma.o \
          energy_integral.o \
          curtis_rad.o \
          linint.o \
          logint.o

.SUFFIXES : .o .f90 .f

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<
 
clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/math_constants.mod
