include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = cgyro

OBJECTS = cgyro_globals.o \
          cgyro_io.o \
          cgyro_equilibrium.o \
          cgyro_check.o \
          cgyro_field.o \
          cgyro_flux.o \
          cgyro_freq.o \
          cgyro_hsym.o \
          cgyro_interface.o \
          cgyro_init.o \
          cgyro_init_h.o \
          cgyro_init_arrays.o \
          cgyro_init_collision.o \
          cgyro_step_collision.o \
	  cgyro_implicit_gk.o \
	  cgyro_implicit_partial_gk.o \
          cgyro_kernel.o \
          cgyro_mpi_grid.o \
	  cgyro_experimental_globals.o \
	  cgyro_experimental_alloc.o \
	  cgyro_experimental_profiles.o \
	  cgyro_experimental_map.o \
	  cgyro_make_profiles.o \
	  cgyro_nl.o \
          cgyro_read_input.o \
          cgyro_run.o \
          cgyro_shear.o \
	  cgyro_step_gk.o \
          cgyro_write_initdata.o \
          cgyro_read_restart.o \
          cgyro_write_restart.o \
          cgyro_write_timedata.o

.SUFFIXES : .o .f90

all: cgyro_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o cgyro_lib.a $(EXTRA_LIBS) $(LMATH)

cgyro_lib.a: $(OBJECTS)
	$(ARCH) cgyro_lib.a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -rf *.o cgyro_lib.a $(EXEC) *~
