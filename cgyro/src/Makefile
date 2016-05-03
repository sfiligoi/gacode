include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

EXEC = cgyro

OBJECTS = cgyro_globals.o \
          cgyro_io.o \
          cgyro_equilibrium.o \
          cgyro_check.o \
	  cgyro_check_memory.o \
          cgyro_error_estimate.o \
          cgyro_field.o \
          cgyro_flux.o \
          cgyro_freq.o \
          cgyro_upwind.o \
          cgyro_init_h.o \
          cgyro_init_arrays.o \
	  cgyro_init_collision_simple.o \
          cgyro_init_collision.o \
          cgyro_init_manager.o \
	  cgyro_step_collision_simple.o \
          cgyro_step_collision.o \
	  cgyro_implicit_gk.o \
          cgyro_kernel.o \
          cgyro_mpi_grid.o \
	  cgyro_experimental_globals.o \
	  cgyro_experimental_alloc.o \
	  cgyro_experimental_profiles.o \
	  cgyro_experimental_map.o \
	  cgyro_make_profiles.o \
	  cgyro_nl_direct.o \
	  cgyro_nl_fftw.o \
          cgyro_read_input.o \
          cgyro_shear_interface.o \
          cgyro_shear.o \
          cgyro_shear_dft.o \
          cgyro_shear_pt.o \
	  cgyro_step_gk.o \
          cgyro_write_initdata.o \
          cgyro_read_restart.o \
          cgyro_write_restart.o \
          cgyro_write_timedata.o

.SUFFIXES : .o .f90 .F90

all: cgyro_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o cgyro_lib.a $(EXTRA_LIBS) $(LMATH)

cgyro_lib.a: $(OBJECTS)
	$(ARCH) cgyro_lib.a $(OBJECTS)

.f90.o .F90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -rf *.o cgyro_lib.a $(EXEC) *~ *.cub *.ptx

-include ../install/make.ext.${GACODE_PLATFORM}
