include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = tgyro_main

OBJECTS = tgyro_globals.o \
          tgyro_iteration_variables.o \
          tgyro_allocate_globals.o \
          tgyro_catch_error.o \
          tgyro_comm_setup.o \
          tgyro_comm_sync.o \
          tgyro_flux.o \
          tgyro_flux_vector.o \
          tgyro_flux_vector_dense.o \
          tgyro_gyro_map.o \
          tgyro_init_profiles.o \
          tgyro_iteration_driver.o \
          tgyro_iteration_standard.o \
          tgyro_iteration_pppl.o \
          tgyro_iteration_serial.o \
          tgyro_iteration_parallel.o \
          tgyro_multi_driver.o \
          tgyro_neo_map.o \
          tgyro_profile_functions.o \
          tgyro_qfm_map.o \
          tgyro_quasineutral.o \
          tgyro_quasigrad.o \
          tgyro_read_input.o \
          tgyro_readbc_int.o \
          tgyro_readbc_real.o \
	  tgyro_residual.o \
          tgyro_restart.o \
          tgyro_source.o \
          tgyro_stab_driver.o \
          tgyro_target_vector.o \
          tgyro_tglf_map.o \
          tgyro_volume_int.o \
          tgyro_write_input.o \
          tgyro_write_intermediate.o \
          tgyro_write_data.o \
          $(EXEC).o

.SUFFIXES : .o .f90 

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f $(EXEC)
	rm -f *.o ${GACODE_ROOT}/modules/tgyro*.mod


