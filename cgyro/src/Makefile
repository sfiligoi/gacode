include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

ifdef FFTW_INC
   FFLAGS += -I${FFTW_INC}
endif

FC += ${FOMP} ${FACC}

EXEC = cgyro

OBJECTS = cgyro_globals.o \
	  cgyro_timer_lib.o \
	  cgyro_math.o \
	  cgyro_globals_math.o \
          cgyro_parallel_lib.o \
          cgyro_step.o \
          cgyro_io.o \
          cgyro_restart.o \
          cgyro_advect_wavenumber.o \
	  cgyro_init_rotation.o \
          cgyro_equilibrium.o \
          cgyro_check.o \
	  cgyro_check_memory.o \
          cgyro_error_estimate.o \
          cgyro_field.o \
          cgyro_field_coefficients.o \
          cgyro_flux.o \
          cgyro_freq.o \
          cgyro_global.o \
          cgyro_upwind.o \
          cgyro_init.o \
          cgyro_init_h.o \
          cgyro_init_arrays.o \
	  cgyro_init_collision_simple.o \
	  cgyro_init_collision_landau.o \
          cgyro_init_collision.o \
          cgyro_init_manager.o \
	  cgyro_init_kernel.o \
          cgyro_step_collision.o \
	  cgyro_kernel.o \
	  cgyro_cleanup.o \
          cgyro_final_kernel.o \
          cgyro_mpi_grid.o \
	  cgyro_make_profiles.o \
	  cgyro_nl_comm.o \
	  cgyro_nl_fftw.o \
          cgyro_read_input.o \
          cgyro_rhs.o \
          cgyro_run.o \
          cgyro_shear_hammett.o \
          cgyro_source.o \
	  cgyro_step_gk.o \
	  cgyro_step_gk_ck.o \
	  cgyro_step_gk_bs5.o \
	  cgyro_step_gk_v76.o \
          cgyro_write_initdata.o \
          cgyro_write_timedata.o \
          cgyro_write_hosts.o \
          cgyro_zftest_em.o 

.SUFFIXES : .o .f90 .F90

all: cgyro_lib.a $(EXEC).o $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o cgyro_lib.a $(EXTRA_LIBS) $(LMATH)

cgyro_lib.a: $(OBJECTS)
	$(ARCH) cgyro_lib.a $(OBJECTS)

.f90.o .F90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o cgyro_lib.a $(EXEC) *~ *.cub *.ptx *.gpu *.h *lst *.s
	cd ${GACODE_ROOT}/modules ; rm -f timer_lib.mod

-include ../install/make.ext.${GACODE_PLATFORM}

# cgyro_globals in cgyro_globals.F90
cgyro.o cgyro_advect_wavenumber.o cgyro_check.o cgyro_check_memory.o cgyro_cleanup.o cgyro_equilibrium.o cgyro_error_estimate.o cgyro_field.o cgyro_field_coefficients.o cgyro_final_kernel.o cgyro_flux.o cgyro_freq.o cgyro_global.o cgyro_globals_math.o cgyro_init.o cgyro_init_arrays.o cgyro_init_collision.o cgyro_init_collision_landau.o cgyro_init_collision_simple.o cgyro_init_h.o cgyro_init_kernel.o cgyro_init_manager.o cgyro_init_rotation.o cgyro_io.o cgyro_kernel.o cgyro_make_profiles.o cgyro_mpi_grid.o cgyro_nl_comm.o cgyro_nl_fftw.o cgyro_read_input.o cgyro_restart.o cgyro_rhs.o cgyro_run.o cgyro_shear_hammett.o cgyro_source.o cgyro_step_collision.o cgyro_step_gk.o cgyro_step_gk_bs5.o cgyro_step_gk_ck.o cgyro_step_gk_v76.o cgyro_upwind.o cgyro_write_hosts.o cgyro_write_initdata.o cgyro_write_timedata.o cgyro_zftest_em.o: cgyro_globals.o

# cgyro_globals_math in cgyro_globals_math.f90
cgyro_step_gk.o cgyro_step_gk_bs5.o cgyro_step_gk_ck.o cgyro_step_gk_v76.o: cgyro_globals_math.o

# cgyro_init_collision_landau in cgyro_init_collision_landau.F90
cgyro_init_collision.o: cgyro_init_collision_landau.o

# cgyro_io in cgyro_io.f90
cgyro.o cgyro_check.o cgyro_equilibrium.o cgyro_error_estimate.o cgyro_final_kernel.o cgyro_freq.o cgyro_init_collision.o cgyro_init_h.o cgyro_init_kernel.o cgyro_init_manager.o cgyro_init_rotation.o cgyro_kernel.o cgyro_make_profiles.o cgyro_mpi_grid.o cgyro_restart.o cgyro_step_gk_bs5.o cgyro_step_gk_ck.o cgyro_step_gk_v76.o cgyro_write_hosts.o cgyro_write_timedata.o cgyro_zftest_em.o: cgyro_io.o

# cgyro_math in cgyro_math.F90
cgyro_globals_math.o: cgyro_math.o

# cgyro_nl_comm in cgyro_nl_comm.F90
cgyro_nl_fftw.o cgyro_rhs.o: cgyro_nl_comm.o

# parallel_lib in cgyro_parallel_lib.F90
cgyro_cleanup.o cgyro_field.o cgyro_init_arrays.o cgyro_mpi_grid.o cgyro_nl_comm.o cgyro_nl_fftw.o cgyro_rhs.o cgyro_step_collision.o cgyro_upwind.o: cgyro_parallel_lib.o

# cgyro_restart in cgyro_restart.F90
cgyro_init_h.o cgyro_kernel.o: cgyro_restart.o

# cgyro_step in cgyro_step.f90
cgyro_init_kernel.o cgyro_kernel.o cgyro_step_gk_bs5.o cgyro_step_gk_ck.o cgyro_step_gk_v76.o cgyro_write_timedata.o: cgyro_step.o

# timer_lib in cgyro_timer_lib.f90
cgyro.o cgyro_advect_wavenumber.o cgyro_error_estimate.o cgyro_field.o cgyro_final_kernel.o cgyro_global.o cgyro_init.o cgyro_init_collision.o cgyro_init_collision_simple.o cgyro_init_kernel.o cgyro_init_manager.o cgyro_init_rotation.o cgyro_kernel.o cgyro_mpi_grid.o cgyro_nl_comm.o cgyro_nl_fftw.o cgyro_rhs.o cgyro_run.o cgyro_shear_hammett.o cgyro_source.o cgyro_step_collision.o cgyro_step_gk.o cgyro_step_gk_bs5.o cgyro_step_gk_ck.o cgyro_step_gk_v76.o cgyro_upwind.o cgyro_write_timedata.o: cgyro_timer_lib.o
