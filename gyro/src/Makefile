include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif
      
LLIB = gyro_lib

OBJECTS = gyro_globals.o \
          gyro_profile_exp.o \
          gyro_pointers.o \
          gyro_nl_private.o \
          gyro_collision_private.o \
          gyro_fieldeigen_private.o \
          gyro_maxwell_private.o \
          gyro_poisson_private.o \
	  gyro_poissonaperp_private.o \
          gyro_mumps_private.o \
          gyro_banana_private.o \
          gyro.o \
          gyro_adaptive_source.o \
          gyro_interface.o \
          gyro_init.o \
          gyro_run.o \
          gyro_do.o \
          gyro_alloc_profile_sim.o \
          gyro_alloc_big.o \
          gyro_alloc_orbit.o \
          gyro_alloc_velocity.o \
          gyro_alloc_distrib.o \
          gyro_alloc_nl.o \
          gyro_ballooning_mode.o \
          gyro_banana_init.o \
          gyro_banana_integrate_tau.o \
          gyro_banana_getlambda.o \
          gyro_banana_operators.o \
          gyro_banana_uniform_taugrid.o \
          gyro_banana_s2lambda.o \
          gyro_bessel_operator.o \
          gyro_bessel_stencils.o \
          gyro_bounce_points.o \
          gyro_build_sparse_ampere.o \
          gyro_build_sparse_poisson.o \
          gyro_build_sparse_poissonaperp.o \
          gyro_build_sparse_maxwell.o \
          gyro_cleanup.o \
          gyro_collect_routines.o \
          gyro_collision.o \
          gyro_collision_setup.o \
          gyro_collision_kernel.o \
	  gyro_collision_ebelli.o \
	  gyro_collision_setup_ebelli.o \
	  gyro_collision_kernel_ebelli.o \
	  gyro_conserve_number_ebelli.o \
          gyro_collision_jc.o \
          gyro_collision_grid_write.o \
          gyro_collision_krook.o \
          gyro_conserve_number.o \
          gyro_conserve_all.o \
          gyro_diffusivity.o \
          gyro_entropy.o \
          gyro_field_time_derivative.o \
          gyro_field_interpolation.o \
	  gyro_fieldeigen.o \
	  gyro_fieldeigen_kernel.o \
	  gyro_fieldeigen_df.o \
          gyro_fulladvance.o \
          gyro_get_he_implicit.o \
          gyro_gbflux.o \
          gyro_initial_condition.o \
          gyro_initialize_arrays.o \
          gyro_initialize_timestep.o \
          gyro_input_check.o \
          gyro_kxky_spectrum.o \
          gyro_lambda_grid.o \
          gyro_map_experimental_profiles.o \
          gyro_message_routines.o \
          gyro_moments_plot.o \
          gyro_moments_plot_wedge.o \
          gyro_mpi_grid.o \
          gyro_memory_usage.o \
          gyro_nl_setup.o \
          gyro_nl_direct.o \
          gyro_nl_fft.o \
          gyro_nonlinear_flux.o \
          gyro_nonlinear_transfer.o \
          gyro_omegas.o \
          gyro_profile_init.o \
          gyro_radial_operators.o \
          gyro_radial_simulation_box.o \
          gyro_read_input.o \
          gyro_read_input_extra.o \
          gyro_read_experimental_profiles.o \
          gyro_read_restart.o \
          gyro_set_exit_status.o \
          gyro_write_restart.o \
          gyro_rhs_nonlinear.o \
          gyro_rhs_total.o \
          gyro_select_methods.o \
          gyro_set_blend_arrays.o \
          gyro_set_phase_space.o \
          gyro_set_pointers.o \
          gyro_set_pointer_dim.o \
          gyro_tau_derivative.o \
          gyro_theta_grid.o \
          gyro_timestep_error.o \
          gyro_timestep_explicit.o \
          gyro_timestep_implicit.o \
          gyro_to_geo.o \
          gyro_write_efficiency.o \
          gyro_write_initdata.o \
          gyro_write_input.o \
          gyro_write_precision.o \
          gyro_write_radial_op.o \
          gyro_write_timedata.o \
          sparse_solve_umfpack.o \
          sparse_solve_mumps.o \
          get_vel_sum_p.o \
          get_vel_sum_a.o \
	  get_vel_sum_aperp.o \
          catch_error.o \
          catch_blowup.o \
          catch_halt_signal.o \
          parallel_dim.o \
          alloc_add.o \
          allocate_profile_exp.o \
          make_maxwell_matrix.o \
          make_geometry_arrays.o \
          make_implicit_advect.o \
          make_electron_current.o \
	  make_electron_current_perp.o \
          make_ampere_matrix.o \
          make_poisson_matrix.o \
          make_ampere_blend.o \
          make_poisson_blend.o \
	  make_ampereperp_blend.o \
	  make_poissonaperp_matrix.o \
          get_ampere_coll.o \
          get_gyro_h.o \
	  get_gyro_h_aperp.o \
          get_maxwell_solution.o \
          get_poisson_solution.o \
          get_poisson_explicit.o \
          get_ampere_explicit.o \
	  get_maxwell3_explicit.o \
          get_field_plot.o \
          get_field_r0_plot.o \
          get_nonlinear_flux_velocity.o \
          get_field_fluxave.o \
          get_field_explicit.o \
          get_kinetic_advance.o \
          get_delta_he.o \
          get_phi_squared.o \
          get_g_squared.o \
          write_step.o \
          write_freq.o \
          write_h.o \
          write_timing.o \
          write_error.o \
          write_matrix_stat.o \
          gyro_dump_input.o \
          gyro_dump_interface.o

.SUFFIXES : .o .f90 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/gyro*.mod

include ../install/make.ext.${GACODE_PLATFORM}
