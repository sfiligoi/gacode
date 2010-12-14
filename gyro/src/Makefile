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
          gyro_umfpack_private.o \
          ORB_private.o \
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
          gyro_ave.o \
          gyro_ave_rot.o \
	  gyro_ave_aperp.o \
          gyro_ballooning_mode.o \
          gyro_build_sparse_ampere.o \
          gyro_build_sparse_poisson.o \
          gyro_build_sparse_poissonaperp.o \
          gyro_build_sparse_maxwell.o \
          gyro_cleanup.o \
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
          gyro_entropy.o \
          gyro_field_time_derivative.o \
          gyro_field_interpolation.o \
	  gyro_fieldeigen.o \
	  gyro_fieldeigen_kernel.o \
	  gyro_fieldeigen_df.o \
          gyro_fulladvance.o \
          gyro_memory_usage.o \
          gyro_nonlinear_flux.o \
          gyro_omegas.o \
          gyro_to_geo.o \
          gyro_profile_init.o \
          gyro_diffusivity.o \
          gyro_gbflux.o \
          gyro_bdoubleave.o \
          gyro_initialize_arrays.o \
          gyro_radial_operators.o \
          gyro_map_experimental_profiles.o \
          gyro_moments_plot.o \
          gyro_select_methods.o \
          gyro_bounce_points.o \
          gyro_radial_simulation_box.o \
          gyro_read_input.o \
          gyro_read_input_extra.o \
          gyro_read_experimental_profiles.o \
          gyro_rhs_total.o \
          gyro_write_input.o \
          gyro_write_precision.o \
          gyro_write_units.o \
          gyro_write_master.o \
          sparse_solve_umfpack.o \
          sparse_solve_mumps.o \
          get_vel_sum_p.o \
          get_vel_sum_a.o \
	  get_vel_sum_aperp.o \
          initialize_timestep.o \
          catch_error.o \
          catch_blowup.o \
          send_message_real.o \
          send_message.o \
          send_line.o \
          catch_halt_signal.o \
          parallel_dim.o \
          collect_complex.o \
          collect_real.o \
          collect_integer.o \
          alloc_add.o \
          allocate_profile_exp.o \
          timestep_SSP_322.o \
          timestep_explicit.o \
          make_maxwell_matrix.o \
          make_MPI_grid.o \
          make_blend_arrays.o \
          make_theta_grid.o \
          make_geometry_arrays.o \
          make_theta_operators.o \
          make_gyro.o \
          make_lambda_grid.o \
          make_phase_space.o \
          make_pointers.o \
          make_pointer_dimensions.o \
          make_dharmonic.o \
          make_initial_h.o \
          make_implicit_advect.o \
          make_electron_current.o \
	  make_electron_current_perp.o \
          make_ampere_matrix.o \
          make_poisson_matrix.o \
          make_ampere_blend.o \
          make_poisson_blend.o \
	  make_ampereperp_blend.o \
	  make_poissonaperp_matrix.o \
          make_nl.o \
          get_nonlinear_advance.o \
          get_nonlinear_transfer.o \
          get_ampere_coll.o \
          get_gyro_h.o \
	  get_gyro_h_aperp.o \
          get_maxwell_solution.o \
          get_poisson_solution.o \
          get_poisson_explicit.o \
          get_ampere_explicit.o \
	  get_maxwell3_explicit.o \
          get_field_spectrum.o \
          get_field_plot.o \
          get_field_r0_plot.o \
          get_nonlinear_flux_velocity.o \
          get_field_fluxave.o \
          get_field_explicit.o \
          get_kinetic_advance.o \
          get_delta_he.o \
          get_he.o \
          get_error.o \
          get_phi_squared.o \
          get_g_squared.o \
          ORB_do.o \
          ORB_init.o \
          ORB_lambda.o \
          ORB_s2lambda.o \
          ORB_s.o \
          do_dtau.o \
          do_nlfast_p.o \
          do_nlfast_match.o \
          do_nlfft_p.o \
          do_nlfft_match.o \
          read_restart.o \
          set_exit_status.o \
          write_step.o \
          write_theta_operators.o \
          write_radial_operators.o \
          write_freq.o \
          write_h.o \
          write_efficiency.o \
          write_profile_sim.o \
          write_profile_vugyro.o \
          write_timing.o \
          write_restart.o \
          write_error.o \
          write_geometry_arrays.o \
          write_matrix_stat.o \
          write_distributed_complex.o \
          write_distributed_real.o \
          write_local_real.o \
          cleanup_mumps.o 

.SUFFIXES : .o .f90 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a ${GACODE_ROOT}/modules/gyro*.mod

include ../install/make.ext.${GACODE_PLATFORM}
