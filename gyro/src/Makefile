include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

ifdef FFTW_INC
   FFLAGS += -I${FFTW_INC}
endif

FC += ${FOMP}

LLIB = gyro_lib

OBJECTS += gyro_globals.o \
          gyro_profile_exp.o \
          gyro_pointers.o \
          gyro_nl_private.o \
          gyro_collision_private.o \
          gyro_fieldeigen_private.o \
          gyro_maxwell_private.o \
          gyro_poisson_private.o \
          gyro_poissonaperp_private.o \
          gyro_banana_private.o \
          gyro_module_ompdata.o \
          gyro.o \
          gyro_adaptive_source.o \
          gyro_interface.o \
          gyro_init.o \
          gyro_init_ompdata.o \
          gyro_run.o \
          gyro_do.o \
          gyro_alloc_profile_exp.o \
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
          gyro_blend_ampere.o \
          gyro_blend_ampereperp.o \
          gyro_blend_poisson.o \
          gyro_bounce_points.o \
          gyro_build_sparse_ampere.o \
          gyro_build_sparse_poisson.o \
          gyro_build_sparse_poissonaperp.o \
          gyro_build_sparse_maxwell.o \
          gyro_catch_error_routines.o \
          gyro_cleanup.o \
          gyro_collect_routines.o \
          gyro_collision_setup.o \
          gyro_collision_kernel.o \
          gyro_collision_main.o \
          gyro_collision_grid_write.o \
          gyro_conserve_number.o \
          gyro_conserve_all.o \
          gyro_entropy.o \
          gyro_field_fluxave.o \
          gyro_field_time_derivative.o \
          gyro_field_interpolation.o \
          gyro_field_plot.o \
          gyro_field_solve_explicit.o \
          gyro_field_solve_implicit.o \
          gyro_fieldeigen.o \
          gyro_fieldeigen_kernel.o \
          gyro_fieldeigen_df.o \
          gyro_fulladvance.o \
          gyro_geometry_arrays.o \
          gyro_get_he_implicit.o \
          gyro_get_delta_he.o \
          gyro_gbflux.o \
          gyro_initial_condition.o \
          gyro_initialize_arrays.o \
          gyro_initialize_timestep.o \
          gyro_kxky_spectrum.o \
          gyro_kinetic_advance.o \
          gyro_lambda_grid.o \
          gyro_make_maxwell_matrix.o \
          gyro_make_poisson_matrix.o \
          gyro_make_ampere_matrix.o \
          gyro_make_poissonaperp_matrix.o \
          gyro_make_implicit_advect.o \
          gyro_make_jpar_elec.o \
          gyro_make_jperp_elec.o \
          gyro_map_experimental_profiles.o \
          gyro_message_routines.o \
          gyro_moments_plot.o \
          gyro_mpi_grid.o \
          gyro_memory_usage.o \
          gyro_miscellaneous_math.o \
          gyro_nl_setup.o \
          gyro_nl_direct.o \
          gyro_nl_fft.o \
          gyro_nonlinear_flux.o \
          gyro_omegas.o \
          gyro_operators_on_h.o \
          gyro_phi_kp_squared.o \
          gyro_profile_init.o \
          gyro_radial_operators.o \
          gyro_radial_simulation_box.o \
          gyro_read_input.o \
          gyro_read_input_extra.o \
          gyro_read_experimental_profiles.o \
          gyro_read_restart.o \
          gyro_set_exit_status.o \
          gyro_sparse_solve_umfpack.o \
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
          gyro_timer.o \
          gyro_timestep_error.o \
          gyro_timestep_explicit.o \
          gyro_timestep_implicit.o \
          gyro_to_geo.o \
          gyro_write_efficiency.o \
          gyro_write_error.o \
          gyro_write_freq.o \
          gyro_write_h.o \
          gyro_write_initdata.o \
          gyro_write_input.o \
          gyro_write_precision.o \
          gyro_write_radial_op.o \
          gyro_write_step.o \
          gyro_write_timedata.o \
          gyro_write_timers.o \
          gyro_velocity_sum.o \
          i0.o


.SUFFIXES : .o .f90

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o .F90.o:
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o *.gpu *.h $(LLIB).a gyro
	rm -f ${GACODE_ROOT}/modules/gyro*.mod ${GACODE_ROOT}/modules/ompdata.mod

-include ../install/make.ext.${GACODE_PLATFORM}

# all module induced dependencies in directory gyro/src

# module gyro_banana_private in gyro_banana_private.f90
gyro_banana_getlambda.o gyro_banana_init.o gyro_banana_integrate_tau.o gyro_banana_s2lambda.o gyro_banana_uniform_taugrid.o: gyro_banana_private.o

# module gyro_collision_private in gyro_collision_private.f90
gyro_alloc_orbit.o gyro_build_sparse_ampere.o gyro_collision_grid_write.o gyro_collision_kernel.o gyro_collision_setup.o gyro_fieldeigen.o gyro_make_ampere_matrix.o gyro_memory_usage.o: gyro_collision_private.o

# module gyro_fieldeigen_private in gyro_fieldeigen_private.f90
gyro_fieldeigen.o gyro_fieldeigen_df.o gyro_fieldeigen_kernel.o gyro_run.o: gyro_fieldeigen_private.o

# module gyro_globals in gyro_globals.f90
gyro.o gyro_adaptive_source.o gyro_alloc_big.o gyro_alloc_distrib.o gyro_alloc_nl.o gyro_alloc_orbit.o gyro_alloc_profile_exp.o gyro_alloc_profile_sim.o gyro_alloc_velocity.o gyro_ballooning_mode.o gyro_banana_operators.o gyro_bessel_operator.o gyro_bessel_stencils.o gyro_blend_ampere.o gyro_blend_ampereperp.o gyro_blend_poisson.o gyro_build_sparse_ampere.o gyro_build_sparse_maxwell.o gyro_build_sparse_poisson.o gyro_build_sparse_poissonaperp.o gyro_catch_error_routines.o gyro_cleanup.o gyro_collect_routines.o gyro_collision_grid_write.o gyro_collision_kernel.o gyro_collision_main.o gyro_collision_setup.o gyro_conserve_all.o gyro_conserve_number.o gyro_do.o gyro_entropy.o gyro_field_fluxave.o gyro_field_interpolation.o gyro_field_plot.o gyro_field_solve_explicit.o gyro_field_solve_implicit.o gyro_field_time_derivative.o gyro_fieldeigen.o gyro_fieldeigen_df.o gyro_fieldeigen_kernel.o gyro_fulladvance.o gyro_gbflux.o gyro_geometry_arrays.o gyro_get_delta_he.o gyro_get_he_implicit.o gyro_init.o gyro_init_ompdata.o gyro_initial_condition.o gyro_initialize_arrays.o gyro_initialize_timestep.o gyro_interface.o gyro_kinetic_advance.o gyro_kxky_spectrum.o gyro_lambda_grid.o gyro_make_ampere_matrix.o gyro_make_implicit_advect.o gyro_make_jpar_elec.o gyro_make_jperp_elec.o gyro_make_maxwell_matrix.o gyro_make_poisson_matrix.o gyro_make_poissonaperp_matrix.o gyro_map_experimental_profiles.o gyro_memory_usage.o gyro_message_routines.o gyro_moments_plot.o gyro_mpi_grid.o gyro_nl_direct.o gyro_nl_fft.o gyro_nl_setup.o gyro_nonlinear_flux.o gyro_omegas.o gyro_operators_on_h.o gyro_phi_kp_squared.o gyro_profile_init.o gyro_radial_operators.o gyro_radial_simulation_box.o gyro_read_experimental_profiles.o gyro_read_input.o gyro_read_input_extra.o gyro_read_restart.o gyro_rhs_nonlinear.o gyro_rhs_total.o gyro_run.o gyro_select_methods.o gyro_set_blend_arrays.o gyro_set_exit_status.o gyro_set_phase_space.o gyro_set_pointer_dim.o gyro_set_pointers.o gyro_sparse_solve_umfpack.o gyro_tau_derivative.o gyro_theta_grid.o gyro_timer.o gyro_timestep_error.o gyro_timestep_explicit.o gyro_timestep_implicit.o gyro_to_geo.o gyro_velocity_sum.o gyro_write_efficiency.o gyro_write_error.o gyro_write_freq.o gyro_write_h.o gyro_write_initdata.o gyro_write_input.o gyro_write_precision.o gyro_write_radial_op.o gyro_write_restart.o gyro_write_step.o gyro_write_timedata.o gyro_write_timers.o: gyro_globals.o

# module gyro_interface in gyro_interface.f90
gyro_cleanup.o gyro_do.o gyro_init.o gyro_run.o: gyro_interface.o

# module gyro_maxwell_private in gyro_maxwell_private.f90
gyro_build_sparse_maxwell.o gyro_field_solve_implicit.o gyro_make_maxwell_matrix.o: gyro_maxwell_private.o

# module ompdata in gyro_module_ompdata.f90
gyro.o gyro_collision_main.o gyro_conserve_number.o gyro_field_time_derivative.o gyro_get_delta_he.o gyro_get_he_implicit.o gyro_init.o gyro_init_ompdata.o gyro_kxky_spectrum.o gyro_moments_plot.o gyro_nl_direct.o gyro_nl_fft.o gyro_nl_setup.o gyro_nonlinear_flux.o gyro_phi_kp_squared.o gyro_rhs_total.o gyro_select_methods.o gyro_tau_derivative.o gyro_timestep_error.o gyro_velocity_sum.o gyro_write_input.o: gyro_module_ompdata.o

# module gyro_nl_private in gyro_nl_private.f90
gyro_alloc_nl.o gyro_nl_direct.o gyro_nl_fft.o gyro_nl_setup.o: gyro_nl_private.o

# module gyro_pointers in gyro_pointers.f90
gyro_adaptive_source.o gyro_alloc_big.o gyro_alloc_distrib.o gyro_bessel_stencils.o gyro_blend_ampere.o gyro_blend_ampereperp.o gyro_blend_poisson.o gyro_collision_kernel.o gyro_collision_main.o gyro_collision_setup.o gyro_conserve_all.o gyro_conserve_number.o gyro_do.o gyro_entropy.o gyro_field_fluxave.o gyro_field_interpolation.o gyro_field_plot.o gyro_fieldeigen.o gyro_fieldeigen_df.o gyro_fieldeigen_kernel.o gyro_fulladvance.o gyro_get_delta_he.o gyro_get_he_implicit.o gyro_initial_condition.o gyro_make_implicit_advect.o gyro_make_jpar_elec.o gyro_make_jperp_elec.o gyro_memory_usage.o gyro_moments_plot.o gyro_mpi_grid.o gyro_nl_direct.o gyro_nl_fft.o gyro_nonlinear_flux.o gyro_omegas.o gyro_operators_on_h.o gyro_phi_kp_squared.o gyro_read_restart.o gyro_rhs_nonlinear.o gyro_rhs_total.o gyro_set_blend_arrays.o gyro_set_pointer_dim.o gyro_set_pointers.o gyro_tau_derivative.o gyro_timestep_error.o gyro_timestep_explicit.o gyro_timestep_implicit.o gyro_velocity_sum.o gyro_write_h.o gyro_write_input.o gyro_write_precision.o gyro_write_restart.o: gyro_pointers.o

# module gyro_poisson_private in gyro_poisson_private.f90
gyro_build_sparse_poisson.o gyro_make_poisson_matrix.o: gyro_poisson_private.o

# module gyro_poissonaperp_private in gyro_poissonaperp_private.f90
gyro_build_sparse_poissonaperp.o gyro_make_poissonaperp_matrix.o: gyro_poissonaperp_private.o

# module gyro_profile_exp in gyro_profile_exp.f90
gyro_alloc_profile_exp.o gyro_map_experimental_profiles.o gyro_profile_init.o gyro_read_experimental_profiles.o gyro_write_initdata.o gyro_write_input.o: gyro_profile_exp.o
