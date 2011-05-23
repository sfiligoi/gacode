## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 180 2010-09-06 20:50:31Z kruger $
##
## ############################################################################

libdir = $(prefix)/lib
noinst_LTLIBRARIES = libgacodegyro_SRC.la

if HAVE_MUMPS
theCleanupMumps=cleanup_mumps.mumps.f90
theMumpsPrivate=gyro_mumps_private.mumps.f90
theSparseSolveMumps=sparse_solve_mumps.mumps.f90
else
theCleanupMumps=cleanup_mumps.f90
theMumpsPrivate=gyro_mumps_private.f90
theSparseSolveMumps=sparse_solve_mumps.f90
endif


if HAVE_MPIIO
theReadRestart=read_restart.mpiio.f90
theWriteRestart=write_restart.mpiio.f90
else
theReadRestart=read_restart.f90
theWriteRestart=write_restart.f90
endif

theDoNlfastMatch=do_nlfast_match.f90
theDoNlfastP=do_nlfast_p.f90

if HAVE_FFTW
theDoNlfftMatch=do_nlfft_match.fftw.f90
theDoNlfftP=do_nlfft_p.fftw.f90
theMakeNl=make_nl.fftw.f90
else
theDoNlfftMatch=do_nlfft_match.stub.f90
theDoNlfftP=do_nlfft_p.stub.f90
theMakeNl=make_nl.f90
endif

if HAVE_ESSL
theDoNlfftMatch=do_nlfft_match.essl.f90
theDoNlfftP=do_nlfft_p.essl.f90
theMakeNl=make_nl.essl.f90
else
theDoNlfftMatch=do_nlfft_match.stub.f90
theDoNlfftP=do_nlfft_p.stub.f90
theMakeNl=make_nl.f90
endif

gyro_includes= -I $(top_builddir)/shared/math \
	-I $(top_builddir)/shared/GEO \
	-I $(top_builddir)/gyro/BLEND \
	-I $(top_builddir)/gyro/SSUB \
	-I $(top_builddir)/gyro/TRANSPOSE \
	-I $(top_builddir)/shared/UMFPACK \
	-I $(top_builddir)/shared/EXPRO  \
	-I $(top_builddir)/shared/vshdf5

# To get to work with optimization reducing compiler wrappers
IOINC=@HDF5_FC_INC@ @HDF5_INC@
IOLIBS=@HDF5_ALIBS@ @NETCDF_ALIBS@
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.


libgacodegyro_SRC_la_SOURCES = \
	gyro_globals.f90 \
	gyro_profile_exp.f90 \
	gyro_pointers.f90 \
	gyro_nl_private.f90 \
	gyro_collision_private.f90 \
	$(theMumpsPrivate) \
	gyro_interface.f90 \
	gyro_alloc_profile_sim.f90 \
	gyro_alloc_big.f90 \
	gyro_alloc_orbit.f90 \
	gyro_alloc_velocity.f90 \
	gyro_alloc_distrib.f90 \
	gyro_alloc_nl.f90 \
	gyro_ballooning_mode.f90 \
  gyro_banana_private.f90 \
  gyro_banana_init.f90 \
  gyro_banana_integrate_tau.f90 \
  gyro_banana_operators.f90 \
  gyro_banana_getlambda.f90 \
  gyro_banana_s2lambda.f90 \
  gyro_banana_uniform_taugrid.f90 \
  gyro_bessel_operator.f90 \
  gyro_collect_routines.f90 \
	gyro_cleanup.f90 \
  gyro_omegas.f90 \
  gyro_rhs_total.f90 \
  gyro_tau_derivative.f90 \
  gyro_to_geo.f90 \
  gyro_write_initdata.f90 \
  gyro_write_input.f90 \
  gyro_write_master.f90 \
  gyro_write_precision.f90 \
  gyro_poisson_private.f90 \
  gyro_poissonaperp_private.f90 \
  gyro_profile_init.f90 \
  gyro_read_experimental_profiles.f90 \
  gyro_map_experimental_profiles.f90 \
	gyro_collision.f90 \
	gyro_collision_setup.f90 \
	gyro_collision_kernel.f90 \
	gyro_collision_ebelli.f90 \
	gyro_collision_setup_ebelli.f90 \
	gyro_collision_kernel_ebelli.f90 \
	gyro_conserve_number_ebelli.f90 \
	gyro_collision_grid_write.f90 \
	gyro_collision_krook.f90 \
	gyro_conserve_number.f90 \
	gyro_conserve_all.f90 \
	gyro_field_time_derivative.f90 \
	gyro_fieldeigen.f90 \
	gyro_memory_usage.f90 \
	gyro_nonlinear_flux.f90 \
	gyro_diffusivity.f90 \
	gyro_gbflux.f90 \
	gyro_initialize_arrays.f90 \
	gyro_radial_operators.f90 \
	gyro_moments_plot.f90 \
  gyro_moments_plot_wedge.f90 \
	gyro_select_methods.f90 \
	gyro_bounce_points.f90 \
	gyro_radial_simulation_box.f90 \
	gyro_read_input.f90 \
	gyro_read_input_extra.f90 \
	sparse_solve_umfpack.f90 \
	$(theSparseSolveMumps) \
	get_vel_sum_p.f90 \
	get_vel_sum_a.f90 \
	get_vel_sum_aperp.f90 \
	gyro_initialize_timestep.f90 \
	catch_error.f90 \
	catch_blowup.f90 \
	send_message_real.f90 \
	send_message.f90 \
	send_line.f90 \
	catch_halt_signal.f90 \
	parallel_dim.f90 \
	alloc_add.f90 \
	allocate_profile_exp.f90 \
	timestep_SSP_322.f90 \
	timestep_explicit.f90 \
  gyro_maxwell_private.f90 \
	make_maxwell_matrix.f90 \
	make_theta_grid.f90 \
	make_geometry_arrays.f90 \
	make_gyro.f90 \
	make_initial_h.f90 \
	make_implicit_advect.f90 \
	make_electron_current.f90 \
	make_electron_current_perp.f90 \
	make_ampere_matrix.f90 \
	make_poisson_matrix.f90 \
	make_ampere_blend.f90 \
	make_poisson_blend.f90 \
	make_ampereperp_blend.f90 \
	make_poissonaperp_matrix.f90 \
	$(theMakeNl) \
	get_nonlinear_advance.f90 \
	gyro_nonlinear_transfer.f90 \
	get_ampere_coll.f90 \
	get_gyro_h.f90 \
	get_gyro_h_aperp.f90 \
	get_maxwell_solution.f90 \
	get_poisson_solution.f90 \
	get_poisson_explicit.f90 \
	get_ampere_explicit.f90 \
	get_maxwell3_explicit.f90 \
	get_field_plot.f90 \
	get_field_r0_plot.f90 \
	get_nonlinear_flux_velocity.f90 \
	get_field_fluxave.f90 \
	get_field_explicit.f90 \
	get_kinetic_advance.f90 \
	get_delta_he.f90 \
	get_he.f90 \
	get_phi_squared.f90 \
	get_g_squared.f90 \
	$(theDoNlfastP) \
	$(theDoNlfastMatch) \
	$(theDoNlfftP) \
	$(theDoNlfftMatch) \
	$(theReadRestart) \
	set_exit_status.f90 \
	write_step.f90 \
	write_radial_operators.f90 \
	write_freq.f90 \
  write_banana_operators.f90 \
	write_h.f90 \
	write_efficiency.f90 \
	write_timing.f90 \
	$(theWriteRestart) \
	write_error.f90 \
	write_matrix_stat.f90 \
  gyro_adaptive_source.f90 \
	gyro_fieldeigen_kernel.f90 \
	gyro_collision_jc.f90 \
	gyro_build_sparse_ampere.f90 \
	gyro_build_sparse_poissonaperp.f90 \
  gyro_fieldeigen_df.f90 \
	gyro_fieldeigen_private.f90 \
	gyro_field_interpolation.f90 \
	gyro_build_sparse_maxwell.f90 \
	gyro_entropy.f90 \
	gyro_build_sparse_poisson.f90 \
  gyro_fulladvance.f90 \
  gyro_kxky_spectrum.f90 \
  gyro_lambda_grid.f90 \
  gyro_mpi_grid.f90 \
  gyro_set_pointers.f90 \
  gyro_timestep_error.f90 \
  gyro_set_phase_space.f90 \
  gyro_set_pointer_dim.f90 \
  gyro_set_blend_arrays.f90 \
	$(HDF5IO_FILES) \
	$(theCleanupMumps) \
  gyro_init.f90 \
	gyro_run.f90 \
  gyro_dump_input.f90 \
  gyro_dump_interface.f90 \
	gyro_do.f90  

if TX_FORTRAN_MODNAMEISCAP
nodist_include_HEADERS = GYRO_GLOBALS.$(TX_FORTRAN_MODEXT) GYRO_INTERFACE.$(TX_FORTRAN_MODEXT)
GYRO_GLOBALS.$(TX_FORTRAN_MODEXT):   gyro_globals.f90
GYRO_INTERFACE.$(TX_FORTRAN_MODEXT): gyro_interface.f90
else
nodist_include_HEADERS = gyro_globals.$(TX_FORTRAN_MODEXT) gyro_interface.$(TX_FORTRAN_MODEXT)
gyro_globals.$(TX_FORTRAN_MODEXT):   gyro_globals.f90
gyro_interface.$(TX_FORTRAN_MODEXT): gyro_interface.f90
endif


CLEANFILES = *.mod *.f90 *.lo *.lst *.la

gyro_fieldeigen.o: gyro_fieldeigen_private.o
gyro_fieldeigen.lo: gyro_fieldeigen_private.lo
