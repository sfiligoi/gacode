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

Z_LIBS	= @Z_LIBS@
Z_ALIBS	= @Z_ALIBS@

libdir = $(prefix)/lib/$(FC_LIBSUBDIR)
bin_PROGRAMS     = BigScience
noinst_LTLIBRARIES = libGYRO_SRC.la
lib_LTLIBRARIES = libGYRO.la

if HAVE_MUMPS
theCleanupMumps=cleanup_mumps.mumps.f90
theMumpsPrivate=mumps_private.mumps.f90
theSparseSolveMumps=sparse_solve_mumps.mumps.f90
else
theCleanupMumps=cleanup_mumps.f90
theMumpsPrivate=mumps_private.f90
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

if HAVE_LIBSCI
theDoNlfftMatch=do_nlfft_match.libsci.f90
theDoNlfftP=do_nlfft_p.libsci.f90
theMakeNl=make_nl.libsci.f90
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

if HAVE_HDF5
HDF5IO_FILES=hdf5_mod.f90 write_hdf5.F90 write_distributed_complex_h5.f90 write_distributed_real_h5.f90
else
HDF5IO_FILES=hdf5_mod.f90 write_hdf5.F90
endif
gyro_includes= -I $(top_builddir)/gyro/math \
	-I $(top_builddir)/gyro/GEO \
	-I $(top_builddir)/gyro/BLEND \
	-I $(top_builddir)/gyro/SSUB \
	-I $(top_builddir)/gyro/TRANSPOSE \
	-I $(top_builddir)/gyro/UMFPACK \
	-I $(top_builddir)/gyro/EXPRO \
	-I $(top_builddir)/gyro/EIGEN



# To get to work with optimization reducing compiler wrappers
IOINC=@HDF5_FC_INC@ @FCIOWRAPPERS_INC@
#IOINC=@HDF5_INC@ @FCIOWRAPPERS_INC@
IOLIBS=@FCIOWRAPPERS_ALIBS@ @HDF5_ALIBS@ @NETCDF_ALIBS@
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.


libGYRO_SRC_la_SOURCES = \
	gyro_globals.f90 \
	gyro_profile_exp.f90 \
	gyro_pointers.f90 \
	gyro_nl_private.f90 \
	gyro_collision_private.f90 \
	maxwell_private.f90 \
	poisson_private.f90 \
	poissonaperp_private.f90 \
	$(theMumpsPrivate) \
	umfpack_private.f90 \
	ORB_private.f90 \
	gyro_interface.f90 \
	gyro_init.f90 \
	gyro_run.f90 \
	gyro_do.f90 \
	gyro_alloc_profile_sim.f90 \
	gyro_alloc_big.f90 \
	gyro_alloc_orbit.f90 \
	gyro_alloc_velocity.f90 \
	gyro_alloc_distrib.f90 \
	gyro_alloc_nl.f90 \
	gyro_ballooning_mode.f90 \
	gyro_cleanup.f90 \
	gyro_neo_theory.f90 \
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
	gyro_ave.f90 \
	gyro_ave_rot.f90 \
	gyro_ave_aperp.f90 \
	gyro_bdoubleave.f90 \
	gyro_initialize_arrays.f90 \
	gyro_radial_operators.f90 \
	gyro_moments_plot.f90 \
	gyro_moments_plot_fine.f90 \
	gyro_select_methods.f90 \
	gyro_bounce_points.f90 \
	gyro_radial_simulation_box.f90 \
	gyro_read_input.f90 \
	gyro_dump_input.f90 \
	gyro_dump_interface.f90 \
	gyro_read_input_extra.f90 \
	sparse_solve_umfpack.f90 \
	$(theSparseSolveMumps) \
	get_vel_sum_p.f90 \
	get_vel_sum_a.f90 \
	get_vel_sum_aperp.f90 \
	initialize_timestep.f90 \
	set_profile_miller.f90 \
	filelog.f90 \
	catch_error.f90 \
	catch_blowup.f90 \
	send_message_real.f90 \
	send_message.f90 \
	send_line.f90 \
	catch_halt_signal.f90 \
	parallel_dim.f90 \
	collect_complex.f90 \
	collect_real.f90 \
	collect_integer.f90 \
	alloc_add.f90 \
	allocate_profile_exp.f90 \
	timestep_SSP_322.f90 \
	timestep_explicit.f90 \
	maxwell_p.f90 \
	maxwell_match.f90 \
	poisson_p.f90 \
	poisson_match.f90 \
	poissonaperp_p.f90 \
	poissonaperp_match.f90 \
	map_experimental_profiles.f90 \
	make_maxwell_matrix.f90 \
	make_MPI_grid.f90 \
	make_blend_arrays.f90 \
	make_theta_grid.f90 \
	make_omegas.f90 \
	make_profiles.f90 \
	make_profile_1.f90 \
	make_profile_3.f90 \
	make_experimental_profiles.f90 \
	make_experimental_neo_profiles.f90 \
	make_geometry_arrays.f90 \
	make_theta_operators.f90 \
	make_gyro.f90 \
	make_lambda_grid.f90 \
	make_phase_space.f90 \
	make_pointers.f90 \
	make_pointer_dimensions.f90 \
	make_dharmonic.f90 \
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
	make_phi_doppler.f90 \
	$(theMakeNl) \
	get_nonlinear_advance.f90 \
	get_nonlinear_transfer.f90 \
	get_ampere_coll.f90 \
	get_gyro_h.f90 \
	get_gyro_h_aperp.f90 \
	get_maxwell_solution.f90 \
	get_poisson_solution.f90 \
	get_poisson_explicit.f90 \
	get_ampere_explicit.f90 \
	get_maxwell3_explicit.f90 \
	get_RHS.f90 \
	get_field_spectrum.f90 \
	get_field_plot.f90 \
	get_field_r0_plot.f90 \
	get_nonlinear_flux.f90 \
	get_nonlinear_flux_velocity.f90 \
	get_field_fluxave.f90 \
	get_diffusivity.f90 \
	get_field_explicit.f90 \
	get_kinetic_advance.f90 \
	get_delta_he.f90 \
	get_he.f90 \
	get_error.f90 \
	get_phi_squared.f90 \
	get_g_squared.f90 \
	ORB_do.f90 \
	ORB_init.f90 \
	ORB_lambda.f90 \
	ORB_s2lambda.f90 \
	ORB_s.f90 \
	do_adaptive_source.f90 \
	do_dtau.f90 \
	do_fulladvance.f90 \
	$(theDoNlfastP) \
	$(theDoNlfastMatch) \
	$(theDoNlfftP) \
	$(theDoNlfftMatch) \
	readbc_int.f90 \
	readbc_real.f90 \
	read_experimental_profiles.f90 \
	$(theReadRestart) \
	set_exit_status.f90 \
	write_big.f90 \
	write_step.f90 \
	write_input.f90 \
	write_theta_operators.f90 \
	write_radial_operators.f90 \
	write_freq.f90 \
	write_ballooning_mode.f90 \
	write_sp_diffusivity.f90 \
	write_diffusivity_i_ch.f90 \
	write_sp_diffusivity_i.f90 \
	write_h.f90 \
	write_efficiency.f90 \
	write_profile_sim.f90 \
	write_profile_vugyro.f90 \
	write_timing.f90 \
	$(theWriteRestart) \
	write_error.f90 \
	write_prec.f90 \
	write_geometry_arrays.f90 \
	write_matrix_stat.f90 \
	write_distributed_complex.f90 \
	write_distributed_real.f90 \
	write_local_real.f90 \
	write_units.f90 \
	gyro_fieldeigen_kernel.f90 \
	gyro_collision_jc.f90 \
	gyro_build_sparse_ampere.f90 \
	gyro_build_sparse_poissonaperp.f90 \
	gyro_fieldeigen_private.f90 \
	gyro_field_interpolation.f90 \
	gyro_build_sparse_maxwell.f90 \
	gyro_entropy.f90 \
	write_geometry_arrays.f90.ps \
	gyro_build_sparse_poisson.f90 \
	$(HDF5IO_FILES) \
	$(theCleanupMumps) \
	../EIGEN/EIGEN_do.stub.f90

libGYRO_la_SOURCES =
libGYRO_la_LIBADD = \
	libGYRO_SRC.la \
	$(top_builddir)/gyro/BLEND/libBLEND.la \
	$(top_builddir)/gyro/math/libgyromath.la \
	$(top_builddir)/gyro/TRANSPOSE/libTRANSPOSE.la \
	$(top_builddir)/gyro/SSUB/libSSUB.la \
	$(top_builddir)/gyro/GEO/libGEO.la \
	$(top_builddir)/gyro/EXPRO/libEXPRO.la \
	$(top_builddir)/gyro/UMFPACK/libUMFPACK.la


BigScience_SOURCES = BigScience.f90
BigScience_LDADD  = $(top_builddir)/gyro/src/libGYRO.la \
	@FFTW_LIBS@  @LAPACK_LIBS@  @BLAS_LIBS@ $(IOLIBS) $(Z_LIBS) \
	$(top_builddir)/gyro/EXPRO/libEXPRO.la

BigScience_INCLUDES = \
	-I $(top_builddir)/gyro/src \
	-I $(top_builddir)/gyro/BLEND \
	-I $(top_builddir)/gyro/math \
	-I $(top_builddir)/gyro/TRANSPOSE \
	-I $(top_builddir)/gyro/SSUB \
	-I $(top_builddir)/gyro/GEO \
	-I $(top_builddir)/gyro/UMFPACK \
	-I $(top_builddir)/gyro/EXPRO

if TX_FORTRAN_MODNAMEISCAP
include_HEADERS = GYRO_GLOBALS.$(TX_FORTRAN_MODEXT) GYRO_INTERFACE.$(TX_FORTRAN_MODEXT)
GYRO_GLOBALS.$(TX_FORTRAN_MODEXT):   gyro_globals.f90
GYRO_INTERFACE.$(TX_FORTRAN_MODEXT): gyro_interface.f90
else
include_HEADERS = gyro_globals.$(TX_FORTRAN_MODEXT) gyro_interface.$(TX_FORTRAN_MODEXT)
gyro_globals.$(TX_FORTRAN_MODEXT):   gyro_globals.f90
gyro_interface.$(TX_FORTRAN_MODEXT): gyro_interface.f90
endif


CLEANFILES = *.mod *.f90 *.lo *.lst *.la

write_hdf5.o: hdf5_mod.o
write_hdf5.lo: hdf5_mod.lo
gyro_fieldeigen.o: gyro_fieldeigen_private.o
gyro_fieldeigen.lo: gyro_fieldeigen_private.lo
