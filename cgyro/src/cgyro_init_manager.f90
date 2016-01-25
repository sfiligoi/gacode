!-----------------------------------------------------------------
! cgyro_init_manager.f90
!
! PURPOSE:
!  Manage initialization of arrays and other setup issues.
!  In particular,
!
!  1. collisionless streaming arrays (timed setup)
!  2. collision arrays (timed)
!  3. write initial data
!  4. Construct initial distributions
!-----------------------------------------------------------------

subroutine cgyro_init_manager

  use mpi
  use timer_lib

  use cgyro_globals
  use GEO_interface

  implicit none

  include 'fftw3.f03'

  if (hiprec_flag == 1) then
     fmtstr    ='(es13.6)'
     fmtstr2   ='(2(es13.6,1x))'
     fmtstrn   ='(10(es13.6,1x))'
  endif

  !------------------------------------------------------
  ! Initialize startup timers 
  !  NOTE: All "runtime" timers are initialized 
  !        in cgyro_write_timedata
  !------------------------------------------------------
  call timer_lib_init('str_init')
  call timer_lib_init('coll_init')

  !----------------------------------------------------
  ! Initialize GLOBAL arrays
  !----------------------------------------------------

  allocate(energy(n_energy))
  allocate(w_e(n_energy))
  allocate(e_deriv1_mat(n_energy,n_energy))
  allocate(e_deriv2_mat(n_energy,n_energy))
  ! Construct energy nodes and weights
  if (e_method == 1) then
     call pseudo_maxwell(n_energy,nint(e_max),energy,w_e,e_deriv1_mat,e_deriv2_mat)
  else
     call pseudo_maxwell_new(n_energy,e_max,energy,w_e,e_deriv1_mat,e_deriv2_mat,trim(path)//'out.cgyro.egrid')
  endif

  allocate(xi(n_xi))
  allocate(w_xi(n_xi))
  allocate(xi_lor_mat(n_xi,n_xi))
  allocate(xi_deriv_mat(n_xi,n_xi))
  allocate(xi_upderiv_mat(n_xi,n_xi))
  ! Construct xi (pitch-angle) nodes and weights
  call pseudo_legendre(n_xi,xi,w_xi,xi_deriv_mat,xi_lor_mat,xi_upderiv_mat)

  allocate(theta(n_theta))
  allocate(thetab(n_radial/box_size,n_theta))
  allocate(w_theta(n_theta))
  allocate(bmag(n_theta))
  allocate(k_perp(n_theta,n_radial))
  allocate(omega_stream(n_theta,n_species))
  allocate(omega_trap(n_theta,n_species))
  allocate(omega_rdrift(n_theta,n_species))
  allocate(omega_adrift(n_theta,n_species))
  allocate(omega_aprdrift(n_theta,n_species))
  allocate(omega_cdrift(n_theta,n_species))
  allocate(omega_gammap(n_theta))

  if (test_flag == 0) then

     !----------------------------------------------------
     ! Initialize DISTRIBUTED arrays
     !----------------------------------------------------

     call timer_lib_in('str_init')

     ! Global (undistributed) arrays
     allocate(field(n_radial,n_theta,n_field))
     allocate(field_loc(n_radial,n_theta,n_field))
     allocate(field_old(n_radial,n_theta,n_field))
     allocate(field_old2(n_radial,n_theta,n_field))
     allocate(field_old3(n_radial,n_theta,n_field))
     allocate(moment_loc(n_radial,n_species))
     allocate(moment(n_radial,n_species))
     allocate(f_balloon(n_radial/box_size,n_theta))
     allocate(    flux(n_radial,n_species,2))
     allocate(flux_loc(n_radial,n_species,2))
     allocate(recv_status(MPI_STATUS_SIZE))

     allocate(thcyc(1-n_theta:2*n_theta))
     allocate(rcyc(n_radial,n_theta,-nup_theta:nup_theta))
     allocate(dtheta(n_radial,n_theta,-nup_theta:nup_theta))
     allocate(dtheta_up(n_radial,n_theta,-nup_theta:nup_theta))

     ! Velocity-distributed arrays
     allocate(rhs(4,nc,nv_loc))
     allocate(h_x(nc,nv_loc))
     allocate(g_x(nc,nv_loc))
     allocate(psi(nc,nv_loc))
     allocate(h0_x(nc,nv_loc))
     allocate(cap_h_c(nc,nv_loc))
     allocate(cap_h_ct(nv_loc,nc))
     allocate(omega_cap_h(nc,nv_loc))
     allocate(omega_h(nc,nv_loc))
     allocate(omega_s(n_field,nc,nv_loc))
     allocate(jvec_c(n_field,nc,nv_loc))
     allocate(jvec_v(n_field,nc_loc,nv))

     ! Real-space distributed arrays
     allocate(cap_h_v(nc_loc,nv))
     allocate(cap_h_v_prime(nc_loc,nv))

     ! Nonlinear arrays
     if (nonlinear_method == 1) then
        allocate(f_nl(nc,nsplit,n_toroidal))
        allocate(g_nl(nc,nsplit,n_toroidal))
     else
        allocate(f_nl(n_radial,nsplit,n_toroidal))
        allocate(g_nl(n_radial,nsplit,n_toroidal))
     endif

  endif

  ! Compute equilibrium quantities (even in test mode)
  GEO_model_in    = geo_numeq_flag
  GEO_ntheta_in   = geo_ntheta
  GEO_nfourier_in = geo_ny
  call GEO_alloc(1)
  call cgyro_equilibrium

  if (test_flag == 0) then

     call cgyro_init_arrays
     call cgyro_init_implicit_gk

     call timer_lib_out('str_init')

     call timer_lib_in('coll_init')

     allocate(cmat(nv,nv,nc_loc))
     allocate(cvec(nv))
     allocate(bvec(nv))

     call cgyro_init_collision

     call timer_lib_out('coll_init')

  endif

  ! Write initial data
  call cgyro_write_initdata

  if (test_flag == 1) then
     call MPI_FINALIZE(i_err)
     stop
  endif

  ! Initialize h (via restart or analytic IC)
  call timer_lib_in('str_init')
  call cgyro_init_h
  call timer_lib_out('str_init')

  ! Initialize nonlinear dimensions and arrays 
  if (nonlinear_method == 1) then

     ny0 = n_toroidal-1
     nx0 = n_radial/2
     ny = int(1.5*ny0)+1
     nx = int(1.5*nx0)+1

  else
     ! 2D FFT lengths 
     nx0 = n_radial
     ny0 = 2*n_toroidal-1

     ! 3/2-rule for dealiasing the nonlinear product
     nx = (3*nx0)/2
     ny = (3*ny0)/2
     ! Allocate and deallocate these every time.
     allocate(fx(0:ny/2,0:nx-1))
     allocate(gx(0:ny/2,0:nx-1))
     allocate(fy(0:ny/2,0:nx-1))
     allocate(gy(0:ny/2,0:nx-1))

     allocate(ux(0:ny-1,0:nx-1))
     allocate(vx(0:ny-1,0:nx-1))
     allocate(uy(0:ny-1,0:nx-1))
     allocate(vy(0:ny-1,0:nx-1))
     allocate(uv(0:ny-1,0:nx-1))

     ! Create plans once and for all, with global arrays fx,ux
     plan_c2r = fftw_plan_dft_c2r_2d(nx,ny,fx,ux,FFTW_MEASURE)
     plan_r2c = fftw_plan_dft_r2c_2d(nx,ny,ux,fx,FFTW_MEASURE)

  endif

  call GEO_alloc(0)

!$acc enter data copyin(energy,xi)
end subroutine cgyro_init_manager
