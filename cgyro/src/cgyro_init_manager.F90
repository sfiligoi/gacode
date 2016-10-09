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

#ifdef _OPENACC
  use precision_m, only : singlePrecision
  use cufft_m, only : cufftPlanMany, &
       CUFFT_C2R,CUFFT_Z2D,CUFFT_R2C,CUFFT_D2Z
#endif
  implicit none

  real :: b

  include 'fftw3.f03'

#ifdef _OPENACC
  integer :: howmany,istatus
  integer, parameter :: irank = 2
  integer, dimension(irank) :: ndim,inembed,onembed
  integer :: idist,odist,istride,ostride

#endif

  if (hiprec_flag == 1) then
     fmtstr    ='(es16.9)'
     fmtstr2   ='(2(es16.9,1x))'
     fmtstrn   ='(10(es16.9,1x))'
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
  allocate(vel(n_energy))
  allocate(w_e(n_energy))
  allocate(e_deriv1_mat(n_energy,n_energy))
  allocate(e_deriv2_mat(n_energy,n_energy))
  ! Construct energy nodes and weights
  if (e_method == 1) then
     call pseudo_maxwell(n_energy,&
          nint(e_max),&
          energy,&
          w_e,&
          e_deriv1_mat,&
          e_deriv2_mat)
  else
     call pseudo_maxwell_new(n_energy,&
          e_max,&
          energy,&
          w_e,&
          e_deriv1_mat,&
          e_deriv2_mat,&
          trim(path)//'out.cgyro.egrid')
  endif
  vel(:) = sqrt(energy(:))

  !----------------------------------------------------------------------
  ! Correction factor for missing energy interval to ensure sum(w)=1.0
  ! NOTE: without this we have poor grid-convergence for small e_max
  !
  b = sqrt(e_max)
  !w_e(n_energy) = w_e(n_energy)+2.0*exp(-e_max)*b/sqrt(pi)+erfc(b)
  !----------------------------------------------------------------------

  allocate(xi(n_xi))
  allocate(w_xi(n_xi))
  allocate(xi_lor_mat(n_xi,n_xi))
  allocate(xi_deriv_mat(n_xi,n_xi))
  ! Construct xi (pitch-angle) nodes and weights
  call pseudo_legendre(n_xi,xi,w_xi,xi_deriv_mat,xi_lor_mat)
  w_xi = 0.5*w_xi

  allocate(theta(n_theta))
  allocate(thetab(n_radial/box_size,n_theta))
  allocate(w_theta(n_theta))
  allocate(bmag(n_theta))
  allocate(k_perp(nc))
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
     allocate(fcoef(n_field,nc))
     if (n_field < 3) then
        allocate(gcoef(n_field,nc))
     else
        allocate(gcoef(5,nc))
     endif
     allocate(field(n_field,nc))
     allocate(field_loc(n_field,nc))
     allocate(field_old(n_field,nc))
     allocate(field_old2(n_field,nc))
     allocate(field_old3(n_field,nc))
     allocate(    moment(n_radial,n_species,2))
     allocate(moment_loc(n_radial,n_species,2))
     allocate(    flux(n_radial,n_species,2))
     allocate(flux_loc(n_radial,n_species,2))
     allocate(f_balloon(n_radial/box_size,n_theta))
     allocate(recv_status(MPI_STATUS_SIZE))

     allocate(icd_c(nc,-nup_theta:nup_theta))
     allocate(dtheta(nc,-nup_theta:nup_theta))
     allocate(dtheta_up(nc,-nup_theta:nup_theta))

     ! Velocity-distributed arrays
     allocate(rhs(nc,nv_loc,4))
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
     if (nonlinear_flag == 1) then
        if (nonlinear_method == 1) then
           allocate(f_nl(nc,nsplit,n_toroidal))
           allocate(g_nl(nc,nsplit,n_toroidal))
        else
           allocate(f_nl(n_radial,nsplit,n_toroidal))
           allocate(g_nl(n_radial,nsplit,n_toroidal))
        endif
     endif

     if (collision_model == 5) then
        allocate(cmat_simple(n_xi,n_xi,n_energy,n_species,n_theta))
     else
        allocate(cmat(nv,nv,nc_loc))
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

     call cgyro_init_collision

     call timer_lib_out('coll_init')

  endif

  call cgyro_check_memory(trim(path)//runfile_memory)

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

     ! Direct convolution

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
     plan_c2r = fftw_plan_dft_c2r_2d(nx,ny,fx,ux,FFTW_PATIENT)
     plan_r2c = fftw_plan_dft_r2c_2d(nx,ny,ux,fx,FFTW_PATIENT)

#ifdef _OPENACC
     howmany = nsplit
     allocate( fxmany(0:ny/2,0:nx-1,howmany) )
     allocate( fymany(0:ny/2,0:nx-1,howmany) )
     allocate( gxmany(0:ny/2,0:nx-1,howmany) )
     allocate( gymany(0:ny/2,0:nx-1,howmany) )

     allocate( uxmany(0:ny-1,0:nx-1,howmany) )
     allocate( uymany(0:ny-1,0:nx-1,howmany) )
     allocate( vxmany(0:ny-1,0:nx-1,howmany) )
     allocate( vymany(0:ny-1,0:nx-1,howmany) )
     allocate( uvmany(0:ny-1,0:nx-1,howmany) )
!!$acc enter data create(fxmany,fymany,gxmany,gymany)
!!$acc enter data create(uxmany,uymany,vxmany,vymany)
!!$acc enter data create(uvmany)


     !   -------------------------------------
     ! 2D
     !   input[ b*idist + (x * inembed[1] + y) * istride ]
     !   output[ b*odist + (x * onembed[1] + y)*ostride ]
     !   isign is the sign of the exponent in the formula that defines
     !   Fourier transform  -1 == FFTW_FORWARD
     !                       1 == FFTW_BACKWARD
     !   -------------------------------------

     ndim(1) = nx
     ndim(2) = ny
     idist = size( fxmany,1)*size(fxmany,2)
     odist = size( uxmany,1)*size(uxmany,2)
     istride = 1
     ostride = 1
     inembed = size(fxmany,1)
     onembed = size(uxmany,1)

     istatus = cufftPlanMany( cu_plan_c2r_many,                         &
          &                                    irank,                         &
          &                                    ndim,                          &
          &                                    inembed,                       &
          &                                    istride,                       &
          &                                    idist,                         &
          &                                    onembed,                       &
          &                                    ostride,                       &
          &                                    odist,                         &
          &  merge(CUFFT_C2R,CUFFT_Z2D,kind(uxmany).eq.singlePrecision),      &
          &                                    howmany )

     idist = size(uxmany,1)*size(uxmany,2)
     odist = size(fxmany,1)*size(fxmany,2)
     inembed = size(uxmany,1)
     onembed = size(fxmany,1) 
     istride = 1
     ostride = 1
     istatus = cufftPlanMany( cu_plan_r2c_many,                         &
          &                     irank,                                        &
          &                     ndim,                                         &
          &                     inembed,                                      &
          &                     istride,                                      &
          &                     idist,                                        &
          &                     onembed,                                      &
          &                     ostride,                                      &
          &                     odist,                                        &
          & merge(CUFFT_R2C,CUFFT_D2Z,kind(uxmany).eq.singlePrecision),       &
          &                     howmany )
#endif

  endif

end subroutine cgyro_init_manager
