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
  use half_hermite

  use cgyro_io
#ifdef _OPENACC
  use cufft, only : cufftPlanMany, &
       CUFFT_C2R,CUFFT_Z2D,CUFFT_R2C,CUFFT_D2Z
#endif

  implicit none

#ifndef _OPENACC
  include 'fftw3.f03'
#endif

#ifdef _OPENACC
  integer :: howmany,istatus
  integer, parameter :: irank = 2
  integer, dimension(irank) :: ndim,inembed,onembed
  integer :: idist,odist,istride,ostride
  integer, parameter :: singlePrecision = selected_real_kind(6,30)
#endif

  character(len=128) :: msg

  if (hiprec_flag == 1) then
     fmtstr  = '(es16.9)'
     fmtstr_len = 17
     fmtstrn = '(10(es16.9,1x))'
  endif

  !------------------------------------------------------
  ! Initialize startup timers 
  !  NOTE: "Runtime" timers are initialized in cgyro_write_timedata,
  !  and input timer is initialized by read_input.
  !------------------------------------------------------
  call timer_lib_init('str_init')
  call timer_lib_init('nl_init')
  call timer_lib_init('coll_init')
  call timer_lib_init('io_init')

  !----------------------------------------------------
  ! Initialize GLOBAL arrays
  !----------------------------------------------------

  call timer_lib_in('str_init')
  allocate(energy(n_energy))
  allocate(vel(n_energy))
  allocate(w_e(n_energy))
  allocate(e_deriv1_mat(n_energy,n_energy))
  allocate(e_deriv1_rot_mat(n_energy,n_energy))

  ! Construct energy nodes and weights
  if (e_method<=2) then
     call pseudo_maxwell_new(n_energy,&
          e_max,&
          energy,&
          w_e,&
          e_deriv1_mat,&
          trim(path)//'out.cgyro.egrid')
  else if (e_method==3) then
     ! interface function in module half_hermite
     if (i_proc==0) then
        call pseudo_maxwell_pliocene(n_energy,&
             e_max,&
             energy,&
             w_e,&
             e_deriv1_mat,&
             alpha_poly,& ! weight fct=x^alpha_poly*exp(-x**2)
             trim(path)//'out.cgyro.egrid')
     else
        call pseudo_maxwell_pliocene(n_energy,&
             e_max,&
             energy,&
             w_e,&
             e_deriv1_mat,&
             alpha_poly) ! only write results on i_proc zero.
     end if
  end if
     

  vel(:) = sqrt(energy(:))

  e_deriv1_rot_mat(:,:) = e_deriv1_mat(:,:)
  if(e_fix == 2) then
     e_deriv1_rot_mat(n_energy,:) = 0.0
  endif
  if(e_fix == 3) then
     e_deriv1_rot_mat(n_energy,:) = 0.0
     e_deriv1_mat(n_energy,:) = 0.0
  endif
  
  allocate(xi(n_xi))
  allocate(w_xi(n_xi))
  allocate(xi_lor_mat(n_xi,n_xi))
  allocate(xi_deriv_mat(n_xi,n_xi))
  ! Construct xi (pitch-angle) nodes and weights
  call pseudo_legendre(n_xi,xi,w_xi,xi_deriv_mat,xi_lor_mat)
  w_xi = 0.5*w_xi

  allocate(theta(n_theta))
  allocate(thetab(n_theta,n_radial/box_size))
  allocate(w_theta(n_theta))
  allocate(g_theta(n_theta))
  allocate(g_theta_geo(n_theta))
  allocate(bmag(n_theta))
  allocate(btor(n_theta))
  allocate(bpol(n_theta))
  allocate(k_perp(nc))
  allocate(k_x(nc))
  allocate(bigr(n_theta))
  allocate(bigr_r(n_theta))
  allocate(itp(n_theta))
  allocate(omega_stream(n_theta,n_species))
  allocate(omega_trap(n_theta,n_species))
  allocate(omega_rdrift(n_theta,n_species))
  allocate(omega_adrift(n_theta,n_species))
  allocate(omega_aprdrift(n_theta,n_species))
  allocate(omega_cdrift(n_theta,n_species))
  allocate(omega_cdrift_r(n_theta,n_species))
  allocate(omega_gammap(n_theta))

  allocate(lambda_rot(n_theta,n_species))
  allocate(dlambda_rot(n_theta,n_species))
  allocate(dens_rot(n_theta,n_species))
  allocate(dens_ele_rot(n_theta))
  allocate(dens_avg_rot(n_species))
  allocate(dlnndr_avg_rot(n_species))
  allocate(omega_rot_trap(n_theta,n_species))
  allocate(omega_rot_u(n_theta,n_species))
  allocate(omega_rot_drift(n_theta,n_species))
  allocate(omega_rot_drift_r(n_theta,n_species))
  allocate(omega_rot_edrift(n_theta))
  allocate(omega_rot_edrift_r(n_theta))
  allocate(omega_rot_star(n_theta,n_species))

  if (test_flag == 0) then

     !----------------------------------------------------
     ! Initialize DISTRIBUTED arrays
     !----------------------------------------------------

     ! Global (undistributed) arrays
     allocate(fcoef(n_field,nc))
     if (n_field < 3) then
        allocate(gcoef(n_field,nc))
     else
        allocate(gcoef(5,nc))
     endif
     allocate(field(n_field,nc))
     allocate(field_dot(n_field,nc))
     allocate(field_loc(n_field,nc))
     allocate(field_old(n_field,nc))
     allocate(field_old2(n_field,nc))
     allocate(field_old3(n_field,nc))
     allocate(    moment(n_radial,theta_plot,n_species,3))
     allocate(moment_loc(n_radial,theta_plot,n_species,3))
     allocate(    cflux(n_species,4,n_field))
     allocate(cflux_loc(n_species,4,n_field))
     allocate(    gflux(0:n_global,n_species,4,n_field))
     allocate(gflux_loc(0:n_global,n_species,4,n_field))
     allocate(cflux_tave(n_species,4))
     allocate(gflux_tave(n_species,4))
     
     allocate(recv_status(MPI_STATUS_SIZE))

     allocate(icd_c(-nup_theta:nup_theta, nc))
     allocate(dtheta(-nup_theta:nup_theta, nc))
     allocate(dtheta_up(-nup_theta:nup_theta, nc))
     allocate(source(n_theta,nv_loc))

!$acc enter data create(fcoef,gcoef,field,field_loc,source)

     ! Velocity-distributed arrays

     select case(delta_t_method)
     case(1)
        allocate(h0_old(nc,nv_loc))     
        allocate(rhs(nc,nv_loc,6))
!$acc enter data create(rhs,h0_old)
     case(2)
        allocate(h0_old(nc,nv_loc))
        allocate(rhs(nc,nv_loc,7))
!$acc enter data create(rhs,h0_old)
     case(3)
        allocate(h0_old(nc,nv_loc))
        allocate(rhs(nc,nv_loc,9))
!$acc enter data create(rhs,h0_old)
     case default
        ! Normal timestep
        allocate(rhs(nc,nv_loc,4))
!$acc enter data create(rhs)
     end select 
     
     allocate(h_x(nc,nv_loc))
     allocate(g_x(nc,nv_loc))
     allocate(h0_x(nc,nv_loc))
!$acc enter data create(h_x,g_x,h0_x)

     allocate(cap_h_c(nc,nv_loc))
     allocate(cap_h_ct(nv_loc,nc))
     allocate(cap_h_c_dot(nc,nv_loc))
     allocate(cap_h_c_old(nc,nv_loc))
     allocate(cap_h_c_old2(nc,nv_loc))
     allocate(cap_h_v(nc_loc,nv))
     allocate(omega_cap_h(nc,nv_loc))
     allocate(omega_h(nc,nv_loc))
     allocate(omega_s(n_field,nc,nv_loc))
     allocate(omega_ss(n_field,nc,nv_loc))
     allocate(jvec_c(n_field,nc,nv_loc))
     allocate(jvec_v(n_field,nc_loc,nv))
     allocate(dvjvec_c(n_field,nc,nv_loc))
     allocate(dvjvec_v(n_field,nc_loc,nv))
     allocate(jxvec_c(n_field,nc,nv_loc))
     allocate(upfac1(nc,nv_loc,2))
     allocate(upfac2(nc,nv_loc,2))

!$acc enter data create(cap_h_c,cap_h_ct,cap_h_v,dvjvec_c,dvjvec_v)

     if (upwind_single_flag == 0) then
       allocate(upwind_res_loc(nc,ns1:ns2,2))
       allocate(upwind_res(nc,ns1:ns2,2))
!$acc enter data create(upwind_res,upwind_res_loc)
     else
       allocate(upwind32_res_loc(nc,ns1:ns2,2))
       allocate(upwind32_res(nc,ns1:ns2,2))
!$acc enter data create(upwind32_res,upwind32_res_loc)
     endif

     ! Nonlinear arrays
     if (nonlinear_flag == 1) then
        if (nonlinear_method == 1) then
           call cgyro_error("nonlinear_method==1 has been deprecated")
           return
        else
           allocate(f_nl(n_radial,nsplit,n_toroidal))
           allocate(g_nl(n_field,n_radial,n_jtheta,n_toroidal))
           allocate(fpack(n_radial,nsplit*n_toroidal))
           allocate(gpack(n_field,n_radial,n_jtheta,n_toroidal))
        endif
        allocate(jvec_c_nl(n_field,n_radial,n_jtheta,nv_loc,n_toroidal))
!$acc enter data create(fpack,gpack,f_nl,g_nl,jvec_c_nl)
     endif

     if (collision_model == 5) then
        allocate(cmat_simple(n_xi,n_xi,n_energy,n_species,n_theta))
     else
        if (collision_precision_mode /= 0) then
           allocate(cmat_stripes(-collision_full_stripes:collision_full_stripes,nv,nc_loc))
           allocate(cmat_fp32(nv,nv,nc_loc))

           write (msg, "(A35,I4,A14)") "Using fp32 collision precision with ",collision_full_stripes, " fp64 stripes."
           call cgyro_info(msg)
        else
           allocate(cmat(nv,nv,nc_loc))
        endif
     endif

  endif

  call cgyro_equilibrium

#ifndef _OPENACC
  gpu_bigmem_flag = 0
#endif

  if (test_flag == 0) then

     call cgyro_init_arrays
     call timer_lib_out('str_init')

     call timer_lib_in('coll_init')
     call cgyro_init_collision
     call timer_lib_out('coll_init')

     call timer_lib_in('str_init')
  endif

  call cgyro_check_memory(trim(path)//runfile_memory)

  if (velocity_order == 1) then
    ! traditional ordering
    restart_magic = 140906808
  else
    ! alternative ordering, need different magic
    restart_magic = 140916753
  endif

  call timer_lib_out('str_init')

  ! Write initial data

  call timer_lib_in('io_init')
  call cgyro_write_initdata
  call timer_lib_out('io_init')

  if (test_flag == 1) return

  ! Initialize h (via restart or analytic IC)
  call timer_lib_in('str_init')
  call cgyro_init_h
  if (error_status /=0 ) return
  call timer_lib_out('str_init')

  ! Initialize nonlinear dimensions and arrays 
  call timer_lib_in('nl_init')
  if (nonlinear_method == 1) then
     call cgyro_error("nonlinear_method==1 has been deprecated")
     return
  else

     ! 2D FFT lengths 
     nx0 = n_radial
     ny0 = 2*n_toroidal-1

     ! 3/2-rule for dealiasing the nonlinear product
     nx = (3*nx0)/2
     ny = (3*ny0)/2
     
#ifndef _OPENACC
     allocate(fx(0:ny/2,0:nx-1,n_omp))
     allocate(gx(0:ny/2,0:nx-1,n_omp))
     allocate(fy(0:ny/2,0:nx-1,n_omp))
     allocate(gy(0:ny/2,0:nx-1,n_omp))

     allocate(uxmany(0:ny-1,0:nx-1,nsplit))
     allocate(uymany(0:ny-1,0:nx-1,nsplit))
     allocate(vx(0:ny-1,0:nx-1,n_omp))
     allocate(vy(0:ny-1,0:nx-1,n_omp))
     allocate(uv(0:ny-1,0:nx-1,n_omp))

#ifdef THREADED_FFT
     i_err = fftw_init_threads()
     call fftw_plan_with_nthreads(n_omp)
#endif


     ! Create plans once and for all, with global arrays fx,ux
     plan_c2r = fftw_plan_dft_c2r_2d(nx,ny,gx(:,:,1),vx(:,:,1),FFTW_PATIENT)
     plan_r2c = fftw_plan_dft_r2c_2d(nx,ny,uv(:,:,1),fx(:,:,1),FFTW_PATIENT)
#endif

#ifdef _OPENACC
     call cgyro_info('GPU-aware code triggered.')

     allocate( fxmany(0:ny/2,0:nx-1,nsplit) )
     allocate( fymany(0:ny/2,0:nx-1,nsplit) )
     allocate( gxmany(0:ny/2,0:nx-1,nsplit) )
     allocate( gymany(0:ny/2,0:nx-1,nsplit) )

     allocate( uxmany(0:ny-1,0:nx-1,nsplit) )
     allocate( uymany(0:ny-1,0:nx-1,nsplit) )
     allocate( vxmany(0:ny-1,0:nx-1,nsplit) )
     allocate( vymany(0:ny-1,0:nx-1,nsplit) )
     allocate( uvmany(0:ny-1,0:nx-1,nsplit) )

!$acc enter data create(fxmany,fymany,gxmany,gymany) &
!$acc&           create(uxmany,uymany,vxmany,vymany,uvmany)

     !-------------------------------------------------------------------
     ! 2D
     !   input[ b*idist + (x * inembed[1] + y)*istride ]
     !  output[ b*odist + (x * onembed[1] + y)*ostride ]
     !  isign is the sign of the exponent in the formula that defines
     !  Fourier transform  -1 == FFTW_FORWARD
     !                      1 == FFTW_BACKWARD
     !-------------------------------------------------------------------

     ndim(1) = nx
     ndim(2) = ny
     idist = size(fxmany,1)*size(fxmany,2)
     odist = size(uxmany,1)*size(uxmany,2)
     istride = 1
     ostride = 1
     inembed = size(fxmany,1)
     onembed = size(uxmany,1)

     istatus = cufftPlanMany(&
          cu_plan_c2r_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          merge(CUFFT_C2R,CUFFT_Z2D,kind(uxmany) == singlePrecision), &
          nsplit)

     idist = size(uxmany,1)*size(uxmany,2)
     odist = size(fxmany,1)*size(fxmany,2)
     inembed = size(uxmany,1)
     onembed = size(fxmany,1) 
     istride = 1
     ostride = 1
     istatus = cufftPlanMany(&
          cu_plan_r2c_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          merge(CUFFT_R2C,CUFFT_D2Z,kind(uxmany) == singlePrecision), &
          nsplit)
#endif

  endif
  call timer_lib_out('nl_init')

end subroutine cgyro_init_manager
