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
  use cgyro_nl

#if defined(_OPENACC) || defined(OMPGPU)
#define CGYRO_GPU_FFT
#endif

  implicit none

  character(len=128) :: msg
  integer :: ie,ix

  if (hiprec_flag == 1) then
     BYTE   = 8
     fmtstr = '(es16.9)'
     fmtstr_len = 17
     fmtstrn = '(10(es16.9,1x))'
  else
     BYTE   = 4
     fmtstr ='(es11.4)'
     fmtstr_len = 12
     fmtstrn ='(10(es11.4,1x))'
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
  allocate(vel2(n_energy))
  allocate(w_e(n_energy))
  allocate(e_deriv1_mat(n_energy,n_energy))
  allocate(e_deriv1_rot_mat(n_energy,n_energy))

  ! Construct energy nodes and weights (Hallatschek)
  if (i_proc == 0) then
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
  endif
  
  ! Default value for collision data compression
  n_low_energy = 0

  vel(:) = sqrt(energy(:))
  vel2(:) = sqrt(2.0*energy(:))

#if defined(OMPGPU)
!$omp target enter data map(to:vel,vel2)
#elif defined(_OPENACC)
!$acc enter data copyin(vel,vel2)
#endif

  e_deriv1_rot_mat(:,:) = e_deriv1_mat(:,:)
  if (e_fix == 2) then
     e_deriv1_rot_mat(n_energy,:) = 0.0
  endif
  if (e_fix == 3) then
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

  allocate(w_exi(n_energy,n_xi))
  do ix=1,n_xi
     do ie=1,n_energy
        w_exi(ie,ix) = w_e(ie)*w_xi(ix)
     enddo
  enddo

  allocate(theta(n_theta))
  allocate(thetab(n_theta,n_radial/box_size))
  allocate(w_theta(n_theta))
  allocate(g_theta(n_theta))
  allocate(g_theta_geo(n_theta))
  allocate(bmag(n_theta))
  allocate(btor(n_theta))
  allocate(bpol(n_theta))
  allocate(k_perp(nc,nt1:nt2))
  allocate(k_x(nc,nt1:nt2))
  allocate(bigr(n_theta))
  allocate(bigr_r(n_theta))
  allocate(captheta(n_theta))
  allocate(itp(n_theta))
  allocate(omega_stream(n_theta,n_species,nt1:nt2))
  allocate(omega_trap(n_theta,n_species,nt1:nt2))
  allocate(omega_rdrift(n_theta,n_species))
  allocate(omega_adrift(n_theta,n_species))
  allocate(omega_aprdrift(n_theta,n_species))
  allocate(omega_cdrift(n_theta,n_species))
  allocate(omega_cdrift_r(n_theta,n_species))
  allocate(omega_gammap(n_theta))

  allocate(lambda_rot(n_theta,n_species))
  allocate(dlambda_rot(n_theta,n_species))
  allocate(dens_rot(n_theta,n_species))
  allocate(dens2_rot(n_theta,n_species))
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

  allocate(gtime(nt1:nt2))
  allocate(freq(nt1:nt2))
  freq_err = 0

  if (test_flag == 0) then

     !----------------------------------------------------
     ! Initialize DISTRIBUTED arrays
     !----------------------------------------------------

     ! Global (undistributed) arrays
     allocate(fcoef(n_field,nc,nt1:nt2))
     if (n_field < 3) then
        allocate(gcoef(n_field,nc,nt1:nt2))
     else
        allocate(gcoef(5,nc,nt1:nt2))
     endif
     allocate(field(n_field,nc,nt1:nt2))
     allocate(field_dot(n_field,nc,nt1:nt2))
     allocate(field_loc(n_field,nc,nt1:nt2))
     allocate(field_old(n_field,nc,nt1:nt2))
     allocate(field_old2(n_field,nc,nt1:nt2))
     allocate(field_old3(n_field,nc,nt1:nt2))
     allocate(    moment(n_radial,theta_plot,n_species,nt1:nt2,3))
     allocate(moment_loc(n_radial,theta_plot,n_species,nt1:nt2,3))
     allocate(    cflux(n_species,4,n_field,nt1:nt2))
     allocate(cflux_loc(n_species,4,n_field,nt1:nt2))
     allocate(    gflux(0:n_global,n_species,4,n_field,nt1:nt2))
     allocate(gflux_loc(0:n_global,n_species,4,n_field,nt1:nt2))
     allocate(cflux_tave(n_species,4))
     allocate(gflux_tave(n_species,4))

     allocate(recv_status(MPI_STATUS_SIZE))

     allocate(source(n_theta,nv_loc,nt1:nt2))

#if defined(OMPGPU)
!$omp target enter data map(alloc:fcoef,gcoef,field,field_loc,source)
#elif defined(_OPENACC)
!$acc enter data create(fcoef,gcoef,field,field_loc,source)
#endif

     if ((collision_model /= 5) .AND. (collision_field_model == 1)) then
       ! assuming all collision constants are the same between the simulations
       allocate(dvjvec_v(n_field,nv,nt1:nt2,nc_loc_coll))
       ! we do not really need all n_sim
       ! but since n_sim is assumed to be small, the added cost is small
       ! and this drastically similifies the code
       ! But could be improved in the future
       allocate(jvec_v(n_field,nc_loc_coll,nt1:nt2,nv,n_sim))
       ! nc and nc_loc_coll must be last, since it will be collated     
       allocate(field_v(n_field,nt1:nt2,n_sim,nc))
       allocate(field_loc_v(n_field,nt1:nt2,n_sim,nc_cl1:nc_cl2))
#if defined(OMPGPU)
!$omp target enter data map(alloc:field_v,field_loc_v,dvjvec_v)
#elif defined(_OPENACC)
!$acc enter data create(field_v,field_loc_v,dvjvec_v)
#endif
     endif

     ! Velocity-distributed arrays

     select case(delta_t_method)
     case(1)
        allocate(h0_old(nc,nv_loc,nt1:nt2))
        allocate(rhs(nc,nv_loc,nt1:nt2,6))
#if defined(OMPGPU)
!$omp target enter data map(alloc:rhs,h0_old)
#elif defined(_OPENACC)
!$acc enter data create(rhs,h0_old)
#endif
     case(2)
        allocate(h0_old(nc,nv_loc,nt1:nt2))
        allocate(rhs(nc,nv_loc,nt1:nt2,7))
#if defined(OMPGPU)
!$omp target enter data map(alloc:rhs,h0_old)
#elif defined(_OPENACC)
!$acc enter data create(rhs,h0_old)
#endif
     case(3)
        allocate(h0_old(nc,nv_loc,nt1:nt2))
        allocate(rhs(nc,nv_loc,nt1:nt2,9))
#if defined(OMPGPU)
!$omp target enter data map(alloc:rhs,h0_old)
#elif defined(_OPENACC)
!$acc enter data create(rhs,h0_old)
#endif
     case default
        ! Normal timestep
        allocate(rhs(nc,nv_loc,nt1:nt2,4))
#if defined(OMPGPU)
!$omp target enter data map(alloc:rhs)
#elif defined(_OPENACC)
!$acc enter data create(rhs)
#endif
     end select 
     
     allocate(h_x(nc,nv_loc,nt1:nt2))
     allocate(g_x(nc,nv_loc,nt1:nt2))
     allocate(h0_x(nc,nv_loc,nt1:nt2))
#if defined(OMPGPU)
!$omp target enter data map(alloc:h_x,g_x,h0_x)
#elif defined(_OPENACC)
!$acc enter data create(h_x,g_x,h0_x)
#endif

     allocate(cap_h_c(nc,nv_loc,nt1:nt2))
     allocate(cap_h_c_dot(nc,nv_loc,nt1:nt2))
     allocate(cap_h_c_old(nc,nv_loc,nt1:nt2))
     allocate(cap_h_c_old2(nc,nv_loc,nt1:nt2))
     allocate(cap_h_ct(nv_loc,nt1:nt2,nc))
     allocate(cap_h_v(nc_loc_coll,nt1:nt2,nv,n_sim))
     allocate(omega_cap_h(nc,nv_loc,nt1:nt2))
     allocate(omega_h(nc,nv_loc,nt1:nt2))
     allocate(omega_s(n_field,nc,nv_loc,nt1:nt2))
     allocate(omega_ss(n_field,nc,nv_loc,nt1:nt2))
     allocate(omega_sbeta(nc,nv_loc,nt1:nt2))
     allocate(jvec_c(n_field,nc,nv_loc,nt1:nt2))
     allocate(dvjvec_c(n_field,nc,nv_loc,nt1:nt2))
     allocate(jxvec_c(n_field,nc,nv_loc,nt1:nt2))
     allocate(upfac1(nc,nv_loc,nt1:nt2))
     allocate(upfac2(nc,nv_loc,nt1:nt2))
     
#if defined(OMPGPU)
!$omp target enter data map(alloc:cap_h_c,cap_h_ct,cap_h_c_dot,cap_h_c_old,cap_h_c_old2)
!$omp target enter data map(alloc:cap_h_v,dvjvec_c)
#elif defined(_OPENACC)
!$acc enter data create(cap_h_c,cap_h_ct,cap_h_c_dot,cap_h_c_old,cap_h_c_old2)
!$acc enter data create(cap_h_v,dvjvec_c)
#endif

     if (upwind_single_flag == 0) then
       allocate(upwind_res_loc(nc,ns1:ns2,nt1:nt2))
       allocate(upwind_res(nc,ns1:ns2,nt1:nt2))
#if defined(OMPGPU)
!$omp target enter data map(alloc:upwind_res,upwind_res_loc)
#elif defined(_OPENACC)
!$acc enter data create(upwind_res,upwind_res_loc)
#endif
     else
       allocate(upwind32_res_loc(nc,ns1:ns2,nt1:nt2))
       allocate(upwind32_res(nc,ns1:ns2,nt1:nt2))
#if defined(OMPGPU)
!$omp target enter data map(alloc:upwind32_res,upwind32_res_loc)
#elif defined(_OPENACC)
!$acc enter data create(upwind32_res,upwind32_res_loc)
#endif
     endif

     ! Nonlinear arrays
     if (nonlinear_flag == 1) then
        if (nl_single_flag < 2) then
          allocate(fA_nl(n_radial,nt_loc,nsplitA,n_toroidal_procs))
          allocate(g_nl(n_field,n_radial,n_jtheta,n_toroidal))
          allocate(fpackA(n_radial,nt_loc,nsplitA*n_toroidal_procs))
          allocate(gpack(n_field,n_radial,n_jtheta,n_toroidal))
          allocate(jvec_c_nl(n_field,n_radial,n_jtheta,nv_loc,n_toroidal))
#if defined(OMPGPU)
!$omp target enter data map(alloc:fpackA,gpack,fA_nl,g_nl,jvec_c_nl)
#elif defined(_OPENACC)
!$acc enter data create(fpackA,gpack,fA_nl,g_nl,jvec_c_nl)
#endif
        endif
        if (nl_single_flag > 1) then
          allocate(g_nl32(n_field,n_radial,n_jtheta,n_toroidal))
          allocate(gpack32(n_field,n_radial,n_jtheta,n_toroidal))
          allocate(jvec_c_nl32(n_field,n_radial,n_jtheta,nv_loc,n_toroidal))
#if defined(OMPGPU)
!$omp target enter data map(alloc:gpack32,g_nl32,jvec_c_nl32)
#elif defined(_OPENACC)
!$acc enter data create(gpack32,g_nl32,jvec_c_nl32)
#endif
        endif
        if (nl_single_flag .NE. 0) then
          allocate(fA_nl32(n_radial,nt_loc,nsplitA,n_toroidal_procs))
          allocate(fpackA32(n_radial,nt_loc,nsplitA*n_toroidal_procs))
#if defined(OMPGPU)
!$omp target enter data map(alloc:fpackA32,fA_nl32)
#elif defined(_OPENACC)
!$acc enter data create(fpackA32,fA_nl32)
#endif
        endif
        if (nsplitB > 0) then ! nsplitB can be zero at large MPI
          allocate(fB_nl(n_radial,nt_loc,nsplitB,n_toroidal_procs))
          allocate(fpackB(n_radial,nt_loc,nsplitB*n_toroidal_procs))
#if defined(OMPGPU)
!$omp target enter data map(alloc:fpackB,fB_nl)
#elif defined(_OPENACC)
!$acc enter data create(fpackB,fB_nl)
#endif
          if (nl_single_flag .NE. 0) then
          allocate(fB_nl32(n_radial,nt_loc,nsplitB,n_toroidal_procs))
          allocate(fpackB32(n_radial,nt_loc,nsplitB*n_toroidal_procs))
#if defined(OMPGPU)
!$omp target enter data map(alloc:fpackB32,fB_nl32)
#elif defined(_OPENACC)
!$acc enter data create(fpackB32,fB_nl32)
#endif
          endif
        endif
     endif

     if (collision_model == 5) then
        if (nc_loc_coll/=nc_loc) then
           call cgyro_error("CMAT sharing not supported for COLLISION_MODEL 5")
           return
        endif
        allocate(cmat_simple(n_xi,n_xi,n_energy,n_species,n_theta,nt1:nt2))
     else
        if (collision_precision_mode == 1) then
           ! the lowest energy(s) has the most spread, so treat differently
           n_low_energy = 1
           do ie=2,n_energy
             if (energy(ie)<1.0e-2) then
               n_low_energy = ie
             endif
           enddo
           allocate(cmat_fp32(nv,nv,nc_loc_coll,nt1:nt2))
           allocate(cmat_stripes(n_xi,n_species,(n_low_energy+1):n_energy,n_xi,nc_loc_coll,nt1:nt2))
           allocate(cmat_e1(n_xi,n_species,n_low_energy,nv,nc_loc_coll,nt1:nt2))

           write (msg, "(A,I1,A)") "Using fp32 collision precision except e<=",n_low_energy," or same e&s."
           call cgyro_info(msg)
        else if (collision_precision_mode == 32) then
           allocate(cmat_fp32(nv,nv,nc_loc_coll,nt1:nt2))
        else
           allocate(cmat(nv,nv,nc_loc_coll,nt1:nt2))
        endif
     endif

  endif

  call cgyro_equilibrium
  if (error_status > 0) then
     ! something went terribly wrong
     return
  endif

#ifndef CGYRO_GPU_FFT
  gpu_bigmem_flag = 0
#endif

  if (test_flag == 0) then

     call cgyro_init_arrays
     call timer_lib_out('str_init')
     if (error_status > 0) then
        ! something went terribly wrong
        return
     endif

     call timer_lib_in('coll_init')
     call cgyro_init_collision
     call timer_lib_out('coll_init')
     if (error_status > 0) then
        ! something went terribly wrong
        return
     endif

     call timer_lib_in('str_init')
  endif

  ! 2D FFT lengths 
  nx0 = n_radial
  ny0 = 2*n_toroidal-1

  ! +1 to properly handle odd n_radial
  nx2 = (nx0+1)/2

  ! 3/2-rule for dealiasing the nonlinear product
  nx = (3*nx0)/2
  ! old, obsolete defintion: ny = (3*ny0)/2
  ! new definiton: ny = (3*(ny0+1))/2
  ny = 3*n_toroidal

  ny2 = ny/2

  call cgyro_check_memory(trim(path)//runfile_memory)

  call timer_lib_out('str_init')

  ! Write initial data

  call timer_lib_in('io_init')
  call cgyro_write_initdata
  call timer_lib_out('io_init')

  if (test_flag == 1) return

  ! Initialize h (via restart or analytic IC)
  call timer_lib_in('str_init')
  call cgyro_init_h
  if (error_status > 0) return
  call timer_lib_out('str_init')

  call timer_lib_in('nl_init')

  if (nonlinear_flag == 1) then

  ! Initialize nonlinear dimensions and arrays 
#if defined(OMPGPU)
!$omp target enter data map(to:nx0,ny0,nx,ny,nx2,ny2)
#elif defined(_OPENACC)
!$acc enter data copyin(nx0,ny0,nx,ny,nx2,ny2)
#endif

#ifndef CGYRO_GPU_FFT

  if (nl_single_flag > 1) then
    ! Use 32-bit NL
    allocate(fx32(0:ny2,0:nx-1,n_omp))
    allocate(gx32(0:ny2,0:nx-1,n_omp))
    allocate(fy32(0:ny2,0:nx-1,n_omp))
    allocate(gy32(0:ny2,0:nx-1,n_omp))

    ! Note: Assuming nsplitA>=nsplitB
    !       So we can use the same buffers for both
    allocate(vxmany32(0:ny-1,0:nx-1,nsplit))
    allocate(vymany32(0:ny-1,0:nx-1,nsplit))
    allocate(uxmany32(0:ny-1,0:nx-1,nsplitA))
    allocate(uymany32(0:ny-1,0:nx-1,nsplitA))
    allocate(uv32(0:ny-1,0:nx-1,n_omp))
  else
    allocate(fx(0:ny2,0:nx-1,n_omp))
    allocate(gx(0:ny2,0:nx-1,n_omp))
    allocate(fy(0:ny2,0:nx-1,n_omp))
    allocate(gy(0:ny2,0:nx-1,n_omp))

    ! Note: Assuming nsplitA>=nsplitB
    !       So we can use the same buffers for both
    allocate(vxmany(0:ny-1,0:nx-1,nsplit))
    allocate(vymany(0:ny-1,0:nx-1,nsplit))
    allocate(uxmany(0:ny-1,0:nx-1,nsplitA))
    allocate(uymany(0:ny-1,0:nx-1,nsplitA))
    allocate(uv(0:ny-1,0:nx-1,n_omp))
  endif

#endif

#ifdef CGYRO_GPU_FFT
  call cgyro_info('GPU-aware code triggered.')

  ! Note: Assuming nsplitA>=nsplitB
  !       So we can use the same buffers for both
  if (nl_single_flag > 1) then
    ! Use 32-bit NL
    allocate( fxmany32(0:ny2,0:nx-1,nsplitA) )
    allocate( fymany32(0:ny2,0:nx-1,nsplitA) )
    allocate( gxmany32(0:ny2,0:nx-1,nsplit) )
    allocate( gymany32(0:ny2,0:nx-1,nsplit) )

    allocate( uxmany32(0:ny-1,0:nx-1,nsplitA) )
    allocate( uymany32(0:ny-1,0:nx-1,nsplitA) )
    allocate( vxmany32(0:ny-1,0:nx-1,nsplit) )
    allocate( vymany32(0:ny-1,0:nx-1,nsplit) )
    allocate( uvmany32(0:ny-1,0:nx-1,nsplitA) )

#if defined(OMPGPU)
!$omp target enter data map(alloc:fxmany32,fymany32,gxmany32,gymany32) &
!$omp&                  map(alloc:uxmany32,uymany32,vxmany32,vymany32,uvmany32)
#elif defined(_OPENACC)
!$acc enter data create(fxmany32,fymany32,gxmany32,gymany32) &
!$acc&           create(uxmany32,uymany32,vxmany32,vymany32,uvmany32)
#endif
  else
    ! use "regular", 64-bit FP for NL 
    allocate( fxmany(0:ny2,0:nx-1,nsplitA) )
    allocate( fymany(0:ny2,0:nx-1,nsplitA) )
    allocate( gxmany(0:ny2,0:nx-1,nsplit) )
    allocate( gymany(0:ny2,0:nx-1,nsplit) )

    allocate( uxmany(0:ny-1,0:nx-1,nsplitA) )
    allocate( uymany(0:ny-1,0:nx-1,nsplitA) )
    allocate( vxmany(0:ny-1,0:nx-1,nsplit) )
    allocate( vymany(0:ny-1,0:nx-1,nsplit) )
    allocate( uvmany(0:ny-1,0:nx-1,nsplitA) )

#if defined(OMPGPU)
!$omp target enter data map(alloc:fxmany,fymany,gxmany,gymany) &
!$omp&                  map(alloc:uxmany,uymany,vxmany,vymany,uvmany)
#elif defined(_OPENACC)
!$acc enter data create(fxmany,fymany,gxmany,gymany) &
!$acc&           create(uxmany,uymany,vxmany,vymany,uvmany)
#endif
  endif

  write (msg, "(A,I5,A,I5,A,I5)") "NL using FFT batching of ",nsplit,",",nsplitA," and ",nsplitB
  call cgyro_info(msg)

#endif ! CGYRO_GPU_FFT

  call cgyro_nl_fftw_init

  endif ! (nonlinear_flag == 1)


  call timer_lib_out('nl_init')

end subroutine cgyro_init_manager

