!------------------------------------------------
! write_hdf5_master_hdf5.f90 
!
! PURPOSE:
!  Write a bunch of data to hdf5 file
!------------------------------------------------
subroutine write_hdf5_data(datafile,action)
  !------------------------------------------
  !  Data that does not change with time.  
  !  It is equivalent to:
  !    profile_vugyro.out
  !------------------------------------------
  use gyro_globals
  use hdf5_api
  use math_constants
  use GEO_interface
  !------------------------------------------
  implicit none
  include 'mpif.h'
  !
  integer, intent(in) :: action
  character (len=*), intent(in) :: datafile
  !
  integer :: data_loop
  integer :: i_dummy
  real :: dummy
  integer(HID_T) :: fid, rootid
  character(90) :: description, filename
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: n_fine
  integer :: io_mode
  real :: theta
  real :: dr
  real, allocatable :: buffer(:,:,:)
  !------------------------------------------
  ! Do the initialization here.  Might need
  ! better logic here based on action.
  !------------------------------------------
  call vshdf5_fcinit()
  if(i_proc/=0) return
  call vshdf5_inith5vars(h5in, h5err)
  h5in%comm=MPI_COMM_SELF
  h5in%info=MPI_INFO_NULL
  h5in%typeConvert=.true.
  h5in%wrd_type=H5T_NATIVE_REAL
  !h5in%wrd_type=H5T_NATIVE_DOUBLE
  h5in%doTranspose=.true.
  !h5in%vsTime=intime
  h5in%wrVsTime=.false.
  h5in%verbose=.true.

  filename=trim(path)//'gyro_profile.h5'

  !---------------------------------------------------------------------
  ! Write the variables to an hdf5 file
  ! These variables are essentially the write_profile_vugyro.f90 
  !---------------------------------------------------------------------
  description=" "
  call open_newh5file(filename,fid,description,rootid,h5in,h5err)

  h5in%mesh=" "; h5in%units=" "
  call dump_h5(rootid,"n_x",n_x,h5in,h5err)
  call dump_h5(rootid,"n_theta_section",n_theta_section,h5in,h5err)
  call dump_h5(rootid,"n_pass",n_pass,h5in,h5err)
  call dump_h5(rootid,"n_trap",n_trap,h5in,h5err)
  call dump_h5(rootid,"n_energy",n_energy,h5in,h5err)
  call dump_h5(rootid,"n_theta_plot",n_theta_plot,h5in,h5err)
  call dump_h5(rootid,"n0",n0,h5in,h5err)
  call dump_h5(rootid,"n_n",n_n,h5in,h5err)
  call dump_h5(rootid,"d_n",d_n,h5in,h5err)
  call dump_h5(rootid,"nonlinear_flag",nonlinear_flag,h5in,h5err)
  call dump_h5(rootid,"electron_method",electron_method,h5in,h5err)
  call dump_h5(rootid,"n_field",n_field,h5in,h5err)
  call dump_h5(rootid,"n_ion",n_ion,h5in,h5err)
  call dump_h5(rootid,"n_kinetic",n_kinetic,h5in,h5err)
  call dump_h5(rootid,"n_spec",n_spec,h5in,h5err)
  call dump_h5(rootid,"field_r0_flag",field_r0_flag,h5in,h5err)
  call dump_h5(rootid,"field_r0_grid",field_r0_grid,h5in,h5err)
  call dump_h5(rootid,"boundary_method",boundary_method,h5in,h5err)
  call dump_h5(rootid,"r",r,h5in,h5err)
  call dump_h5(rootid,"q",q,h5in,h5err)
  call dump_h5(rootid,"Rmaj_s",rmaj_s,h5in,h5err)
  call dump_h5(rootid,"r_s",r_s,h5in,h5err)
  call dump_h5(rootid,"q_s",q_s,h5in,h5err)
  call dump_h5(rootid,"dlntdr_s",dlntdr_s,h5in,h5err)
  call dump_h5(rootid,"dlnndr_s",dlnndr_s,h5in,h5err)
  call dump_h5(rootid,"tem_s",tem_s,h5in,h5err)
  call dump_h5(rootid,"den_s",den_s,h5in,h5err)
  call dump_h5(rootid,"aspect_s",rmaj_s/r_s,h5in,h5err)
  call dump_h5(rootid,"delta_s",delta_s,h5in,h5err)
  call dump_h5(rootid,"zeta_s",zeta_s,h5in,h5err)
  call dump_h5(rootid,"kappa_s",kappa_s,h5in,h5err)
  call dump_h5(rootid,"drmaj_s",drmaj_s,h5in,h5err)
  call dump_h5(rootid,"shat_s",shat_s,h5in,h5err)
  call dump_h5(rootid,"s_delta_s",s_delta_s,h5in,h5err)
  call dump_h5(rootid,"s_zeta_s",s_zeta_s,h5in,h5err)
  call dump_h5(rootid,"s_kappa_s",s_kappa_s,h5in,h5err)
  call dump_h5(rootid,"zmag_s",zmag_s,h5in,h5err)
  call dump_h5(rootid,"dzmag_s",dzmag_s,h5in,h5err)
  call dump_h5(rootid,"beta_unit_s",beta_unit_s,h5in,h5err)
  call dump_h5(rootid,"gamma_e_s",gamma_e_s,h5in,h5err)
  call dump_h5(rootid,"gamma_p_s",gamma_p_s,h5in,h5err)
  call dump_h5(rootid,"mach_s",mach_s,h5in,h5err)
  call dump_h5(rootid,"b_unit_s",b_unit_s,h5in,h5err)
  call dump_h5(rootid,"dr_eodr",dr_eodr,h5in,h5err)
  call dump_h5(rootid,"z_eff_s",z_eff_s,h5in,h5err)
  call dump_h5(rootid,"nu_s",nu_s,h5in,h5err)
  call dump_h5(rootid,"w0_s",w0_s,h5in,h5err)
  call dump_h5(rootid,"box_multiplier",box_multiplier,h5in,h5err)

  call dump_h5(rootid,"lambda", lambda(ir_norm,:),h5in,h5err)
  call dump_h5(rootid,"energy", energy,h5in,h5err)
  call dump_h5(rootid,"lambda_tp", lambda_tp(ir_norm),h5in,h5err)
  call dump_h5(rootid,"krho_collect", krho_collect(:),h5in,h5err)
  call dump_h5(rootid,"rhos_norm", rhos_norm,h5in,h5err)
  call dump_h5(rootid,"zcharge", z(:),h5in,h5err)
  call dump_h5(rootid,"n_moment", n_moment ,h5in,h5err)

  !---------------------------------------------------------------------
  ! These variables are essentially the write_profile_vugyro.f90 
  !---------------------------------------------------------------------

  n_fine = n_theta_plot*n_theta_mult

  allocate(buffer(14,1:n_x,n_fine))
  do i=1,n_x

     if (flat_profile_flag == 0) then

        ! All profiles are global and radial variation is consistent

        GEO_rmin_in      = r_s(i)
        GEO_rmaj_in      = rmaj_s(i)
        GEO_drmaj_in     = drmaj_s(i)
        GEO_zmag_in      = zmag_s(i)
        GEO_dzmag_in     = dzmag_s(i)
        GEO_q_in         = q_s(i)
        GEO_s_in         = shat_s(i)
        GEO_kappa_in     = kappa_s(i)
        GEO_s_kappa_in   = s_kappa_s(i)
        GEO_delta_in     = delta_s(i)
        GEO_s_delta_in   = s_delta_s(i)
        GEO_zeta_in      = zeta_s(i)
        GEO_s_zeta_in    = s_zeta_s(i)
        GEO_beta_star_in = beta_star_s(i)

     else

        ! Profiles are flat and so some parameters need to be linearly extrapolated.

        dr = r(i)-r(ir_norm)

        GEO_rmin_in  = r(i)
        GEO_rmaj_in  = rmaj_s(ir_norm)+drmaj_s(ir_norm)*dr
        GEO_drmaj_in = drmaj_s(ir_norm)
        GEO_zmag_in  = zmag_s(ir_norm)+dzmag_s(ir_norm)*dr
        GEO_dzmag_in = dzmag_s(ir_norm)
        GEO_q_in     = q(i)
        GEO_s_in     = shat_s(ir_norm)
        GEO_kappa_in = kappa_s(ir_norm)+&
                   kappa_s(ir_norm)*s_kappa_s(ir_norm)/r(ir_norm)*dr
        GEO_s_kappa_in = s_kappa_s(ir_norm)
        GEO_delta_in   = delta_s(ir_norm)+s_delta_s(ir_norm)/r(ir_norm)*dr
        GEO_s_delta_in = s_delta_s(ir_norm)
        GEO_zeta_in    = zeta_s(ir_norm)+s_zeta_s(ir_norm)/r(ir_norm)*dr
        GEO_s_zeta_in  = s_zeta_s(ir_norm)
        GEO_beta_star_in = beta_star_s(ir_norm)

     endif

     GEO_fourier_in(:,:) = a_fourier_geo_s(:,0:n_fourier_geo,i)
     call GEO_do()

     do j=1,n_fine

        theta = -pi+(j-1)*pi_2/n_fine

        ! Test for special case
        if (n_fine == 1) theta = 0.0

        call GEO_interp(theta)

        buffer(1 ,i,j)=GEO_nu
        buffer(2 ,i,j)=GEO_gsin
        buffer(3 ,i,j)=GEO_gcos1
        buffer(4 ,i,j)=GEO_gcos2
        buffer(5 ,i,j)=GEO_usin
        buffer(6 ,i,j)=GEO_ucos
        buffer(7 ,i,j)=GEO_b
        buffer(8 ,i,j)=GEO_g_theta
        buffer(9 ,i,j)=GEO_grad_r
        buffer(10,i,j)=GEO_gq
        buffer(11,i,j)=GEO_captheta

     enddo ! j
       ! Set up nu for plotting and synthetic diagnostic.
       ! Note we need theta=-pi and pi so include 0 index
       do j=0,n_theta_plot
          theta = -pi+real(j)*pi_2/n_theta_plot
          if (n_theta_plot == 1) theta = 0.0 ! Test for special case
          call GEO_interp(theta)
          nu_coarse(j,i)=GEO_nu
       enddo
       do j=1,n_fine
          theta=theta_fine_start+real(j-1)*theta_fine_angle/       &
                      real(n_theta_plot*n_theta_mult-1)
          if (n_fine == 1) theta = 0.0 ! Test for special case
          call GEO_interp(theta)
          nu_fine(j,i)=GEO_nu
       enddo
  enddo ! i


  call dump_h5(rootid,"nu",      buffer(1,:,:),h5in,h5err)
  call dump_h5(rootid,"gsin",    buffer(2,:,:),h5in,h5err)
  call dump_h5(rootid,"gcos1",   buffer(3,:,:),h5in,h5err)
  call dump_h5(rootid,"gcos2",   buffer(4,:,:),h5in,h5err)
  call dump_h5(rootid,"usin",    buffer(5,:,:),h5in,h5err)
  call dump_h5(rootid,"ucos",    buffer(6,:,:),h5in,h5err)
  call dump_h5(rootid,"b",       buffer(7,:,:),h5in,h5err)
  call dump_h5(rootid,"g_theta", buffer(8,:,:),h5in,h5err)
  call dump_h5(rootid,"grad_r",  buffer(9,:,:),h5in,h5err)
  call dump_h5(rootid,"gq",      buffer(10,:,:),h5in,h5err)
  call dump_h5(rootid,"captheta",buffer(11,:,:),h5in,h5err)
  deallocate(buffer)

  call close_h5file(fid,rootid,h5err)

 return
 end subroutine write_hdf5_data

!------------------------------------------------------
! write_hdf5_timedata
! PURPOSE:
!  This is an hdf5 version of gyro_write_master.f90
!-----------------------------------------------------

subroutine write_hdf5_timedata(action)
  use gyro_globals
  use hdf5
  use hdf5_api
  use gyro_vshdf5_mod

  !---------------------------------------------------
  implicit none
  include 'mpif.h'
  !
  integer :: mode
  integer, intent(in) :: action
  integer, parameter :: hr4=SELECTED_REAL_KIND(6,37)
  !
  real :: cp0
  real :: cp1
  real :: cp2
  real :: cp3
  real :: cp4
  real :: cp5
  real :: cp6
  real :: cp7
  real :: cp8
  real :: pi=3.141592653589793
  !
  real, dimension(:), allocatable, save :: zeta_phi
  real, dimension(:,:), allocatable :: a2
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(:,:,:), allocatable :: n_plot, e_plot, v_plot
  character(60) :: description
  character(64) :: step_name, tempVarName
  character(128) :: dumpfile
  integer(HID_T) :: dumpGid,dumpFid,gid3D,fid3D
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label

  logical :: write_threed


  !---------------------------------------------------
  ! Determine if the fine meshed files need to be written 
  if (n_alpha_plot > 1 ) then
          write_threed = .true.
  else
          write_threed = .false.
  endif

  !---------------------------------------------------
  ! Initialization
  !---------------------------------------------------
  if (i_proc == 0) then
    call vshdf5_inith5vars(h5in, h5err)
    h5in%comm=MPI_COMM_SELF
    h5in%info=MPI_INFO_NULL
    h5in%wrd_type=H5T_NATIVE_REAL
    h5in%typeConvert=.true.
    !h5in%wrd_type=H5T_NATIVE_DOUBLE
    h5in%doTranspose=.true.
    h5in%verbose=.true.
    h5in%debug=.false.
    h5in%wrVsTime=.true.
    h5in%vsTime=t_current
    h5in%vsStep=step
  
    !---------------------------------------------------
    ! Timestep data:
    !
      number_label=NINT(t_current/dt)
      if (number_label>999999) THEN
        write(step_name,fmt='(i7.7)') number_label
      else if (number_label>99999) THEN
        write(step_name,fmt='(i6.6)') number_label
      else
        write(step_name,fmt='(i5.5)') number_label
      endif
  
    dumpfile=TRIM(path)//"gyro"//TRIM(step_name)//".h5"
    description="GYRO dump file"
    call open_newh5file(dumpfile,dumpFid,description,dumpGid,h5in,h5err)

    if (write_threed) then
      dumpfile=TRIM(path)//"gyro3D"//TRIM(step_name)//".h5"
      description="GYRO 3D plot file"
      call open_newh5file(dumpfile,fid3d,description,gid3D,h5in,h5err)
    endif

    call hdf5_write_coords
  endif
  !---------------------------------------------------

  call proc_time(cp0)

   !--------------------------------------------------
  ! Output of field-like quantities:
  !
  !
  if (plot_u_flag == 1) then
     ! POTENTIALS
     !SEK: Should find out the units here
     !SEK: Need to figure out mesh
     !SEK: h5in%mesh=" "

     !SriV:  This needs to be a dump of phi PLUS both parallel
     ! and perp vector potential
     h5in%units="phi units"
     h5in%mesh="/cartMesh"
     write(*,*) "writing phi"
     call write_distributed_complex_h5("phi",&
          dumpGid,gid3D,&
          n_theta_plot*n_x*n_field,&
          n_theta_plot,n_x,n_field,&
          phi_plot(:,:,1:n_field),&
          write_threed,&
          h5in,h5err)
  endif

    if (plot_n_flag == 1) then
     ! DENSITY
     h5in%units=" "
     h5in%mesh="/cartMesh"
     call write_distributed_complex_h5("density",&
          dumpGid,gid3D,&
           n_theta_plot*n_x*n_kinetic,&
           n_theta_plot,n_x,n_kinetic,&
           moments_plot(:,:,:,1),&
          write_threed,&
          h5in,h5err)
  endif

  if (plot_e_flag == 1) then
     ! ENERGY
     h5in%units="energy units"
     h5in%mesh="/cartMesh"
     call write_distributed_complex_h5("energy",&
          dumpGid,gid3D,&
           n_theta_plot*n_x*n_kinetic,&
           n_theta_plot,n_x,n_kinetic,&
           moments_plot(:,:,:,2),&
          write_threed,&
          h5in,h5err)
  endif

  if (plot_v_flag == 1) then
     ! PARALLEL VELOCITY
     h5in%units="vpar units"
     call write_distributed_complex_h5("v_par",&
          dumpGid,gid3D,&
           n_theta_plot*n_x*n_kinetic,&
           n_theta_plot,n_x,n_kinetic,&
           moments_plot(:,:,:,3),&
          write_threed,&
          h5in,h5err)
     
  endif

  !--------------------------------------------------

  !--------------------------------------------------
  ! Output of field at r=r0:
  !
!-PRE  if (field_r0_flag == 1) then
!-PRE     h5in%units="?"
!-PRE     call write_distributed_complex_h5("field_r0",&
!-PRE          dumpGid,gid3D,&
!-PRE          size(field_r0_plot),&
!-PRE          size(field_r0_plot,1),&
!-PRE          size(field_r0_plot,2),&
!-PRE          size(field_r0_plot,3),&
!-PRE          field_r0_plot,&
!-PRE          write_threed,&
!-PRE          h5in,h5err)
!-PRE  endif
  !--------------------------------------------------

  call proc_time(cp1)
  !Assume gyro_write_master.f90 has calculated this: call get_field_spectrum
  h5in%units="m^-2?"
  h5in%mesh=' '
  call write_distributed_real_h5("kxkyspec",dumpGid,&
       size(kxkyspec),&
       kxkyspec,&
       h5in,h5err)
  call proc_time(cp2)

  if (i_proc == 0) then
     h5in%units="m^-2?"
     h5in%mesh=" "
     call dump_h5(dumpGid,'k_perp_squared',k_perp_squared,h5in,h5err)
  endif


  !-----------------------------
  ! Set remaining timers to zero
  cp3 = 0.0
  cp4 = 0.0
  cp5 = 0.0
  cp6 = 0.0
  cp7 = 0.0
  cp8 = 0.0
  !-----------------------------

  call proc_time(cp3)
  !SEK Assume write_big.f90 has calculated this
  !SEK call get_field_fluxave
  call proc_time(cp4)

  !-------------------------------------------------------------------
  ! Calculation of fundamental nonlinear fluxes and related 
  ! diffusivities
  !
  call gyro_nonlinear_flux
  call proc_time(cp5)
  call gyro_diffusivity
  call gyro_gbflux
  call proc_time(cp6)
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Output specific to linear/nonlinear operation:
  !
  if (nonlinear_flag == 0) then

     !=============
     ! BEGIN LINEAR 
     !=============

     !SEK Worry about this later.
     
     !=============
     ! END LINEAR 
     !=============

  else

     !================
     ! BEGIN NONLINEAR 
     !================

     call proc_time(cp7)

     h5in%units="diff units"
     call write_distributed_real_h5("diff_n",dumpGid,&
          size(diff_n),&
          diff_n,&
          h5in,h5err)

     call write_distributed_real_h5("gbflux_n",dumpGid,&
          size(gbflux_n),&
          gbflux_n,&
          h5in,h5err)

     if (lindiff_method >= 4) then
        call write_distributed_real_h5('phi_squared_QL_n',dumpGid,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n,&
             h5in,h5err)

        call write_distributed_real_h5('g_squared_QL_n',dumpGid,&
             size(g_squared_QL_n),&
             g_squared_QL_n,&
             h5in,h5err)
     endif

     if (nonlinear_transfer_flag == 1) then
        call write_distributed_real_h5('out.gyro.nl_transfer',dumpGid,&
             size(nl_transfer),&
             nl_transfer,&
             h5in,h5err)
     endif

     call proc_time(cp8)

     if (i_proc == 0 .and. lindiff_method > 1) then

        call dump_h5(dumpGid,'field_rms',ave_phi,h5in,h5err)
        call dump_h5(dumpGid,'diff',diff,h5in,h5err)
        call dump_h5(dumpGid,'diff_i',diff_i,h5in,h5err)
        call dump_h5(dumpGid,'gbflux',gbflux,h5in,h5err)
        call dump_h5(dumpGid,'gbflux_mom',gbflux_mom,h5in,h5err)
        call dump_h5(dumpGid,'gbflux_i',gbflux_i,h5in,h5err)

        if (trapdiff_flag == 1) then
           call dump_h5(dumpGid,'diff_trapped',diff_trapped,h5in,h5err)
           call dump_h5(dumpGid,'diff_i_trapped',diff_i_trapped,h5in,h5err)
           call dump_h5(dumpGid,'gbflux_trapped',gbflux_trapped,h5in,h5err)
           call dump_h5(dumpGid,'gbflux_i_trapped',gbflux_i_trapped,h5in,h5err)
        endif

        allocate(a2(3,n_x))
        a2(1,:) = phi_fluxave(:) 
        a2(2,:) = a_fluxave(:)
        a2(3,:) = aperp_fluxave(:)
        call dump_h5(dumpGid,'zerobar',a2,h5in,h5err)
        deallocate(a2)

        allocate(a3(n_kinetic,4,n_x))
        do i=1,n_x
           a3(:,1,i) = h0_n(:,i)
           a3(:,2,i) = h0_e(:,i)
           a3(:,3,i) = source_n(:,i)
           a3(:,4,i) = source_e(:,i)
        enddo
        call dump_h5(dumpGid,'source',a3,h5in,h5err)
        deallocate(a3)

        call dump_h5(dumpGid,'moments_zero',moments_zero_plot,h5in,h5err)
     endif

     !================
     ! END NONLINEAR 
     !================

  endif
  !-------------------------------------------------------------------
  ! Entropy diagnostics
  !
  if (entropy_flag == 1) then
     call gyro_entropy 
     if (i_proc == 0) then 
        call dump_h5(dumpGid,'entropy',entropy,h5in,h5err)
     endif
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Velocity-space diagnostics
  !
  if (velocity_output_flag == 1) then
     !SEK already assumed calculated
     !SEK call get_nonlinear_flux_velocity
     call write_distributed_real_h5('flux_velocity',dumpGid,&
             size(nonlinear_flux_velocity),&
             nonlinear_flux_velocity,&
             h5in,h5err)
  endif
  !------------------------------------------------------------

  !---------------------------------------------------------
  ! Dump restart parameters
  !
  if (i_proc == 0) then
     h5in%mesh=' '
     call dump_h5(dumpGid,'data_step',data_step,h5in,h5err)
     call dump_h5(dumpGid,'t_current',t_current,h5in,h5err)
     call dump_h5(dumpGid,'n_proc',n_proc,h5in,h5err)
  endif

   if (i_proc == 0) then
         call close_h5file(dumpFid,dumpGid,h5err)
         call close_h5file(fid3d,gid3d,h5err)
   endif

   return

  contains
      subroutine hdf5_write_coords
      use GEO_interface
      !------------------------------------------
      !  Write the coordinates out
      !  We want to have same coordinate system as:
      !    allocate(phi_plot(n_theta_plot,n_x,n_field+eparallel_plot_flag))
      !  This should be generalized to include the other GEO options
      !------------------------------------------
       real, dimension(:,:), allocatable :: Rc,Zc,Rf,Zf
       real, dimension(:,:,:,:), allocatable :: buffer
       real, dimension(:,:,:), allocatable :: bufferMesh
       real :: theta, rmajc, zmagc, kappac, deltac, zetac, r_c, dr,xdc
       integer :: iphi, ix, iy, j, ncoarse

       ncoarse = n_theta_plot
       allocate(Rc(0:ncoarse,n_x), Zc(0:ncoarse,n_x))
    
       !----------------------------------------
       ! Calculate the R,Z coordinates.  See write_geometry_arrays.f90
       ! The theta grid needs to correspond to the the theta_plot
       ! array which sets the interpolation arrays in 
       ! gyro_set_blend_arrays.  The theta_plot array is defined as:
       !  do j=1,n_theta_plot
       !          theta_plot(j) = -pi+(j-1)*pi_2/n_theta_plot
       !  enddo
       ! such that theta E [0,2 pi) in gyro_banana_operators.f90
       ! For the 3D arrays, we want the periodic point repeated for
       ! nice plots; i.e., theta E [0,2 pi], but we plot the raw
       ! mode data on theta E [0, 2 pi).  Can be a bit confusing.
       !---------------------------------------- 

       do ix=1,n_x
         r_c=r(ix)
         rmajc = rmaj_s(ix)
         zmagc = zmag_s(ix)
         kappac = kappa_s(ix)
         deltac = delta_s(ix)
         xdc    = asin(deltac)
         zetac  = zeta_s(ix)
         ! Note:  This needs to be generalized for all geometries
         do j=0,ncoarse
             theta = -pi+REAL(j)*pi*2./REAL(ncoarse)
             if(radial_profile_method==1) then
                Rc(j,ix)=rmajc+r_c*cos(theta)
                Zc(j,ix)=zmagc+r_c*sin(theta)
             else
                Rc(j,ix)=rmajc+r_c*cos(theta+xdc*sin(theta))
                Zc(j,ix)=zmagc+kappac*r_c*sin(theta+zetac*sin(2.*theta))
             endif
         enddo
       enddo
      
       !------------------------------------------------
       ! Set up the phi grid.  Only used for coarse grid
       !-------------------------------------------------

       allocate(zeta_phi(n_alpha_plot))
       do iphi=1,n_alpha_plot
          zeta_phi(iphi)=REAL(iphi-1)/REAL(n_alpha_plot-1)*2.*pi
       end do

       !-------------------------------------------------------
       ! Set up the alpha grid
       ! These are set up in a module so no need to recalculate
       !-------------------------------------------------------

       if (.not. allocated(alpha_phi) ) then 
           allocate(alpha_phi(0:ncoarse,n_x,n_alpha_plot))
           do iphi=1,n_alpha_plot
              alpha_phi(:,:,iphi)=zeta_phi(iphi)+nu_coarse(:,:)
           end do
       endif
       !----------------------------------------
       ! Dump the coarse meshes
       !---------------------------------------- 

       call dump_h5(dumpGid,'R',Rc,h5in,h5err)
       call dump_h5(dumpGid,'Z',Zc,h5in,h5err)
       call dump_h5(dumpGid,'zeta_offset',zeta_offset,h5in,h5err)
       call dump_h5(dumpGid,'alpha',alpha_phi,h5in,h5err)

       ! Here we do not repeat the points since this is the grid
       ! that will be used for the mode plots on thete E [0,2 pi)
       allocate(bufferMesh(0:ncoarse,n_x,2))
       bufferMesh(:,:,1)= Rc
       bufferMesh(:,:,2)= Zc
       h5in%units="m"
       h5in%mesh="mesh-structured"
       call dump_h5(dumpGid,'cartMesh',bufferMesh(:,:,:),h5in,h5err)
       h5in%mesh=" "
       deallocate(bufferMesh)

       !----------------------------------------
       ! Dump the coarse mesh(es) in 3D
       !---------------------------------------- 
       if (write_threed) then
         call dump_h5(gid3d,'R',Rc,h5in,h5err)
         call dump_h5(gid3d,'Z',Zc,h5in,h5err)
         call dump_h5(gid3d,'torAngle',zeta_phi,h5in,h5err)
         call dump_h5(gid3d,'alpha',alpha_phi,h5in,h5err)

         allocate(buffer(ncoarse+1,n_x,n_alpha_plot,3))
         do iphi=1,n_alpha_plot
           buffer(:,:,iphi,1)= Rc(:,:)*COS(zeta_phi(iphi))
           buffer(:,:,iphi,2)=-Rc(:,:)*SIN(zeta_phi(iphi))
           buffer(:,:,iphi,3)= Zc(:,:)
         enddo

         h5in%units="m"; h5in%mesh="mesh-structured"
         call dump_h5(gid3d,'cartMesh',buffer,h5in,h5err)
         deallocate(buffer)
      endif

       !----------------------------------------
       ! Dump the fine mesh(es)
       !---------------------------------------- 
       deallocate(Rc, Zc)
       deallocate(zeta_phi)
      return
      end subroutine hdf5_write_coords
 
end subroutine write_hdf5_timedata

!------------------------------------------------------
! write_hdf5_fine_timedata
! PURPOSE:
!  This is like the above, only it is for just the 
!  fine data
!-----------------------------------------------------

subroutine write_hdf5_fine_timedata(action)
  use gyro_globals
  use hdf5
  use hdf5_api
  use gyro_vshdf5_mod

  !---------------------------------------------------
  implicit none
  include 'mpif.h'
  !
  integer :: mode
  integer, intent(in) :: action
  integer, parameter :: hr4=SELECTED_REAL_KIND(6,37)
  !
  real :: cp0
  real :: cp1
  real :: cp2
  real :: cp3
  real :: cp4
  real :: cp5
  real :: cp6
  real :: cp7
  real :: cp8
  real :: pi=3.141592653589793
  !
  real, dimension(:), allocatable, save :: zeta_phi
  real, dimension(:,:), allocatable :: a2
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(:,:,:), allocatable :: n_plot, e_plot, v_plot
  character(60) :: description
  character(64) :: step_name, tempVarName
  character(128) :: dumpfile
  integer(HID_T) :: fidfine,gidfine
  integer :: n_fine
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label


  !---------------------------------------------------
  ! Grid
  !
  n_fine = n_theta_plot*n_theta_mult

  !---------------------------------------------------
  ! SEK: Should I do this every time?
  !---------------------------------------------------
  if (i_proc == 0) then
    call vshdf5_inith5vars(h5in, h5err)
    h5in%comm=MPI_COMM_SELF
    h5in%info=MPI_INFO_NULL
    h5in%wrd_type=H5T_NATIVE_REAL
    h5in%typeConvert=.true.
    !h5in%wrd_type=H5T_NATIVE_DOUBLE
    h5in%doTranspose=.true.
    h5in%verbose=.true.
    h5in%debug=.false.
    h5in%wrVsTime=.true.
    h5in%vsTime=t_current
    h5in%vsStep=step
  
    !---------------------------------------------------
    ! Timestep data:
    !
    number_label=NINT(t_current/dt)
    if (number_label>999999) THEN
      write(step_name,fmt='(i7.7)') number_label
    else if (number_label>99999) THEN
      write(step_name,fmt='(i6.6)') number_label
    else
      write(step_name,fmt='(i5.5)') number_label
    endif

    dumpfile=TRIM(path)//"gyrofine"//TRIM(step_name)//".h5"
    description="GYRO real space file at single phi plane.  On fine mesh"
    call open_newh5file(dumpfile,fidfine,description,gidfine,h5in,h5err)

    call hdf5_write_fine_coords
  endif
  !---------------------------------------------------

  call proc_time(cp0)

   !--------------------------------------------------
  ! Output of field-like quantities:
  !
  !
  if (plot_n_flag == 1) then
     ! DENSITY
     h5in%units=" "
     h5in%mesh="/cartMesh"
     call write_distributed_complex_h5("density",&
        gidfine,gidfine,&
         n_fine*n_x*n_kinetic,&
         n_fine,n_x,n_kinetic,&
         moments_plot_fine(:,:,:,1),&
        .true.,&
        h5in,h5err)
  endif

  if (plot_e_flag == 1) then
     ! ENERGY
     h5in%units="energy units"
     h5in%mesh="/cartMesh"
     call write_distributed_complex_h5("energy",&
         gidfine,gidfine,&
         n_fine*n_x*n_kinetic,&
         n_fine,n_x,n_kinetic,&
         moments_plot_fine(:,:,:,2),&
        .true.,&
        h5in,h5err)
  endif

  if (plot_v_flag == 1) then
     ! PARALLEL VELOCITY
     h5in%units="vpar units"
     h5in%mesh="/cartMesh"
     call write_distributed_complex_h5("v_par",&
        gidfine,gidfine,&
         n_fine*n_x*n_kinetic,&
         n_fine,n_x,n_kinetic,&
         moments_plot_fine(:,:,:,3),&
        .true.,&
        h5in,h5err)
  endif

  !--------------------------------------------------

   if (i_proc == 0) then
      call close_h5file(fidfine,gidfine,h5err)
   endif

   return

  contains
      subroutine hdf5_write_fine_coords
      use GEO_interface
      !------------------------------------------
      !  Write the coordinates out
      !  We want to have same coordinate system as:
      !    allocate(phi_plot(n_theta_plot,n_x,n_field+eparallel_plot_flag))
      !  This should be generalized to include the other GEO options
      !------------------------------------------
       real, dimension(:,:), allocatable :: Rc,Zc,Rf,Zf
       real, dimension(:,:,:), allocatable :: bufferFineMesh
       real :: theta, rmajc, zmagc, kappac, deltac, zetac, r_c, dr,xdc
       real :: zeta_fine
       integer :: iphi, ix, iy, j, ncoarse, nfine

       ncoarse = n_theta_plot
       nfine = n_theta_plot*n_theta_mult
       allocate(Rf(1:nfine,n_x), Zf(1:nfine,n_x))
    
       !----------------------------------------
       ! Calculate the R,Z coordinates.  See write_geometry_arrays.f90
       !---------------------------------------- 

       do ix=1,n_x
         if (flat_profile_flag == 0) then
            r_c=r_s(ix)
         else
            r_c=r(ix)
         endif
         rmajc = rmaj_s(ix)
         zmagc = zmag_s(ix)
         kappac = kappa_s(ix)
         deltac = delta_s(ix)
         xdc    = asin(deltac)
         zetac  = zeta_s(ix)
         do j=1,nfine
             ! This needs to match up with what's in gyro_set_blend_arrays.f90
             theta=theta_fine_start+real(j-1)*theta_fine_angle/       &
                                   real(n_theta_plot*n_theta_mult-1)
             !theta = -pi+REAL(j)*pi*2./REAL(nfine)
             if(radial_profile_method==1) then
                Rf(j,ix)=rmajc+r_c*cos(theta)
                Zf(j,ix)=zmagc+r_c*sin(theta)
             else
                Rf(j,ix)=rmajc+r_c*cos(theta+xdc*sin(theta))
                Zf(j,ix)=zmagc+kappac*r_c*sin(theta+zetac*sin(2.*theta))
             endif
         enddo
       enddo
      
       !-------------------------------------------------------
       ! Set up the alpha grid
       ! These are saved in a module so no need to recalculate
       !-------------------------------------------------------

       if (.not. allocated(alpha_phi_fine) ) then
           allocate(alpha_phi_fine(nfine,n_x,n_alpha_fine))
           do iphi=1,n_alpha_fine
              !Don't store zeta_fine b/c analysis is on a plane by plane basis
              zeta_fine=REAL(iphi-1)/REAL(n_alpha_fine)*2.*pi
              alpha_phi_fine(:,:,iphi)=zeta_offset+zeta_fine+nu_fine(:,:)
           end do
       endif

       !----------------------------------------
       ! Dump the fine meshes
       !---------------------------------------- 

       call dump_h5(gidfine,'R',Rf,h5in,h5err)
       call dump_h5(gidfine,'Z',Zf,h5in,h5err)
       call dump_h5(gidfine,'zeta_offset',zeta_offset,h5in,h5err)
       call dump_h5(gidfine,'alpha',alpha_phi_fine,h5in,h5err)

       ! For ease of use, have a single data set that has R,Z. 
       allocate(bufferFineMesh(nfine,n_x,2))
       bufferFineMesh(:,:,1) = Rf(:,:)
       bufferFineMesh(:,:,2) = Zf(:,:)
       h5in%units="m"
       h5in%mesh="mesh-structured"
       call dump_h5(gidfine,'cartMesh',bufferFineMesh(:,:,:),h5in,h5err)
       h5in%mesh=""
       deallocate(bufferFineMesh)

       !----------------------------------------
       ! 
       !---------------------------------------- 
       deallocate(Rf, Zf)
      return
      end subroutine hdf5_write_fine_coords
 
end subroutine write_hdf5_fine_timedata

  !------------------------------------------------
  ! write_restart
  ! PURPOSE:
  !  File that can be used 
  !------------------------------------------------
subroutine write_hdf5_restart
  use gyro_globals
  use hdf5
  use hdf5_api
  use gyro_vshdf5_mod

  !---------------------------------------------------
  implicit none
  include 'mpif.h'
  !
  real :: pi=3.141592653589793
  character(60) :: description
  character(64) :: step_name, tempVarName
  character(128) :: dumpfile
  integer(HID_T) :: dumpGid,dumpFid,fid3D
  type(hdf5ErrorType) :: errval
  character(4) :: iname
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label

  !---------------------------------------------------
  if (i_proc == 0) then
    call vshdf5_inith5vars(h5in, h5err)
    h5in%comm=MPI_COMM_SELF
    h5in%info=MPI_INFO_NULL
    h5in%wrd_type=H5T_NATIVE_DOUBLE
    h5in%doTranspose=.true.
    h5in%vsTime=t_current
    h5in%wrVsTime=.true.
    h5in%verbose=.true.
  
    !---------------------------------------------------
    ! Timestep data:
    !
      number_label=NINT(t_current/dt)
      if (number_label>999999) THEN
        write(step_name,fmt='(i7.7)') number_label
      else if (data_step>99999) THEN
        write(step_name,fmt='(i6.6)') number_label
      else
        write(step_name,fmt='(i5.5)') number_label
      endif
  
     dumpfile=TRIM(path)//"gyroRestart"//TRIM(step_name)//".h5"
     description="GYRO restart file"
     call open_newh5file(dumpfile,dumpFid,description,dumpGid,h5in,h5err)

     h5in%mesh=' '
     call write_attribute(dumpGid,"data_step",data_step,errval)
     call write_attribute(dumpGid,"n_proc",n_proc,errval)
     call write_attribute(dumpGid,"i_restart",i_restart,errval)
   endif
 
   if (n_proc-1 == 0) then
      call dump_h5(dumpGid,'h_0_real',REAL(h_0),h5in,h5err)
      call dump_h5(dumpGid,'h_0_imag',AIMAG(h_0),h5in,h5err)
      call close_h5file(dumpFid,dumpGid,h5err)
      return
   endif

   do i_proc_w=1,n_proc-1

     if (i_proc == 0) then
        call MPI_RECV(h_0,&
             size(h_0),&
             MPI_DOUBLE_COMPLEX,&
             i_proc_w,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             recv_status,&
             i_err)

        WRITE(iname,fmt='(i4.4)') i_proc_w
        tempVarName="h_0_real"//iname
        call dump_h5(dumpGid,tempVarName,REAL(h_0),h5in,h5err)
        tempVarName="h_0_imag"//iname
        call dump_h5(dumpGid,tempVarName,AIMAG(h_0),h5in,h5err)

     else if (i_proc == i_proc_w) then
        call MPI_SEND(h,&
             size(h),&
             MPI_DOUBLE_COMPLEX,&
             0,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             i_err)
     endif
  enddo

  if (i_proc == 0) call close_h5file(dumpFid,dumpGid,h5err)
 return
 end subroutine write_hdf5_restart


subroutine myhdf5_close
  use gyro_globals
  use hdf5
  use hdf5_api
  use gyro_vshdf5_mod
  integer ierr
  call h5close_f(ierr)
end subroutine myhdf5_close
!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of distributed real array.
!------------------------------------------------------

subroutine write_distributed_real_h5(varName,rGid,n_fn,fn,h5in,h5err)

  use hdf5_api
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: varName
  integer(HID_T), intent(in) :: rGid
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  real :: fn_recv(n_fn)
  character(128) :: n_varName
  character(3) :: n_name
  !------------------------------------------------------


  include 'mpif.h'

     do in=1,n_n

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   i_send,&
                   in,&
                   GYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   0,&
                   in,&
                   GYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        WRITE(n_name,fmt='(i3.3)') in
        n_varName=trim(varName)//"_"//n_name
        if (i_proc == 0) call dump_h5(rGid,n_varName,fn_recv,h5in,h5err)

     enddo ! in

end subroutine write_distributed_real_h5
!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex_h5(vname,rGid,r3Did,&
                     n_fn,n1,n2,n3,fn,plot3d,h5in,h5err)

  use gyro_vshdf5_mod
  use hdf5_api
  use gyro_globals, only : &
       q, &
       n0,&
       d_n,&
       n_n,&
       n_n_1,&
       n_proc_1,&
       n_theta_plot,&
       t_current,&
       debug_flag,&
       recv_status,&
       data_step,&
       w0_s,&
       ir_norm,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err, &
       electron_method,&
       n_alpha_fine,&
       n_alpha_plot

  !------------------------------------------------------
  !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
  !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
  !------------------------------------------------------
  implicit none
  !
  real :: pi=3.141592653589793
  character*(*), intent(in) :: vname
  integer(HID_T), intent(in) :: rGid,r3Did
  integer, intent(in) :: n_fn,n1,n2,n3
  complex, intent(in) :: fn(n_fn)
  logical, intent(in) :: plot3d
  character(128) :: tempVarName , tempVarNameGr
  character(128), dimension(:),allocatable :: vnameArray
  character(3) :: n_name
  character(1) :: ikin_name
  integer(HID_T) :: grGid
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: data_loop
  integer :: i_group_send, ispcs
  integer :: i_send, iphi, istart,nn,i,ikin,in, ix,nphi
  integer :: iloop
  !
  complex :: fn_recv(n_fn), c_i
  complex, dimension(:,:,:,:), allocatable :: buffn
  real, dimension(:,:,:,:), allocatable:: real_buff
  real, dimension(:,:), allocatable:: alpha_loc
  real :: omega_exp
  logical :: iscoarse

  !------------------------------------------------------
  include 'mpif.h'
  c_i=(0,1)

  omega_exp=w0_s(ir_norm)

  if(n1==n_theta_plot) then
     iscoarse=.true.
     allocate(buffn(0:n1,n2,n3,n_n)); buffn=0.
  else
     iscoarse=.false.
     allocate(buffn(0:n1-1,n2,n3,n_n)); buffn=0.
  endif

!when n3=n_kinetic
!electron_method =1 => n3=n_ion (gk ions and addiabtic electrons )
!electron_method =2 => n3=n_spec (gk ions and drift electrons)
!electron_method =3 => n3=1 (gk electrons and addiabtic ions)
!electron_method =4 => n3=n_ion (gk electrons and gk ions)

  if (trim(vname) /= "phi") then 
    ALLOCATE(vnameArray(n3))
    vnameArray=""
    do ikin=1,n3
      if(electron_method==2 .and. ikin==n3 ) THEN
        tempVarName=trim(vname)//"_drift_electron"
      elseif(electron_method==3 .or. (electron_method==4.and.ikin==n3)) THEN
        tempVarName=trim(vname)//"_gk_electron"
      else
        write(ikin_name,fmt='(i1.1)') ikin-1
        tempVarName=trim(vname)//"_ion"//ikin_name
      endif
      vnameArray(ikin)=tempVarName
    enddo
  else
    ALLOCATE(vnameArray(3))
      vnameArray=""
      vnameArray(1)="phi"
      vnameArray(2)="A_par"
      vnameArray(3)="B_par"
  endif

     do in=1,n_n
        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1
        if (i_group_send /= 0) then
           i_send = i_group_send*n_proc_1
           if (i_proc == 0) then
              call MPI_RECV(fn_recv,&
                   n_fn, MPI_DOUBLE_COMPLEX, i_send, in,&
                   GYRO_COMM_WORLD, recv_status, i_err)
           else if (i_proc == i_send) then
              call MPI_SEND(fn,&
                   n_fn, MPI_DOUBLE_COMPLEX, 0, in,&
                   GYRO_COMM_WORLD, i_err)
           endif
        else
           fn_recv(:) = fn(:)
        endif
        
        if (i_proc == 0) then
           buffn(0:n1-1,:,:,in)=reshape(fn_recv,(/n1,n2,n3/))
        endif
            

     enddo ! in

     !-----------------------------------------
     if (i_proc /= 0) return
     !-----------------------------------------
     !-----------------------------------------
     ! Apply boundary conditions
     !-----------------------------------------
     if (iscoarse) then
       do in=1,n_n
         nn=n0+(in-1)*d_n
         !apply theta BC: z_n(r,,2*pi) = z_n(r,0)exp(I*n*(nu(r,2*pi)-nu(r,0)))
         !with nu(r,2*pi) - nu(r,0) = -2*pi*q by definition
         ! phase[*] = EXP(-2*!PI*C_I*n[i_n]*profile_data.q[*])
         do ix=1,n2
           buffn(n1,ix,:,in)=buffn(0,ix,:,in)*exp(-2.*pi*c_i*nn*q(ix))
         enddo
       enddo ! in
     endif

     !-----------------------------------------
     ! Dump each species independently
     !-----------------------------------------
     do ispcs=1,n3
       tempVarNameGr=trim(vnameArray(ispcs))//"_modes"
       call make_group(rGid,trim(tempVarNameGr),grGid,"",h5err)
       tempVarName=trim(vnameArray(ispcs))//"_real"
       call dump_h5(grGid,trim(tempVarName),real(buffn(:,:,ispcs,:)),h5in,h5err)
       tempVarName=trim(vnameArray(ispcs))//"_imag"
       call dump_h5(grGid,trim(tempVarName),aimag(buffn(:,:,ispcs,:)),h5in,h5err)
       CALL close_group(trim(tempVarNameGr),grGid,h5err)
     enddo ! in
     if(.not.plot3d) then
       deallocate(buffn)
       return
     endif

     !-----------------------------------------
     ! Tranform into real space
     !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
     !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
     !-----------------------------------------
     if (iscoarse) then
       nphi=n_alpha_plot
       allocate(real_buff(0:n1,n2,n3,nphi))
       allocate(alpha_loc(0:n1,n2))
     else
       nphi=n_alpha_fine
       allocate(real_buff(0:n1-1,n2,n3,nphi))
       allocate(alpha_loc(0:n1-1,n2))
     endif
     if (n0==0) then
         istart=2
         do iphi=1,nphi
           real_buff(:,:,:,iphi)=real(buffn(:,:,:,1))
         enddo
     else
         istart=1
         real_buff(:,:,:,:)=0.
     endif
     do iphi=1,nphi
       !Get alpha coordinate on either the coarse or fine mesh.
       ! Include doppler shift here
       if (iscoarse) then
         alpha_loc=alpha_phi(:,:,iphi)+omega_exp*t_current
       else
         alpha_loc=alpha_phi_fine(:,:,iphi)+omega_exp*t_current
       endif
       do in=istart,n_n
          nn=n0+(in-1)*d_n
          do ikin=1,n3
            real_buff(:,:,ikin,iphi)=real_buff(:,:,ikin,iphi)&
                  +2.*real(buffn(:,:,ikin,in)*exp(-c_i*nn*alpha_loc(:,:)))
          enddo
       enddo
     enddo
        
     deallocate(buffn,alpha_loc)

     ! Mapping of the variable names to array indices depends on input types
     if (iscoarse) then
       do ikin=1,n3
         tempVarName=trim(vnameArray(ikin))
         call dump_h5(r3Did,trim(tempVarName),real_buff(:,:,ikin,:),h5in,h5err)
       enddo
     else
       ! Dump each phi slice as a separate variable
       do ikin=1,n3
         tempVarNameGr=trim(vnameArray(ikin))//"_toroidal"
         call make_group(r3Did,trim(tempVarNameGr),grGid,"",h5err)
         call dump_h5(grGid,trim(vnameArray(ikin)),real_buff(:,:,ikin,:),h5in,h5err)
         call close_group(trim(tempVarNameGr),grGid,h5err)
       enddo
     endif

     deallocate(real_buff)
     deallocate(vnameArray)

return
end subroutine write_distributed_complex_h5
