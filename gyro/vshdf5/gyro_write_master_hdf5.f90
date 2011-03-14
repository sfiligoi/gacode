!------------------------------------------------
! write_hdf5.f90 [caller write_big]
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
  h5in%doTranspose=.false.
  !h5in%vsTime=intime
  h5in%wrVsTime=.false.
  h5in%verbose=.true.

  filename=trim(path)//'gyro_profile.h5'

  !---------------------------------------------------------------------
  ! Write the variables to an hdf5 file
  !---------------------------------------------------------------------
  description=" "
  !write(*,*) "write_hdf5: opening ", filename, " with h5in%comm = ", h5in%comm
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
  call dump_h5(rootid,"phi_dop_s",phi_dop_s,h5in,h5err)
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
  if (allocated(pgamma_s)) call dump_h5(rootid,"pgamma_s",pgamma_s,h5in,h5err)
  if (allocated(b_unit_s)) call dump_h5(rootid,"b_unit_s",b_unit_s,h5in,h5err)
  call dump_h5(rootid,"dr_eodr",dr_eodr,h5in,h5err)
  call dump_h5(rootid,"z_eff_s",z_eff_s,h5in,h5err)
  call dump_h5(rootid,"nu_s",nu_s,h5in,h5err)
  call dump_h5(rootid,"gamma_eb_s",gamma_eb_s,h5in,h5err)
  call dump_h5(rootid,"w0_s",w0_s,h5in,h5err)
  call dump_h5(rootid,"box_multiplier",box_multiplier,h5in,h5err)

  !SEK: I was trying to get these from the website, but I'm confused
  !SEK: Not sure at all about the names and how they correlate to calculated quantities

  ! chi_i_exp in chi_gb_norm units (main ions)
  if (diff_to_flow(2,1,ir_norm) > 0.0) then
     call dump_h5(rootid,"chi_i_exp",pow_i_s(:)/diff_to_flow(2,1,:),h5in,h5err)
  else 
     call dump_h5(rootid,"chi_i_exp",0.*r(:),h5in,h5err)
  endif

  ! chi_e_exp in chi_gb_norm units (electrons)
  if (diff_to_flow(2,n_spec,ir_norm) > 0.0) then
     call dump_h5(rootid,"chi_e_exp",pow_e_s(:)/diff_to_flow(2,n_spec,:),h5in,h5err)
  else 
     call dump_h5(rootid,"chi_e_exp",0.*r(:),h5in,h5err)
  endif

  call dump_h5(rootid,"diff_to_flow_e1",diff_to_flow(2,1,:),h5in,h5err)
  call dump_h5(rootid,"diff_to_flow_e2",diff_to_flow(2,n_spec,:),h5in,h5err)

  ! Add toroidal viscosity and diff_to_flow for momentum flow
  ! eta_i_exp in chi_gb_norm units plus diff_to_flow
  ! includes convective toroidal velocity flow if mach not zero

  if (abs(diff_to_flow(3,1,ir_norm)) > 0.0) then
    call dump_h5(rootid,"eta_i_tot_exp",flow_mom_s(:)/diff_to_flow(3,1,:),h5in,h5err)
  else
    call dump_h5(rootid,"eta_i_tot_exp",0.*r(:),h5in,h5err)
  endif
  call dump_h5(rootid,"diff_to_flow_mi",0.*r(:),h5in,h5err)
!?SEK  if (allocated(pgamma_s)) then
!?SEK    if (abs(mach_s(1,1)) > 0.0) then 
!?SEK      call dump_h5(rootid,"pgamma_s/mach_s",pgamma_s(1,:)/mach_s(1,:),h5in,h5err)
!?SEK    else
!?SEK      call dump_h5(rootid,"pgamma_s/mach_s",0.*r(:),h5in,h5err)
!?SEK    endif
!?SEK  endif

  ! diff_e_exp in chi_gb_norm units  plus diff_to_flow

  if (diff_to_flow(1,n_spec,ir_norm) > 0.0) then
    call dump_h5(rootid,"diff_ne_exp",powp_s(:)/diff_to_flow(1,n_spec,:),h5in,h5err)
  else 
    call dump_h5(rootid,"diff_ne_exp",0.*r(:),h5in,h5err)
  endif
  call dump_h5(rootid,"diff_to_flow_ne",diff_to_flow(1,n_spec,:),h5in,h5err)

  ! diff_to_flows_heating: 
  !   diff_heating in chi_gb_norm to heating flow in MW
  call dump_h5(rootid,"aolne_exp", &
     diff_to_flow(2,1,:)/(tem_s(1,:)*den_s(1,:)*dlntdr_s(1,:)) ,h5in,h5err)

  call dump_h5(rootid,"lambda", lambda(ir_norm,:),h5in,h5err)
  call dump_h5(rootid,"energy", energy,h5in,h5err)
  call dump_h5(rootid,"lambda_tp", lambda_tp(ir_norm),h5in,h5err)
  call dump_h5(rootid,"krho_collect", krho_collect(:),h5in,h5err)
  call dump_h5(rootid,"rhos_norm", rhos_norm,h5in,h5err)
  call dump_h5(rootid,"zcharge", z(:),h5in,h5err)
  call dump_h5(rootid,"n_moment", n_moment ,h5in,h5err)

  call close_h5file(fid,rootid,h5err)
 return
 end subroutine write_hdf5_data

!------------------------------------------------------
! write_hdf5_timedata
! PURPOSE:
!  This is an hdf5 version of write_big.f90
!-----------------------------------------------------

subroutine write_hdf5_timedata(action)
  use gyro_globals
  use math_constants
  use hdf5
  use hdf5_api
  use hdf5_mod

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
  !
  real, dimension(:), allocatable, save :: zeta_phi
  real, dimension(:,:), allocatable :: a2
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(:,:,:), allocatable :: n_plot, e_plot, v_plot
  character(60) :: description
  character(64) :: step_name, tempVarName
  character(128) :: dumpfile
  integer(HID_T) :: dumpGid,dumpFid,gid3D,fid3D,gridGid,fidfine,gidfine,grdfine,grdcoarse
  integer :: n_fine
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err

  logical :: write_fine


  !---------------------------------------------------
  ! Determine if the fine meshed files need to be written 
  if (n_theta_mult == 1 ) then
          write_fine = .false.
  else
          write_fine = .true.
  endif

  !---------------------------------------------------
  ! Determine file control mode:
  ! mode = 1 -> file create 
  !      = 2 -> file write 
  !      = 3 -> file reposition 
  !      = 4 -> file write on step = 0
  ! action = 1 -> simulation still initializing
  !        = 2 -> simulation running

  if (action == 1) then
     ! File creation or repositioning:
     select case (restart_method)
     case(-1) 
        mode = 1 ! No use of restart facility
        !SEK: This is confusing.  At the beginning, not everything 
        !        ! is allocated and setup, so just return
        return
     case(0,2) 
        mode = 1 ! Start of restartable simulation
        !SEK: This is confusing.  At the beginning, not everything 
        !        ! is allocated and setup, so just return
        return
     case(1)
        mode = 3 ! Continuation of restartable simulation        
     end select
  else
     ! File writing
     if (step == 0) then
        mode = 4
     else
        mode = 2
     endif
  endif
  if (output_flag == 0) mode = -mode

  !---------------------------------------------------
  ! Grid
  !
  n_fine = n_theta_plot*n_theta_mult

  !---------------------------------------------------
  if (i_proc == 0) then
    call vshdf5_inith5vars(h5in, h5err)
    h5in%comm=MPI_COMM_SELF
    h5in%info=MPI_INFO_NULL
    h5in%wrd_type=H5T_NATIVE_REAL
    h5in%typeConvert=.true.
    !h5in%wrd_type=H5T_NATIVE_DOUBLE
    h5in%doTranspose=.false.
    h5in%verbose=.true.
    h5in%debug=.false.
    h5in%wrVsTime=.true.
    h5in%vsTime=t_current
    h5in%vsStep=step
  
    !---------------------------------------------------
    ! Timestep data:
    !
      if (step>999999) THEN
        write(step_name,fmt='(i7.7)') step
      else if (step>99999) THEN
        write(step_name,fmt='(i6.6)') step
      else
        write(step_name,fmt='(i5.5)') step
      endif
  
    dumpfile=TRIM(path)//"gyro"//TRIM(step_name)//".h5"
    description="GYRO dump file"
    call open_newh5file(dumpfile,dumpFid,description,dumpGid,h5in,h5err)

    dumpfile=TRIM(path)//"gyro3D"//TRIM(step_name)//".h5"
    description="GYRO 3D plot file"
    call open_newh5file(dumpfile,fid3d,description,gid3D,h5in,h5err)

    if(write_fine) then
      dumpfile=TRIM(path)//"gyrofine"//TRIM(step_name)//".h5"
      description="GYRO real space file at single phi plane.  On fine mesh"
      call open_newh5file(dumpfile,fidfine,description,gidfine,h5in,h5err)
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
     h5in%mesh="/grid/cartMesh"
     write(*,*) "writing phi"
     call write_distributed_complex_h5("phi",&
          dumpGid,gid3D,&
          n_theta_plot*n_x*n_field,&
          n_theta_plot,n_x,n_field,&
          phi_plot(:,:,1:n_field),&
          .true.,&
          h5in,h5err)
  endif

    if (plot_n_flag == 1) then
     ! DENSITY
     h5in%units=" "
     h5in%mesh="/grid/cartMesh"
     call write_distributed_complex_h5("density",&
          dumpGid,gid3D,&
           n_theta_plot*n_x*n_kinetic,&
           n_theta_plot,n_x,n_kinetic,&
           moments_plot(:,:,:,1),&
          .true.,&
          h5in,h5err)
     if(write_fine) then
       call write_distributed_complex_h5("density",&
          gidfine,gidfine,&
           n_theta_mult*n_theta_plot*n_x*n_kinetic,&
           n_theta_mult*n_theta_plot,n_x,n_kinetic,&
           moments_plot_fine(:,:,:,1),&
          .true.,&
          h5in,h5err)
     endif
  endif

  if (plot_e_flag == 1) then
     ! ENERGY
     h5in%units="energy units"
     h5in%mesh="/grid/cartMesh"
     call write_distributed_complex_h5("energy",&
          dumpGid,gid3D,&
           n_theta_plot*n_x*n_kinetic,&
           n_theta_plot,n_x,n_kinetic,&
           moments_plot(:,:,:,2),&
          .true.,&
          h5in,h5err)

     if(write_fine) then
       call write_distributed_complex_h5("energy",&
           gidfine,gidfine,&
           n_theta_mult*n_theta_plot*n_x*n_kinetic,&
           n_theta_mult*n_theta_plot,n_x,n_kinetic,&
           moments_plot_fine(:,:,:,2),&
          .true.,&
          h5in,h5err)
     endif
  endif

  if (plot_v_flag == 1) then
     ! PARALLEL VELOCITY
     h5in%units="vpar units"
     call write_distributed_complex_h5("v_par",&
          dumpGid,gid3D,&
           n_theta_plot*n_x*n_kinetic,&
           n_theta_plot,n_x,n_kinetic,&
           moments_plot(:,:,:,3),&
          .true.,&
          h5in,h5err)
     
     if(write_fine) then
       call write_distributed_complex_h5("v_par",&
          gidfine,gidfine,&
           n_theta_mult*n_theta_plot*n_x*n_kinetic,&
           n_theta_mult*n_theta_plot,n_x,n_kinetic,&
           moments_plot_fine(:,:,:,3),&
          .true.,&
          h5in,h5err)
     endif
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
!-PRE          .true.,&
!-PRE          h5in,h5err)
!-PRE  endif
  !--------------------------------------------------

  call proc_time(cp1)
  !SEK Assume write_big.f90 has calculated this
  !SEK call get_field_spectrum
  h5in%units="m^-2?"
  call write_distributed_real_h5("kxkyspec",dumpGid,&
       size(kxkyspec),&
       kxkyspec,&
       h5in,h5err)
  call proc_time(cp2)

  if (i_proc == 0) then
     h5in%units="m^-2?"
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
  if (rotation_method == 2) then 
     call get_nonlinear_flux
     call proc_time(cp5)
     call get_diffusivity
     call proc_time(cp6)
  else
     call gyro_nonlinear_flux
     call proc_time(cp5)
     call gyro_diffusivity
     if (rotation_method == 3) call gyro_gbflux
     call proc_time(cp6)
  endif
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

     if (rotation_method == 3) then
        call write_distributed_real_h5("gbflux_n",dumpGid,&
             size(gbflux_n),&
             gbflux_n,&
             h5in,h5err)
     endif

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
        call write_distributed_real_h5('nonlinear_transfer_n',dumpGid,&
             size(Tr_p),&
             Tr_p,&
             h5in,h5err)
        call write_distributed_real_h5('turbulent_energy_n',dumpGid,&
             size(Eng_p),&
             Eng_p,&
             h5in,h5err)
     endif

     call proc_time(cp8)

     if (i_proc == 0) then

        call dump_h5(dumpGid,'field_rms',ave_phi,h5in,h5err)
        call dump_h5(dumpGid,'diff',diff,h5in,h5err)
        call dump_h5(dumpGid,'diff_i',diff_i,h5in,h5err)

        if (rotation_method == 3) then
           call dump_h5(dumpGid,'gbflux',gbflux,h5in,h5err)
           call dump_h5(dumpGid,'gbflux_i',gbflux_i,h5in,h5err)
        endif

        if (rotation_method == 2) then
           call dump_h5(dumpGid,'diff_i_ch',diff_i_ch,h5in,h5err)
           call dump_h5(dumpGid,'sp_diff',sp_diff,h5in,h5err)
           call dump_h5(dumpGid,'s_diff',s_diff,h5in,h5err)
           call dump_h5(dumpGid,'sp_diff_i',sp_diff_i,h5in,h5err)
           call dump_h5(dumpGid,'s_diff_i',s_diff_i,h5in,h5err)
        endif

        if (trapdiff_flag == 1) then
           call dump_h5(dumpGid,'diff_trapped',diff_trapped,h5in,h5err)
           call dump_h5(dumpGid,'diff_i_trapped',diff_i_trapped,h5in,h5err)
           if (rotation_method == 3) then
              call dump_h5(dumpGid,'gbflux_trapped',gbflux_trapped,h5in,h5err)
              call dump_h5(dumpGid,'gbflux_i_trapped',gbflux_i_trapped,h5in,h5err)
           endif
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
     call dump_h5(dumpGid,'data_step',data_step,h5in,h5err)
     call dump_h5(dumpGid,'t_current',t_current,h5in,h5err)
     call dump_h5(dumpGid,'n_proc',n_proc,h5in,h5err)
  endif

   if (i_proc == 0) then
         call close_h5file(dumpFid,dumpGid,h5err)
         call close_h5file(fid3d,gid3d,h5err)
         if(write_fine) then
           call close_h5file(fidfine,gidfine,h5err)
         endif
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
       real :: theta, rmajc, zmagc, kappac, deltac, zetac, r_c, dr,xdc
       real :: zeta_fine
       integer :: iphi, ix, iy, j, ncoarse, nfine

       ncoarse = n_theta_plot
       nfine = n_theta_plot*n_theta_mult
       allocate(Rc(0:ncoarse,n_x), Zc(0:ncoarse,n_x))
       allocate(Rf(0:nfine,n_x), Zf(0:nfine,n_x))
    
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
!SEK: I am totally confused here.  This is in write_geometry arrays, but isn't used
!SEK: by Chris
!            dr = r(ix)-r(ir_norm)
!            rmajc = rmaj_s(ir_norm)+drmaj_s(ir_norm)*dr
!            zmagc = zmag_s(ir_norm)+dzmag_s(ir_norm)*dr
!            kappac = kappa_s(ir_norm)+kappa_s(ir_norm)*s_kappa_s(ir_norm)/r(ir_norm)*dr
!            deltac = delta_s(ir_norm)+s_delta_s(ir_norm)/r(ir_norm)*dr
!            zetac  = zeta_s(ir_norm) +s_zeta_s(ir_norm)/r(ir_norm)*dr
         do j=0,ncoarse
             theta = -pi+REAL(j)*pi_2/REAL(ncoarse)
             if(radial_profile_method==1) then
                Rc(j,ix)=rmajc+r_c*cos(theta)
                Zc(j,ix)=zmagc+r_c*sin(theta)
             else
                Rc(j,ix)=rmajc+r_c*cos(theta+xdc*sin(theta))
                Zc(j,ix)=zmagc+kappac*r_c*sin(theta+zetac*sin(2.*theta))
             endif
         enddo
         do j=0,nfine
             theta = -pi+REAL(j)*pi_2/REAL(nfine)
             if(radial_profile_method==1) then
                Rf(j,ix)=rmajc+r_c*cos(theta)
                Zf(j,ix)=zmagc+r_c*sin(theta)
             else
                Rf(j,ix)=rmajc+r_c*cos(theta+xdc*sin(theta))
                Zf(j,ix)=zmagc+kappac*r_c*sin(theta+zetac*sin(2.*theta))
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
       if (.not. allocated(alpha_phi_fine) ) then
           allocate(alpha_phi_fine(0:nfine,n_x,n_alpha_fine))
           do iphi=1,n_alpha_fine
              !Don't store zeta_fine b/c analysis is on a plane by plane basis
              zeta_fine=REAL(iphi-1)/REAL(n_alpha_fine)*2.*pi
              alpha_phi_fine(:,:,iphi)=zeta_offset+zeta_fine+nu_fine(:,:)
           end do
       endif
       !----------------------------------------
       ! Dump the course mesh(es)
       !---------------------------------------- 

         call make_group(dumpFid,"grid", grdcoarse,"RZ grid for course mesh",h5err)
         call dump_h5(grdcoarse,'R',Rc,h5in,h5err)
         call dump_h5(grdcoarse,'Z',Zc,h5in,h5err)
         call dump_h5(grdcoarse,'zeta_offset',zeta_offset,h5in,h5err)
         call dump_h5(grdcoarse,'alpha',alpha_phi,h5in,h5err)

       ! For ease of use, have a single data set that has R,Z. 
         allocate(buffer(2,0:ncoarse,n_x,1))
         buffer(1,:,:,1)= Rc(:,:)
         buffer(2,:,:,1)= Zc(:,:)
         h5in%units="m"
         h5in%mesh="mesh-structured"
         call dump_h5(grdcoarse,'cartMesh',buffer(:,:,:,1),h5in,h5err)
         h5in%mesh=""
         deallocate(buffer)
         call close_group("grid",grdcoarse,h5err)




       !----------------------------------------
       ! Dump the coarse mesh(es) in 3D
       !---------------------------------------- 
!SEK - still defining
!       call make_mesh_group(fid3d, gridGid,h5in,"grid", "structured",&
!              "R","Z","phi","cylindrical","cylGrid",h5err)
!       call make_mesh_group(fid3d, gridGid,h5in,"grid", "structured",&
!              "R","Z","phi"," ","cylGrid",h5err)

       call make_group(fid3d,"grid", gridGid,"All of the various grids",h5err)
       call dump_h5(gridGid,'R',Rc,h5in,h5err)
       call dump_h5(gridGid,'Z',Zc,h5in,h5err)
       call dump_h5(gridGid,'phi',zeta_phi,h5in,h5err)
       call dump_h5(gridGid,'alpha',alpha_phi,h5in,h5err)

       allocate(buffer(3,0:ncoarse,n_x,n_alpha_plot))
       do iphi=1,n_alpha_plot
         buffer(1,:,:,iphi)= Rc(:,:)*COS(zeta_phi(iphi))
         buffer(2,:,:,iphi)=-Rc(:,:)*SIN(zeta_phi(iphi))
         buffer(3,:,:,iphi)= Zc(:,:)
       enddo
       h5in%units="m"; h5in%mesh="mesh-structured"
       call dump_h5(gridGid,'cartMesh',buffer,h5in,h5err)
       call close_group("grid",gridGid,h5err)
       deallocate(buffer)

       !----------------------------------------
       ! Dump the fine mesh(es)
       !---------------------------------------- 

       if(write_fine) then
         call make_group(fidfine,"grid", grdfine,"RZ grid for fine mesh",h5err)
         call dump_h5(grdfine,'R',Rf,h5in,h5err)
         call dump_h5(grdfine,'Z',Zf,h5in,h5err)
         call dump_h5(grdfine,'zeta_offset',zeta_offset,h5in,h5err)
         call dump_h5(grdfine,'alpha',alpha_phi_fine,h5in,h5err)

       ! For ease of use, have a single data set that has R,Z. 
         allocate(buffer(2,0:nfine,n_x,1))
         buffer(1,:,:,1)= Rf(:,:)
         buffer(2,:,:,1)= Zf(:,:)
         h5in%units="m"
         h5in%mesh="mesh-structured"
         call dump_h5(grdfine,'cartMesh',buffer(:,:,:,1),h5in,h5err)
         h5in%mesh=""
         deallocate(buffer)
         call close_group("grid",grdfine,h5err)
       endif

       !----------------------------------------
       ! 
       !---------------------------------------- 
       deallocate(Rc, Zc)
       deallocate(zeta_phi)
      return
      end subroutine hdf5_write_coords
 
end subroutine write_hdf5_timedata

  !------------------------------------------------
  ! write_restart
  ! PURPOSE:
  !  File that can be used 
  !------------------------------------------------
subroutine write_hdf5_restart
  use gyro_globals
  use math_constants
  use hdf5
  use hdf5_api
  use hdf5_mod

  !---------------------------------------------------
  implicit none
  include 'mpif.h'
  !
  character(60) :: description
  character(64) :: step_name, tempVarName
  character(128) :: dumpfile
  integer(HID_T) :: dumpGid,dumpFid,fid3D,gridGid
  type(hdf5ErrorType) :: errval
  character(4) :: iname
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err

  !---------------------------------------------------
  if (i_proc == 0) then
    call vshdf5_inith5vars(h5in, h5err)
    h5in%comm=MPI_COMM_SELF
    h5in%info=MPI_INFO_NULL
    h5in%wrd_type=H5T_NATIVE_DOUBLE
    h5in%doTranspose=.false.
    h5in%vsTime=t_current
    h5in%wrVsTime=.true.
    h5in%verbose=.true.
  
    !---------------------------------------------------
    ! Timestep data:
    !
      if (step>999999) THEN
        write(step_name,fmt='(i7.7)') step
      else if (data_step>99999) THEN
        write(step_name,fmt='(i6.6)') step
      else
        write(step_name,fmt='(i5.5)') step
      endif
  
     dumpfile=TRIM(path)//"gyroRestart"//TRIM(step_name)//".h5"
     description="GYRO restart file"
     call open_newh5file(dumpfile,dumpFid,description,dumpGid,h5in,h5err)

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
  use math_constants
  use hdf5
  use hdf5_api
  use hdf5_mod
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

  use math_constants
  use hdf5_mod
  use hdf5_api
  use gyro_globals, only : &
       q, &
       n0,&
       d_n,&
       n_n,&
       n_n_1,&
       n_proc_1,&
       n_theta_plot,&
       omega_exp,&
       t_current,&
       debug_flag,&
       recv_status,&
       data_step,&
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
  character*(*), intent(in) :: vname
  integer(HID_T), intent(in) :: rGid,r3Did
  integer, intent(in) :: n_fn,n1,n2,n3
  complex, intent(in) :: fn(n_fn)
  logical, intent(in) :: plot3d
  character(128) :: tempVarName
  character(3) :: n_name
  character(1) :: ikin_name
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send, iphi, istart,nn,i,ikin,in, ix,nphi
  !
  complex :: fn_recv(n_fn), c_i
  complex, dimension(:,:,:,:), allocatable :: buffn
  real, dimension(:,:,:,:), allocatable:: real_buff
  real, dimension(:,:), allocatable:: alpha_loc
  logical :: iscoarse

  !------------------------------------------------------
  include 'mpif.h'
  c_i=(0,1)

  if(n1==n_theta_plot) then
     iscoarse=.true.
  else
     iscoarse=.false.
  endif

  allocate(buffn(0:n1,n2,n3,n_n)); buffn=0.

     do in=1,n_n
  !WRITE(*,*) "in ", in, i_proc
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
     ! Dump each mode in the same format that gyro does
     !-----------------------------------------
     if (iscoarse) then
       do in=1,n_n
         WRITE(n_name,fmt='(i3.3)') in
         tempVarName=trim(vname)//"_real"//n_name
         call dump_h5(rGid,trim(tempVarName),real(buffn(0:n1-1,:,:,in)),h5in,h5err)
         tempVarName=trim(vname)//"_imag"//n_name
         call dump_h5(rGid,trim(tempVarName),aimag(buffn(0:n1-1,:,:,in)),h5in,h5err)
       enddo ! in
     else
!      if (debug_flag == 1) then
        do in=1,n_n
         WRITE(n_name,fmt='(i3.3)') in
         tempVarName=trim(vname)//"fine_real"//n_name
         call dump_h5(rGid,trim(tempVarName),real(buffn(0:n1-1,:,:,in)),h5in,h5err)
         tempVarName=trim(vname)//"fine_imag"//n_name
         call dump_h5(rGid,trim(tempVarName),aimag(buffn(0:n1-1,:,:,in)),h5in,h5err)
       enddo ! in
!      endif
     endif
     if(.not.plot3d) then
       deallocate(buffn)
       return
     endif
     !-----------------------------------------
     ! Apply boundary conditions
     !-----------------------------------------
     do in=1,n_n
       nn=n0+(in-1)*d_n
       !apply theta BC: z_n(r,,2*pi) = z_n(r,0)exp(I*n*(nu(r,2*pi)-nu(r,0)))
       !with nu(r,2*pi) - nu(r,0) = -2*pi*q by definition
       ! phase[*] = EXP(-2*!PI*C_I*n[i_n]*profile_data.q[*])
       do ix=1,n2
         buffn(n1,ix,:,in)=buffn(0,ix,:,in)*exp(-2.*pi*c_i*nn*q(ix))
       enddo
     enddo ! in

     !-----------------------------------------
     ! Tranform into real space
     !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
     !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
     !-----------------------------------------
     if (iscoarse) then
       nphi=n_alpha_plot
     else
       nphi=n_alpha_fine
     endif
     allocate(real_buff(0:n1,n2,n3,nphi))
     allocate(alpha_loc(0:n1,n2))
!sv     Default for n0=30, because of k_rho_s scaling.
!       I think we always want to count over all toroidal modes calculated
!       which means below, going from istart=1 to n_n (toroidal grid).
!     if (n0==0) then
         istart=1
         do iphi=1,nphi
           real_buff(:,:,:,iphi)=real(buffn(:,:,:,1))
         enddo
!     else
!         istart=2
!         real_buff(:,:,:,:)=0.
!         !where does real_buff come from in thise case??/
!     endif
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
         if (trim(vname) /= "phi") then 
           ! See gyro_select_methods for understanding this logic
           if(electron_method==2 .and. ikin==n3) THEN
              tempVarName=trim(vname)//"_electron"
           elseif(electron_method==3) THEN
              tempVarName=trim(vname)//"_electron"
           else
              write(ikin_name,fmt='(i1.1)') ikin
              tempVarName=trim(vname)//"_ion"//ikin_name
           endif
         else
           if(ikin==1) tempVarName="phi"
           if(ikin==2) tempVarName="A_par"
           if(ikin==3) tempVarName="B_par"
         endif
         call dump_h5(r3Did,trim(tempVarName),real_buff(:,:,ikin,:),h5in,h5err)
       enddo
     else
       ! Dump each phi slice as a separate variable
       do iphi=1,nphi
         do ikin=1,n3
           if (trim(vname) /= "phi") then 
             if(electron_method==2 .and. ikin==n3) THEN
                tempVarName=trim(vname)//"_electron"
             elseif(electron_method==3) THEN
                tempVarName=trim(vname)//"_electron"
             else
                write(ikin_name,fmt='(i1.1)') ikin-1
                tempVarName=trim(vname)//"_ion"//ikin_name
             endif
           else
             if(ikin==1) tempVarName="phi"
             if(ikin==2) tempVarName="A_par"
             if(ikin==3) tempVarName="B_par"
           endif
           write(n_name,fmt='(i2.2)') iphi
           tempVarName=trim(tempVarName)//"_phi"//TRIM(n_name)
           call dump_h5(r3Did,trim(tempVarName),real_buff(:,:,ikin,iphi),h5in,h5err)
         enddo
       enddo
     endif

     deallocate(real_buff)

return
end subroutine write_distributed_complex_h5
