!------------------------------------------------
! gyro_write_timedata_hdf5.f90 
!
! PURPOSE:
!  This is an hdf5 version of gyro_write_timedata.f90
!-----------------------------------------------------

subroutine gyro_write_timedata_hdf5

  use gyro_globals
  use hdf5
  use hdf5_api
  use gyro_vshdf5_mod
  use mpi

  !---------------------------------------------------
  implicit none
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
  !
  ! HDF5 variables
  !
  integer, parameter :: hr4=SELECTED_REAL_KIND(6,37)
  character(60) :: description
  character(64) :: step_name, tempVarName
  character(128) :: dumpfile
  character(20)   :: openmethod
  integer(HID_T) :: dumpGid,dumpFid,gid3D,fid3D
  integer(HID_T) :: dumpTGid,dumpTFid
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label
  logical :: write_threed
  !---------------------------------------------------

  !---------------------------------------------------
  ! Determine if the 3D files need to be written 
  if (n_torangle_3d > 1 ) then
     write_threed = .true.
  else
     write_threed = .false.
  endif

  ! io_control: 0=no I/O, 1=open, 2=append, 3=rewind
  select case (io_control)
  case(0)
     return
  case(1)
     openmethod='overwr'
  case(2)
     openmethod='append'
  case(3)
     openmethod='append'
  end select

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

     ! open the timedata file (incremental)
     dumpfile=TRIM(path)//"out.gyro.timedata.h5" 
     description="GYRO scalar time data file"

     call open_h5file(trim(openmethod),dumpfile,dumpTFid,description,dumpTGid,h5in,h5err)
     if (h5err%errBool) call catch_error(h5err%errorMsg)

     dumpfile=TRIM(path)//"gyro"//TRIM(step_name)//".h5"
     description="GYRO field file"
     call open_newh5file(dumpfile,dumpFid,description,dumpGid,h5in,h5err)

     if (write_threed) then
        dumpfile=TRIM(path)//"gyro3D"//TRIM(step_name)//".h5"
        description="GYRO 3D field file"
        call open_newh5file(dumpfile,fid3d,description,gid3D,h5in,h5err)
     endif

     call hdf5_write_coords
  endif
  !---------------------------------------------------

  if (io_control == 1 .or. io_control == 3) return

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
     call write_distributed_complex_h5("phi",&
          dumpGid,gid3D,&
          n_theta_plot*n_x*n_field,&
          n_theta_plot,n_x,n_field,&
          phi_plot(:,:,1:n_field),&
          write_threed,&
          .false., &
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
          .false., &
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
          .false., &
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
          .false., &
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
  call write_distributed_real_h5("kxkyspec",dumpTGid,&
       n_x,1,1,&
       size(kxkyspec),&
       kxkyspec,&
       h5in,h5err)
  call proc_time(cp2)

  if (i_proc == 0) then
     h5in%units="m^-2?"
     h5in%mesh=" "
     call add_h5(dumpTGid,'k_perp_squared',k_perp_squared,h5in,h5err)
  endif

  WRITE(*,*) "passed k_perp squared"
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
     !=============
     ! END LINEAR 
     !=============

  else

     !================
     ! BEGIN NONLINEAR 
     !================

     call proc_time(cp7)

     h5in%units="diff units"
     call write_distributed_real_h5("diff_n",dumpTGid,&
          n_kinetic,n_field,2,&
          size(diff_n),&
          diff_n,&
          h5in,h5err)

     call write_distributed_real_h5("gbflux_n",dumpTGid,&
          n_kinetic,n_field,4,&
          size(gbflux_n),&
          gbflux_n,&
          h5in,h5err)

     if (lindiff_method >= 4) then
        call write_distributed_real_h5('phi_squared_QL_n',dumpTGid,&
             1,1,1,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n,&
             h5in,h5err)

        call write_distributed_real_h5('g_squared_QL_n',dumpTGid,&
             size(g_squared_QL_n),&
             3,1,1,&
             g_squared_QL_n,&
             h5in,h5err)
     endif

     if (nonlinear_transfer_flag == 1) then
        call write_distributed_real_h5('out.gyro.nl_transfer',dumpTGid,&
             n_x,2,1,&
             size(nl_transfer),&
             nl_transfer,&
             h5in,h5err)
     endif

     call proc_time(cp8)

     if (i_proc == 0 ) then

        call add_h5(dumpTGid,'field_rms',ave_phi,h5in,h5err)
        call add_h5(dumpTGid,'diff',diff,h5in,h5err)
        call add_h5(dumpTGid,'diff_i',diff_i,h5in,h5err)
        call add_h5(dumpTGid,'gbflux',gbflux,h5in,h5err)
        call add_h5(dumpTGid,'gbflux_mom',gbflux_mom,h5in,h5err)
        call add_h5(dumpTGid,'gbflux_i',gbflux_i,h5in,h5err)

        if (trapdiff_flag == 1) then
           call add_h5(dumpTGid,'diff_trapped',diff_trapped,h5in,h5err)
           call add_h5(dumpTGid,'diff_i_trapped',diff_i_trapped,h5in,h5err)
           call add_h5(dumpTGid,'gbflux_trapped',gbflux_trapped,h5in,h5err)
           call add_h5(dumpTGid,'gbflux_i_trapped',gbflux_i_trapped,h5in,h5err)
        endif

        allocate(a2(3,n_x))
        a2(1,:) = phi_fluxave(:) 
        a2(2,:) = a_fluxave(:)
        a2(3,:) = aperp_fluxave(:)
        call add_h5(dumpTGid,'zerobar',a2,h5in,h5err)
        deallocate(a2)

        allocate(a3(n_kinetic,4,n_x))
        do i=1,n_x
           a3(:,1,i) = h0_n(:,i)
           a3(:,2,i) = h0_e(:,i)
           a3(:,3,i) = source_n(:,i)
           a3(:,4,i) = source_e(:,i)
        enddo
        call add_h5(dumpTGid,'source',a3,h5in,h5err)
        deallocate(a3)

        call add_h5(dumpTGid,'moments_zero',moments_zero_plot,h5in,h5err)
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
        call add_h5(dumpTGid,'entropy',entropy,h5in,h5err)
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
     call add_h5(dumpTGid,'data_step',data_step,h5in,h5err)
     call add_h5(dumpTGid,'t_current',t_current,h5in,h5err)

     ! dump in the field and 3d files
     call dump_h5(dumpGid,'data_step',data_step,h5in,h5err)
     call dump_h5(dumpGid,'t_current',t_current,h5in,h5err)
     call dump_h5(dumpGid,'n_proc',n_proc,h5in,h5err)
  endif

  if (i_proc == 0) then
     call close_h5file(dumpTFid,dumpTGid,h5err)
     call close_h5file(dumpFid,dumpGid,h5err)
     if (write_threed) call close_h5file(fid3d,gid3d,h5err)
  endif

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
    integer :: iphi, ix, iy, j, ncoarse,nphi

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
    ! Set up the toroidal grid.  Only used for coarse grid
    !-------------------------------------------------

    allocate(zeta_phi(n_torangle_3d))
    do iphi=1,n_torangle_3d
       zeta_phi(iphi)=REAL(iphi-1)/REAL(n_torangle_3d-1)*2.*pi
    end do

    !----------------------------------------
    ! Dump the coarse meshes
    !---------------------------------------- 

    h5in%units=""
    call dump_h5(dumpGid,'Rgyro',Rc,h5in,h5err)
    call dump_h5(dumpGid,'Zgyro',Zc,h5in,h5err)
    h5in%units="m"
    call dump_h5(dumpGid,'R',Rc*a_meters,h5in,h5err)
    call dump_h5(dumpGid,'Z',Zc*a_meters,h5in,h5err)

    ! Here we do not repeat the points since this is the grid
    ! that will be used for the mode plots on thete E [0,2 pi)
    allocate(bufferMesh(0:ncoarse,n_x,2))
    bufferMesh(:,:,1)= Rc*a_meters
    bufferMesh(:,:,2)= Zc*a_meters
    h5in%units="m"
    h5in%mesh="mesh-structured"
    call dump_h5(dumpGid,'cartMesh',bufferMesh(:,:,:),h5in,h5err)
    h5in%mesh=" "
    deallocate(bufferMesh)

    !----------------------------------------
    ! Dump the coarse mesh(es) in 3D
    !---------------------------------------- 
    if (write_threed) then
       !-------------------------------------------------------
       ! Set up the alpha grid
       ! These are set up in a module so no need to recalculate
       !-------------------------------------------------------

       if (.not. allocated(alpha_phi) ) then 
          nphi=1
          if (n_torangle_3d > 0 ) nphi=n_torangle_3d 
          allocate(alpha_phi(0:ncoarse,n_x,nphi))
          do iphi=1,n_torangle_3d
             alpha_phi(:,:,iphi)=zeta_phi(iphi)+nu_coarse(:,:)
          end do
       endif

       h5in%units="m"
       call dump_h5(gid3d,'R',Rc*a_meters,h5in,h5err)
       call dump_h5(gid3d,'Z',Zc*a_meters,h5in,h5err)
       h5in%units="radians"
       call dump_h5(gid3d,'torAngle',zeta_phi,h5in,h5err)
       call dump_h5(gid3d,'torangle_offset',torangle_offset,h5in,h5err)
       call dump_h5(gid3d,'alpha',alpha_phi,h5in,h5err)

       allocate(buffer(ncoarse+1,n_x,n_torangle_3d,3))
       do iphi=1,n_torangle_3d
          buffer(:,:,iphi,1)= Rc(:,:)*COS(zeta_phi(iphi))
          buffer(:,:,iphi,2)=-Rc(:,:)*SIN(zeta_phi(iphi))
          buffer(:,:,iphi,3)= Zc(:,:)
       enddo

       h5in%units="m"; h5in%mesh="mesh-structured"
       call dump_h5(gid3d,'cartMesh',buffer*a_meters,h5in,h5err)
       deallocate(buffer)
    endif

    !----------------------------------------
    ! Dump the wedge mesh(es)
    !---------------------------------------- 
    deallocate(Rc, Zc)
    deallocate(zeta_phi)

  end subroutine hdf5_write_coords

end subroutine gyro_write_timedata_hdf5

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

subroutine write_distributed_real_h5(varName,rGid,n1,n2,n3,n_fn,fn,h5in,h5err)

  use hdf5_api
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err,&
       electron_method

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: varName
  integer(HID_T), intent(in) :: rGid
  integer, intent(in) :: n_fn,n1,n2,n3
  real, intent(in) :: fn(n_fn)
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: ifld,ikin,imom,i
  integer :: in

  real, dimension(:,:,:,:), allocatable :: buffn
  character(128) :: tempVarName 
  character(128), dimension(:,:,:),allocatable :: vnameArray
  character(3) :: n_name
  character(1) :: ikin_name

  !
  real :: fn_recv(n_fn)
  !------------------------------------------------------


  include 'mpif.h'
  !------------------------------------------------------
  !  set up the names for setting species and "field" name
  ALLOCATE(vnameArray(n3,3,4))
  vnameArray=" "
  do ikin=1,n3
     if(electron_method==2 .and. ikin==n3 ) THEN
        tempVarName=trim(varName)//"_electron"
     elseif(electron_method==3 .or. (electron_method==4.and.ikin==n3)) THEN
        tempVarName=trim(varName)//"_electron"
     else
        write(ikin_name,fmt='(i1.1)') ikin-1
        tempVarName=trim(varName)//"_ion"//ikin_name
     endif
     do i=1,4
        vnameArray(ikin,1,i)=trim(tempVarName)//"_phi"
        vnameArray(ikin,2,i)=trim(tempVarName)//"_Apar"
        vnameArray(ikin,3,i)=trim(tempVarName)//"_Bpar"
     enddo
     do i=1,3
        vnameArray(ikin,i,1)=trim(vnameArray(ikin,i,1))//"_density"
        vnameArray(ikin,i,2)=trim(vnameArray(ikin,i,1))//"_energy"
        vnameArray(ikin,i,3)=trim(vnameArray(ikin,i,1))//"_momentum"
        vnameArray(ikin,i,4)=trim(vnameArray(ikin,i,1))//"_energyExchange"
     enddo
  enddo

  !Rectangular array (n_kinetic,n_field,2=i,n_x,n_time)
  !  allocate(kxkyspec(n_x))
  !  allocate(gbflux_n(n_kinetic,n_field,p_moment))
  !  allocate(diff_n(n_kinetic,n_field,n_moment))
  !  real, dimension(1) :: phi_squared_QL_n
  !  real, dimension(3) :: g_squared_QL_n

  ! n1 = n_kinetic; n2 = n_field, n3=n_moment, n4=n_n
  if (i_proc==0) then
     allocate(buffn(n1,n2,n3,n_n)); buffn=0.
  endif
  !-----------------------------------------


  !------------------------------------------------------
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



     if (i_proc == 0) then
        !         WRITE(*,*) "varName=",varName 
        !         WRITE(*,*) " n_fn=",n_fn," and size of fn_recv=",size(fn_recv)
        !         WRITE(*,*) "shape of buffn =", shape(buffn)        
        !         WRITE(*,*) "n1=",n1," n2=",n2," n3=",n3
        buffn(:,:,:,in)=reshape(fn_recv,(/n1,n2,n3/))
     endif
  enddo ! in


  !-----------------------------------------
  if (i_proc /= 0) return


  if (n3==1) then
     if(n2==1) then
        call add_h5(rGid,trim(varName),buffn(:,1,1,:),h5in,h5err)
     else  
        call add_h5(rGid,trim(varName),buffn(:,:,1,:),h5in,h5err)
     endif
  else
     ! n1 = n_kinetic; n2 = n_field, n3=n_moment, n4=n_n
     do ikin=1,n1
        do ifld=1,n2
           do imom=1,n3
              tempVarName=trim(vnameArray(ikin,ifld,imom))
              call add_h5(rGid,trim(tempVarName),buffn(ikin,ifld,imom,:),h5in,h5err)
           enddo
        enddo
     enddo
  endif


  deallocate(buffn)
  deallocate(vnameArray)




end subroutine write_distributed_real_h5
!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex_h5(vname,rGid,r3Did,&
     n_fn,n1,n2,n3,fn,plot3d,plotwedge,h5in,h5err)

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
       n_torangle_wedge,&
       n_torangle_3d

  !------------------------------------------------------
  !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
  !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
  !------------------------------------------------------
  implicit none
  !
  real :: pi=3.141592653589793
  character(*), intent(in) :: vname
  integer(HID_T), intent(in) :: rGid,r3Did
  integer, intent(in) :: n_fn,n1,n2,n3
  complex, intent(in) :: fn(n_fn)
  logical, intent(in) :: plot3d, plotwedge
  character(128) :: tempVarName , tempVarNameGr
  character(128), dimension(:),allocatable :: vnameArray,  pType
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
  logical :: iswedge

  !------------------------------------------------------
  include 'mpif.h'
  c_i=(0,1)

  omega_exp=w0_s(ir_norm)



  if(.not. plotwedge) then
     iswedge=.false.
     allocate(buffn(0:n1,n2,n3,n_n)); buffn=0.
  else
     iswedge=.true.
     allocate(buffn(0:n1-1,n2,n3,n_n)); buffn=0.
  endif

  !when n3=n_kinetic
  !electron_method =1 => n3=n_ion (gk ions and addiabtic electrons )
  !electron_method =2 => n3=n_spec (gk ions and drift electrons)
  !electron_method =3 => n3=1 (gk electrons and addiabtic ions)
  !electron_method =4 => n3=n_ion (gk electrons and gk ions)

  if (trim(vname) /= "phi") then 
     ALLOCATE(vnameArray(n3))
     vnameArray=" "
     do ikin=1,n3
        if(electron_method==2 .and. ikin==n3 ) THEN
           tempVarName=trim(vname)//"_electron"
        elseif(electron_method==3 .or. (electron_method==4.and.ikin==n3)) THEN
           tempVarName=trim(vname)//"_electron"
        else
           write(ikin_name,fmt='(i1.1)') ikin-1
           tempVarName=trim(vname)//"_ion"//ikin_name
        endif
        vnameArray(ikin)=tempVarName
     enddo
  else
     ALLOCATE(vnameArray(3))
     vnameArray=" "
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
  if (.not.iswedge) then
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
  nphi = 1 
  if (.not. iswedge) then
     if (n_torangle_3d > 0) nphi=n_torangle_3d
     allocate(real_buff(0:n1,n2,n3,nphi))
     allocate(alpha_loc(0:n1,n2))
  else
     if (n_torangle_wedge > 0) nphi=n_torangle_wedge
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
     !Get alpha coordinate on either the coarse or wedge mesh.
     ! Include doppler shift here

     if (.not. iswedge) then
        alpha_loc=alpha_phi(:,:,iphi)+omega_exp*t_current
     else
        alpha_loc=alpha_phi_wedge(:,:,iphi)+omega_exp*t_current
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
  if (.not. iswedge) then
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

end subroutine write_distributed_complex_h5
