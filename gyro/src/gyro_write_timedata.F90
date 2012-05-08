!------------------------------------------------------
! gyro_write_timedata.F90
!
! PURPOSE:
!  This is the master file controlling output of
!  data on a per-timestep basis.  This file also 
!  contains the MPI IO routines 
!
!  - write_distributed_real
!  - write_distributed_complex
!  - write_local_real
!-----------------------------------------------------

subroutine gyro_write_timedata

  use gyro_globals
#ifdef HAVE_HDF5
  use hdf5_api
  use gyro_vshdf5_mod
  use mpi
#endif

  !---------------------------------------------------
  implicit none
  !
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: n_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: e_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: v_plot
  !
  real :: pi=3.141592653589793
#ifdef HAVE_HDF5
  integer, parameter :: hr4=SELECTED_REAL_KIND(6,37)
  character(60) :: description
  character(64) :: step_name
  character(128) :: dumpfile
  character(20)   :: openmethod
  integer(HID_T) :: dumpGid,dumpFid,gid3D,fid3D
  integer(HID_T) :: dumpTGid,dumpTFid
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label
  logical :: write_threed
  !integer, INTENT(IN) :: h5_control
#endif

  !---------------------------------------------------

  !---------------------------------------------------
  ! Timestep data:
  !
  if (i_proc == 0) then
     call gyro_write_step(trim(path)//'out.gyro.t',1)
  endif
  !---------------------------------------------------
#ifdef HAVE_HDF5
!---------------------------------------------------
  ! Determine if the 3D files need to be written 
  if (n_torangle_3d > 1 ) then
     write_threed = .true.
  else
     write_threed = .false.
  endif
#endif
  !--------------------------------------------------
  ! Output of field-like quantities:
  !
  if (plot_n_flag+plot_e_flag+plot_v_flag > 0) then
     n_plot(:,:,:) = moments_plot(:,:,:,1)
     e_plot(:,:,:) = moments_plot(:,:,:,2)
     v_plot(:,:,:) = moments_plot(:,:,:,3)
  endif
  !
  if (plot_u_flag == 1) then

     ! POTENTIALS

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_u',&
          10,&
          size(phi_plot(:,:,1:n_field)),&
          phi_plot(:,:,1:n_field))

  endif

  if (plot_epar_flag == 1) then

     ! PARALLEL ELECTRIC FIELD

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_epar',&
          10,&
          size(phi_plot(:,:,n_field+1)),&
          phi_plot(:,:,n_field+1))

  endif

  if (plot_n_flag == 1) then

     ! DENSITY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_n',&
          10,&
          size(n_plot),&
          n_plot)

  endif

  if (plot_e_flag == 1) then

     ! ENERGY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_e',&
          10,&
          size(e_plot),&
          e_plot)

  endif

  if (plot_v_flag == 1) then

     ! PARALLEL VELOCITY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_v',&
          10,&
          size(v_plot),&
          v_plot)

  endif

  !--------------------------------------------------

  !--------------------------------------------------
  ! Output of field at r=r0:
  !
  if (field_r0_flag == 1) then
     call write_distributed_complex(&
          trim(path)//'out.gyro.field_r0',&
          10,&
          size(field_r0_plot),&
          field_r0_plot)
  endif
  !--------------------------------------------------

  call gyro_kxky_spectrum
  call write_distributed_real(&
       trim(path)//'out.gyro.kxkyspec',&
       10,&
       size(kxkyspec),&
       kxkyspec)

  if (i_proc == 0) then
     call write_local_real(&
          trim(path)//'out.gyro.k_perp_squared',&
          10,&
          size(k_perp_squared),&
          k_perp_squared)
  endif

  call gyro_field_fluxave

  !-------------------------------------------------------------------
  ! Calculation of fundamental nonlinear fluxes and related 
  ! diffusivities
  !
  call gyro_nonlinear_flux
  call gyro_diffusivity
  call gyro_gbflux
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Output specific to linear/nonlinear operation:
  !
  if (nonlinear_flag == 0) then

     !=============
     ! BEGIN LINEAR 
     !=============

     call gyro_write_freq(trim(path)//'out.gyro.freq',10)

     if (plot_u_flag == 1) then        

        ! PHI
        call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_phi',10,1,0)

        if (n_field > 1) then
           ! A_PARALLEL 
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_a',10,2,0)
        endif

        if (n_field > 2) then
           ! B_PARALLEL 
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_aperp',10,3,0)
        endif

        ! E_PARALLEL
        if (eparallel_plot_flag == 1) then
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_epar',10,n_field+1,0)
        endif

     endif

     if (plot_n_flag == 1) then

        ! DENSITY
        if (electron_method /= 3) then
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_ion',10,5,1)
        endif
        if (electron_method > 1) then 
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_elec',10,5,indx_e)
        endif
     endif

     if (plot_e_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_ion',10,6,1)
        endif
        if (electron_method > 1) then 
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_elec',10,6,indx_e)
        endif
     endif

     if (plot_v_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_ion',10,7,1)
        endif
        if (electron_method > 1) then 
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_elec',10,7,indx_e)
        endif
     endif

     !-----------------------------------------------------------------
     ! Distribution function data:
     !
     if (n_proc == 1 .and. n_n == 1 .and. dist_print == 1) then
        call gyro_write_h(trim(path)//'out.gyro.hp',trim(path)//'out.gyro.ht',10,11)
     endif
     !-----------------------------------------------------------------

     if (i_proc == 0 .and. lindiff_method > 1) then

        call write_local_real( &
             trim(path)//'out.gyro.diff',10,size(diff),diff) 
        call write_local_real( &
             trim(path)//'out.gyro.diff_i',10,size(diff_i),diff_i)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'out.gyro.diff_trapped',&
                10,size(diff_trapped),diff_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.diff_i_trapped',&
                10,size(diff_i_trapped),diff_i_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_trapped',&
                10,size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_i_trapped',&
                10,size(gbflux_i_trapped),gbflux_i_trapped)
        endif

     endif

     if (lindiff_method >= 4) then

        call write_distributed_real(&
             trim(path)//'out.gyro.diff_n',&
             10,&
             size(diff_n),&
             diff_n)

        call write_distributed_real(&
             trim(path)//'out.gyro.phi_squared_QL_n',&
             10,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n)

        call write_distributed_real(&
             trim(path)//'out.gyro.g_squared_QL_n',&
             10,&
             size(g_squared_QL_n),&
             g_squared_QL_n)

        call write_distributed_real(&
             trim(path)//'out.gyro.gbflux_n',&
             10,&
             size(gbflux_n),&
             gbflux_n)

     endif

     !=============
     ! END LINEAR 
     !=============

  else

     !================
     ! BEGIN NONLINEAR 
     !================

     call write_distributed_real(&
          trim(path)//'out.gyro.diff_n',&
          10,&
          size(diff_n),&
          diff_n)

     call write_distributed_real(&
          trim(path)//'out.gyro.gbflux_n',&
          10,&
          size(gbflux_n),&
          gbflux_n)

     if (lindiff_method >= 4) then
        call write_distributed_real(&
             trim(path)//'out.gyro.phi_squared_QL_n',&
             10,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n)
        call write_distributed_real(&
             trim(path)//'out.gyro.g_squared_QL_n',&
             10,&
             size(g_squared_QL_n),&
             g_squared_QL_n)
     endif

     if (nonlinear_transfer_flag == 1) then
        call write_distributed_real(&
             trim(path)//'out.gyro.nl_transfer',&
             10,&
             size(nl_transfer),&
             nl_transfer)
     endif

     if (i_proc == 0) then

        call write_local_real(trim(path)//'out.gyro.field_rms',10,size(ave_phi),ave_phi)

        call write_local_real( &
             trim(path)//'out.gyro.diff',10,size(diff),diff)
        call write_local_real( &
             trim(path)//'out.gyro.diff_i',10,size(diff_i),diff_i)

        call write_local_real( &
             trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'out.gyro.diff_trapped',&
                10,size(diff_trapped),diff_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.diff_i_trapped',&
                10,size(diff_i_trapped),diff_i_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_trapped',10,&
                size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_i_trapped',10,&
                size(gbflux_i_trapped),gbflux_i_trapped)
        endif

        call write_local_real( &
             trim(path)//'out.gyro.zerobar',10,&
             size(field_fluxave),transpose(field_fluxave))

        allocate(a3(n_kinetic,4,n_x))
        do i=1,n_x
           a3(:,1,i) = h0_n(:,i)
           a3(:,2,i) = h0_e(:,i)
           a3(:,3,i) = source_n(:,i)
           a3(:,4,i) = source_e(:,i)
        enddo
        call write_local_real( &
             trim(path)//'out.gyro.source',10,size(a3),a3)
        deallocate(a3)

        call write_local_real( &
             trim(path)//'out.gyro.moment_zero',10,&
             size(moments_zero_plot),moments_zero_plot)

     endif

     !================
     ! END NONLINEAR 
     !================

  endif
  !-------------------------------------------------------------------

  call gyro_write_error(trim(path)//'out.gyro.error',10)

  !------------------------------------------------------------
  ! Entropy diagnostics
  !
  if (entropy_flag == 1) then
     call gyro_entropy 
     if (i_proc == 0) then 
        call write_local_real(&
             trim(path)//'out.gyro.entropy.out',10,size(entropy),entropy)
     endif
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Velocity-space diagnostics
  !
  if (velocity_output_flag == 1) then
     call gyro_nonlinear_flux_velocity
     call write_distributed_real(&
          trim(path)//'out.gyro.flux_velocity',&
          10,&
          size(nonlinear_flux_velocity),&
          nonlinear_flux_velocity)
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Write precision-monitoring data
  !
  call gyro_write_precision(10,sum(abs(gbflux)))
  !------------------------------------------------------------

  !------------------------------------------------------------
  call gyro_write_timers(trim(path)//'out.gyro.timing',10)
  !------------------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_timedata done]'

10 format(t2,a,t24,es9.3)

end subroutine gyro_write_timedata

!===========================================================================

!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of distributed real array.
!------------------------------------------------------

subroutine write_distributed_real(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err,&
       io_control, &
       fmtstr

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  real :: fn_recv(n_fn)
  !------------------------------------------------------

  select case (io_control)

  case(0)

     return 

  case(1)
 
     ! Initial open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     ! Output

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

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

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_real

!===========================================================================

!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err, &
       io_control, &
       fmtstr

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  complex :: fn_recv(n_fn)
  !------------------------------------------------------


  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     ! Append

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

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
                   MPI_DOUBLE_COMPLEX,&
                   i_send,&
                   in,&
                   GYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
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

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Rewind

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_complex

!===========================================================================

!------------------------------------------------------
! write_local_real.f90
!
! PURPOSE:
!  This routine write a vector of nondistributed reals.
!------------------------------------------------------

subroutine write_local_real(datafile,io,n_fn,fn)

  use gyro_globals, only : &
       data_step, &
       io_control, &
       fmtstr

  !---------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  real :: dummy(n_fn)
  !---------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2)

     ! Append

     open(unit=io,file=datafile,status='old',position='append')
     write(io,fmtstr)  fn(:)
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,fmtstr) dummy(:)
     enddo

     endfile(io)
     close(io)

  end select

end subroutine write_local_real

