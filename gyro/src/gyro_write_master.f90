!------------------------------------------------------
! gyro_write_master.f90
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

subroutine gyro_write_master(action)

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  integer :: mode
  integer, intent(in) :: action
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
  real, dimension(:,:), allocatable :: a2
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: n_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: e_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: v_plot
  !---------------------------------------------------

  !---------------------------------------------------
  ! Determine file control mode:
  !
  ! mode = 1 -> file create 
  !      = 2 -> file write 
  !      = 3 -> file reposition 
  !      = 4 -> file write on step = 0
  !
  ! action = 1 -> simulation still initializing
  !        = 2 -> simulation running
  !
  if (action == 1) then

     ! File creation or repositioning:

     select case (restart_method)

     case(-1) 

        ! No use of restart facility
        mode = 1

     case(0,2) 

        ! Start of restartable simulation
        mode = 1

     case(1)

        ! Continuation of restartable simulation        
        mode = 3

     end select

  else

     ! File writing

     if (step == 0) then
        mode = 4
     else
        mode = 2
     endif

  endif

  if (output_flag == 0) then
     mode = -mode
  endif
  !---------------------------------------------------

  !---------------------------------------------------
  ! Timestep data:
  !
  if (i_proc == 0) then
     call write_step(trim(path)//'t.out',10,mode)
  endif
  !---------------------------------------------------

  call proc_time(cp0)

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
          trim(path)//'u.out',&
          10,&
          mode,&
          size(phi_plot(:,:,1:n_field)),&
          phi_plot(:,:,1:n_field))

  endif

  if (plot_n_flag == 1) then

     ! DENSITY

     call write_distributed_complex(&
          trim(path)//'moment_n.out',&
          10,&
          mode,&
          size(n_plot),&
          n_plot)

  endif

  if (plot_e_flag == 1) then

     ! ENERGY

     call write_distributed_complex(&
          trim(path)//'moment_e.out',&
          10,&
          mode,&
          size(e_plot),&
          e_plot)

  endif

  if (plot_v_flag == 1) then

     ! PARALLEL VELOCITY

     call write_distributed_complex(&
          trim(path)//'moment_v.out',&
          10,&
          mode,&
          size(v_plot),&
          v_plot)

  endif

  !--------------------------------------------------

  !--------------------------------------------------
  ! Output of field at r=r0:
  !
  if (field_r0_flag == 1) then
     call write_distributed_complex(&
          trim(path)//'field_r0.out',&
          10,&
          mode,&
          size(field_r0_plot),&
          field_r0_plot)
  endif
  !--------------------------------------------------

  call proc_time(cp1)
  call get_field_spectrum
  call write_distributed_real(&
       trim(path)//'kxkyspec.out',&
       10,&
       mode,&
       size(kxkyspec),&
       kxkyspec)
  call proc_time(cp2)

  if (i_proc == 0) then
     call write_local_real(&
          trim(path)//'k_perp_squared.out',&
          10,&
          mode,&
          size(k_perp_squared),&
          k_perp_squared)
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
  call get_field_fluxave
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

     call write_freq(trim(path)//'freq.out',10,mode)

     if (plot_u_flag == 1) then        

        ! PHI
        call gyro_ballooning_mode(trim(path)//'balloon_phi.out',10,mode,1,0)

        if (n_field > 1) then
           ! A_PARALLEL 
           call gyro_ballooning_mode(trim(path)//'balloon_a.out',10,mode,2,0)
        endif

        if (n_field > 2) then
           ! B_PARALLEL 
           call gyro_ballooning_mode(trim(path)//'balloon_aperp.out',10,mode,3,0)
        endif

        ! E_PARALLEL
        if (eparallel_plot_flag == 1) then
           call gyro_ballooning_mode(trim(path)//'balloon_epar.out',10,mode,n_field+1,0)
        endif

     endif

     if (plot_n_flag == 1) then

        ! DENSITY
        if (electron_method /= 3) then
           call gyro_ballooning_mode(trim(path)//'balloon_n_ion.out',10,mode,5,1)
        endif
        if (electron_method > 1) then 
           call gyro_ballooning_mode(trim(path)//'balloon_n_elec.out',10,mode,5,indx_e)
        endif
     endif

     if (plot_e_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then
           call gyro_ballooning_mode(trim(path)//'balloon_e_ion.out',10,mode,6,1)
        endif
        if (electron_method > 1) then 
           call gyro_ballooning_mode(trim(path)//'balloon_e_elec.out',10,mode,6,indx_e)
        endif
     endif

     if (plot_v_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then
           call gyro_ballooning_mode(trim(path)//'balloon_v_ion.out',10,mode,7,1)
        endif
        if (electron_method > 1) then 
           call gyro_ballooning_mode(trim(path)//'balloon_v_elec.out',10,mode,7,indx_e)
        endif
     endif

     !-----------------------------------------------------------------
     ! Distribution function data:
     !
     if (n_proc == 1 .and. n_n == 1 .and. dist_print == 1) then
        call write_h(trim(path)//'hp.out',trim(path)//'ht.out',10,11,mode)
     endif
     !-----------------------------------------------------------------

     if (i_proc == 0 .and. lindiff_method > 1) then

        call write_local_real( &
             trim(path)//'diff.out',10,mode,size(diff),diff) 
        call write_local_real( &
             trim(path)//'diff_i.out',10,mode,size(diff_i),diff_i)
        call write_local_real( &
             trim(path)//'gbflux.out',10,mode,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'gbflux_mom.out',10,mode,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'gbflux_i.out',10,mode,size(gbflux_i),gbflux_i)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'diff_trapped.out',&
                10,mode,size(diff_trapped),diff_trapped)
           call write_local_real( &
                trim(path)//'diff_i_trapped.out',&
                10,mode,size(diff_i_trapped),diff_i_trapped)
           call write_local_real( &
                trim(path)//'gbflux_trapped.out',&
                10,mode,size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'gbflux_i_trapped.out',&
                10,mode,size(gbflux_i_trapped),gbflux_i_trapped)
        endif

     endif

     if (lindiff_method >= 4) then

        call write_distributed_real(&
             trim(path)//'diff_n.out',&
             10,&
             mode,&
             size(diff_n),&
             diff_n)

        call write_distributed_real(&
             trim(path)//'phi_squared_QL_n.out',&
             10,&
             mode,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n)

        call write_distributed_real(&
             trim(path)//'g_squared_QL_n.out',&
             10,&
             mode,&
             size(g_squared_QL_n),&
             g_squared_QL_n)

        call write_distributed_real(&
             trim(path)//'gbflux_n.out',&
             10,&
             mode,&
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

     call proc_time(cp7)

     call write_distributed_real(&
          trim(path)//'diff_n.out',&
          10,&
          mode,&
          size(diff_n),&
          diff_n)

     call write_distributed_real(&
          trim(path)//'gbflux_n.out',&
          10,&
          mode,&
          size(gbflux_n),&
          gbflux_n)

     if (lindiff_method >= 4) then
        call write_distributed_real(&
             trim(path)//'phi_squared_QL_n.out',&
             10,&
             mode,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n)
        call write_distributed_real(&
             trim(path)//'g_squared_QL_n.out',&
             10,&
             mode,&
             size(g_squared_QL_n),&
             g_squared_QL_n)
     endif

     if (nonlinear_transfer_flag == 1) then
        call write_distributed_real(&
             trim(path)//'nonlinear_transfer_n.out',&
             10,&
             mode,&
             size(Tr_p),&
             Tr_p)
        call write_distributed_real(&
             trim(path)//'turbulent_energy_n.out',&
             10,&
             mode,&
             size(Eng_p),&
             Eng_p)
     endif

     call proc_time(cp8)

     if (i_proc == 0) then

        call write_local_real(trim(path)//'field_rms.out',10,mode,size(ave_phi),ave_phi)

        call write_local_real( &
             trim(path)//'diff.out',10,mode,size(diff),diff)
        call write_local_real( &
             trim(path)//'diff_i.out',10,mode,size(diff_i),diff_i)

        call write_local_real( &
             trim(path)//'gbflux.out',10,mode,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'gbflux_mom.out',10,mode,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'gbflux_i.out',10,mode,size(gbflux_i),gbflux_i)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'diff_trapped.out',&
                10,mode,size(diff_trapped),diff_trapped)
           call write_local_real( &
                trim(path)//'diff_i_trapped.out',&
                10,mode,size(diff_i_trapped),diff_i_trapped)
           call write_local_real( &
                trim(path)//'gbflux_trapped.out',10,mode,&
                size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'gbflux_i_trapped.out',10,mode,&
                size(gbflux_i_trapped),gbflux_i_trapped)
        endif

        allocate(a2(3,n_x))
        a2(1,:) = phi_fluxave(:) 
        a2(2,:) = a_fluxave(:)
        a2(3,:) = aperp_fluxave(:)
        call write_local_real( &
             trim(path)//'zerobar.out',10,mode,size(a2),a2)
        deallocate(a2)

        allocate(a3(n_kinetic,4,n_x))
        do i=1,n_x
           a3(:,1,i) = h0_n(:,i)
           a3(:,2,i) = h0_e(:,i)
           a3(:,3,i) = source_n(:,i)
           a3(:,4,i) = source_e(:,i)
        enddo
        call write_local_real( &
             trim(path)//'source.out',10,mode,size(a3),a3)
        deallocate(a3)

        call write_local_real( &
             trim(path)//'moment_zero.out',10,mode,&
             size(moments_zero_plot),moments_zero_plot)

     endif

     !================
     ! END NONLINEAR 
     !================

  endif
  !-------------------------------------------------------------------

  call write_error(trim(path)//'error.out',10,mode)

  !------------------------------------------------------------
  ! Entropy diagnostics
  !
  if (entropy_flag == 1) then
     call gyro_entropy 
     if (i_proc == 0) then 
        call write_local_real(&
             trim(path)//'entropy.out',10,mode,size(entropy),entropy)
     endif
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Velocity-space diagnostics
  !
  if (velocity_output_flag == 1) then
     call get_nonlinear_flux_velocity
     call write_distributed_real(&
          trim(path)//'flux_velocity.out',&
          10,&
          mode,&
          size(nonlinear_flux_velocity),&
          nonlinear_flux_velocity)
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Write precision-monitoring data
  !
  call gyro_write_precision(10,sum(abs(gbflux)))
  !------------------------------------------------------------

  !--------------------------------------
  ! ** Timer diagnostics last **
  !
  call proc_time(CPU_diag_outp)
  CPU_diag_b = CPU_diag_b + (CPU_diag_outp - CPU_diag_mid)
  call write_timing(trim(path)//'timing.out',10,mode)
  CPU_diag_mid = CPU_diag_outp
  !--------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_master done]'

10 format(t2,a,t24,es9.3)

end subroutine gyro_write_master

!===========================================================================

!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of distributed real array.
!------------------------------------------------------

subroutine write_distributed_real(datafile,io,action,n_fn,fn)

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
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: action
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

  include 'mpif.h'

  !-----------------------------
  ! action = 1 -> open file
  !        = 2 -> write to file
  !        = 3 -> reposition
  !-----------------------------

  select case (action)

  case(1)
 
     ! Initial open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2,4)

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

           write(io,10) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,10) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

10 format(es11.4,1x)

end subroutine write_distributed_real

!===========================================================================

!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex(datafile,io,action,n_fn,fn)

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
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: action
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

  include 'mpif.h'

  !-----------------------------
  ! action = 1 -> open file
  !        = 2 -> write to file
  !        = 3 -> reposition
  !-----------------------------

  select case (action)

  case(1)

     ! Initial open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2,4)

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

           write(io,10) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,10) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

10 format(2(es11.4,1x))

end subroutine write_distributed_complex

!===========================================================================

!------------------------------------------------------
! write_local_real.f90
!
! PURPOSE:
!  This routine write a vector of nondistributed reals.
!------------------------------------------------------

subroutine write_local_real(datafile,io,action,n_fn,fn)

  use gyro_globals, only : &
       data_step

  !---------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: action
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  real :: dummy(n_fn)
  !---------------------------------------------------

  select case (action)

  case(1)

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2,4)

     open(unit=io,file=datafile,status='old',position='append')
     write(io,'(20(1pe15.8,1x))')  fn(:)
     close(io)

  case(3)

     ! Reposition after restart

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,*) dummy(:)
     enddo

     endfile(io)
     close(io)

  end select

end subroutine write_local_real

