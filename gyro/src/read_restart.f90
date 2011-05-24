!--------------------------------------------------------
! read_restart.f90 [caller: BigScience]
!
! PURPOSE:
!  This is the master file controlling the restart.
!--------------------------------------------------------

subroutine read_restart

  use mpi
  use gyro_globals
  use gyro_pointers

  !------------------------------------------------------------
  ! Local variables:
  !
  implicit none
  !
  integer :: n_proc_old
  integer :: io
  !
  complex, dimension(n_stack,n_x,n_nek_loc_1,n_kinetic) :: h_in
  !------------------------------------------------------------

  io = io_restart

  !-------------------------------------------
  ! Set the step number to zero for gyro_write_master, 
  ! even if this is a restart:
  !
  step = 0
  !-------------------------------------------

  !---------------------------------------------------------------------
  ! Check to see if we are asking for a restart
  ! when there is no restart data available:
  !
  if (i_proc == 0) then

     open(unit=io,&
          file=trim(path)//file_tag_restart(restart_new_flag),&
          status='old',iostat=i_err)

     if (i_err /= 0 .and. restart_method > 0) then
        call send_message(&
             'INFO: Restart data not available.  Reseting restart_method.')
        restart_method = 0
     endif

  endif
  !---------------------------------------------------------------------

  ! Ensure sync.
  call MPI_BCAST(restart_method,1,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,i_err)

  select case(restart_method)

  case (-1)

     !------------------------
     ! Ignore restart facility
     !------------------------

     call make_initial_h

     t_current = 0.0
     i_restart = 0
     data_step = 0

  case (0)

     !---------------------
     ! New simulation block
     !---------------------

     call make_initial_h

     t_current = 0.0
     i_restart = 0
     data_step = 0

     ! i_restart is always tagged to most recent output 
     ! files.  If restart_new_flag = 0, then we want 1-i_restart

  case (1,2)

     !-------------------------------------------
     ! Restart block
     !
     ! 1=restart, continue writing restart files
     ! 2=restart, but do not write restart files
     !-------------------------------------------

     ! Get restart values from last run

     if (i_proc == 0) then

        open(unit=io,&
             file=trim(path)//file_tag_restart(restart_new_flag),&
             status='old')

        read(io,*) data_step
        read(io,10) t_current
        read(io,*) n_proc_old
        read(io,*) i_restart
        close(io)

     endif

     call MPI_BCAST(data_step,&
          1,MPI_INTEGER,0,GYRO_COMM_WORLD,i_err)

     call MPI_BCAST(t_current,&
          1,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,i_err)

     call MPI_BCAST(n_proc_old,&
          1,MPI_INTEGER,0,GYRO_COMM_WORLD,i_err)

     call MPI_BCAST(i_restart,&
          1,MPI_INTEGER,0,GYRO_COMM_WORLD,i_err)

     ! Reset time if running TGYRO local method
     if (transport_method == 2) then
        t_current = 0.0
        data_step = 0
     endif

     ! Trap error for incorrect number of processors
     if (n_proc_old /= n_proc) then 
        call catch_error('Processor number changed.')
     endif

     if (i_proc == 0) then
        open(unit=io,file=trim(path)//file_restart(i_restart),status='old')
        print '(t2,a)','Restarting: '
        read(io,10) h
        print '(a,i1,a,$)','[',0,']'
     endif

     do i_proc_w=1,n_proc-1

        if (i_proc == 0) then

           if (i_proc_w < 10) then
              print '(a,i1,a,$)','[',i_proc_w,']'
           else if (i_proc_w < 100) then
              print '(a,i2,a,$)','[',i_proc_w,']'
           else
              print '(a,i3,a,$)','[',i_proc_w,']'
           endif
           if (modulo(i_proc_w,8) == 0) print *

           read(io,10) h_in

           call MPI_SEND(h_in,&
                size(h_in),&
                MPI_DOUBLE_COMPLEX,&
                i_proc_w,&
                i_proc_w,&
                GYRO_COMM_WORLD,&
                i_err)

        else if (i_proc == i_proc_w) then

           call MPI_RECV(h,&
                size(h),&
                MPI_DOUBLE_COMPLEX,&
                0,&
                i_proc_w,&
                GYRO_COMM_WORLD,&
                recv_status,&
                i_err)

        endif

     enddo

     if (i_proc == 0) then
        print *
        close(io)
     endif

  case default

     call catch_error('Bad value for restart_method.')

  end select

  ! ** Regenerate fields:

  call get_field_explicit

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[read_restart done]'
  endif

  ! ** Keep this consistent with write_restart.f90

10 format(2(es11.4,1x))

end subroutine read_restart
