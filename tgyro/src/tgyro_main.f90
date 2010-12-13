!-----------------------------------------------------------------
! Copyright 1999-2009 by Jeff Candy, Ron Waltz and Mark Fahey.
!-----------------------------------------------------------------

program tgyro_main

  use tgyro_globals
  use gyro_globals, only : path, GYRO_COMM_WORLD

  implicit none

  include 'mpif.h'

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc_global,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc_global,ierr)
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Get INPUT paths and broadcast them
  !
  call tgyro_read_input   
  !
  ! At this stage, information is only on process 0.
  !
  call tgyro_comm_setup
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Set GYRO path to local path (per instance). 
  path = lpath
  !
  ! Set GYRO communicator to local communicator (per instance).
  GYRO_COMM_WORLD = gyro_comm
  !-----------------------------------------------------------------

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  select case (transport_method)

  case (1)

     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile),position='append')
        write(1,*) 'TRANSPORT_METHOD=1 deprecated'
        close(1)
     endif

  case (2)

     ! Local transport model

     call tgyro_iteration_driver 

     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile),position='append')
        write(1,*) error_msg
        close(1)
        call system('[ -f tgyro_osx_exec ] && chmod +x tgyro_osg_exec && tgyro_osg_exec')  
     endif

  case (3)

     ! Multi-job utility

     call tgyro_multi_driver

  end select

  call MPI_FINALIZE(ierr)

end program tgyro_main
