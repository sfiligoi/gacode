!------------------------------------------------
! cgyro_write_restart.f90
!
! PURPOSE:
!  This is the master file controlling output of
!  restart data using MPI-IO.
!------------------------------------------------

subroutine cgyro_write_restart

  use mpi
  use cgyro_globals
  use cgyro_io

  !----------------------------------------------
  implicit none
  !----------------------------------------------
  !
  ! Required for MPI-IO: 
  !
  integer :: filemode
  integer :: finfo
  integer :: fhv
  integer :: fstatus(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer(kind=MPI_OFFSET_KIND) :: offset1
  !----------------------------------------------

  ! Print this data on restart steps only; otherwise exit now
  if (mod(i_time,restart_step*print_step) /= 0) return

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_RDWR,MPI_MODE_CREATE)
  disp     = 0

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
       trim(path)//runfile_restart,&
       filemode,&
       finfo,&
       fhv,&
       i_err)

  call MPI_FILE_SET_VIEW(fhv,&
       disp,&
       MPI_COMPLEX16,&
       MPI_COMPLEX16,&
       'native',&
       finfo,&
       i_err)

  offset1 = size(h_x)*i_proc

  call MPI_FILE_WRITE_AT(fhv,&
       offset1,&
       h_x,&
       size(h_x),&
       MPI_COMPLEX16,&
       fstatus,&
       i_err)

  call MPI_FILE_SYNC(fhv,i_err)
  call MPI_FILE_CLOSE(fhv,i_err)

  !---------------------------------------------------------
  ! Dump restart parameters
  !
  if (i_proc == 0) then

     print *,'[Saving restart data]'

     open(unit=io,file=trim(path)//runfile_restart_tag,status='replace')
     write(io,*) i_current
     write(io,fmtstr) t_current
     close(io)

  endif


end subroutine cgyro_write_restart
