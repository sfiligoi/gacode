!------------------------------------------------------
! cgyro_read_restart.f90
!
! PURPOSE:
!  This is the master file controlling the restart
!  via MPI-IO.
!------------------------------------------------------

subroutine cgyro_read_restart

  use mpi
  use cgyro_globals
  use cgyro_io

  !---------------------------------------------------
  implicit none
  !
  ! Required for MPI-IO: 
  !
  integer :: filemode
  integer :: finfo
  integer :: fhv
  integer :: fstatus(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer(kind=MPI_OFFSET_KIND) :: offset1
  !---------------------------------------------------

  call cgyro_info('Reading MPI-IO restart data.')

  !---------------------------------------------------------
  ! Read restart parameters from ASCII file.
  !
  call cgyro_info('Writing MPI-IO restart data.')
  if (i_proc == 0) then

     open(unit=io,&
          file=trim(path)//runfile_restart_tag,&
          status='old')

     read(io,*) i_current
     read(io,fmtstr) t_current
     close(io)

  endif

  ! Broadcast to all cores.

  call MPI_BCAST(i_current,&
       1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

  call MPI_BCAST(t_current,&
       1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)
  !---------------------------------------------------------

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

  call MPI_FILE_READ_AT(fhv,&
       offset1,&
       h_x,&
       size(h_x),&
       MPI_COMPLEX16,&
       fstatus,&
       i_err)

  call MPI_FILE_CLOSE(fhv,i_err)

end subroutine cgyro_read_restart
