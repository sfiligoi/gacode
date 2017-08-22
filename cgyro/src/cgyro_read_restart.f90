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

  !---------------------------------------------------------
  ! Read restart parameters from ASCII file.
  !
  if (i_proc == 0) then

     open(unit=io,&
          file=trim(path)//runfile_restart_tag,&
          status='old')

     read(io,*) i_current
     read(io,fmtstr) t_current
     close(io)

     open(unit=io,&
          file=trim(path)//runfile_restart_tag_version,&
          iostat=i_err,&
          status='old')

     if (i_err==0) then
          read(io,*) input_restart_format
          close(io)
          call cgyro_info('Restart version found.')
     else
          ! v1 had no version file
          input_restart_format = 1
          call cgyro_info('Restart version not found. Assuming v1.')
     endif

  endif

  ! Broadcast to all cores.

  call MPI_BCAST(i_current,&
       1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

  call MPI_BCAST(t_current,&
       1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)

  call MPI_BCAST(input_restart_format,&
       1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

  !---------------------------------------------------------

  if (input_restart_format==2) then
    if (mpiio_num_files==1) then
       call cgyro_read_restart_one
    else
       call cgyro_read_restart_v2
    endif
  elseif (input_restart_format==1) then
    if (n_chunk==1) then
       call cgyro_read_restart_one
    else
       call cgyro_read_restart_v1
    endif
  else
    call cgyro_error('ERROR: (CGYRO) Invalid restart_format found.')
  endif

end subroutine cgyro_read_restart

subroutine cgyro_read_restart_one

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

  filemode = MPI_MODE_RDONLY
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*i_proc
  if (offset1<0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_read_restart_one')
     return
  endif

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart,&
          filemode,&
          finfo,&
          fhv,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_OPEN in cgyro_read_restart_one failed')
     return
  endif

  call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_COMPLEX16,&
          MPI_COMPLEX16,&
          'native',&
          finfo,&
          i_err)

  call MPI_FILE_READ_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_READ_AT in cgyro_read_restart_one failed')
     return
  endif

  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

end subroutine cgyro_read_restart_one

subroutine cgyro_read_restart_v2

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

  filemode = MPI_MODE_RDONLY
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*i_proc_restart_io

  if (offset1<0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_read_restart_v2')
     return
  endif

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_FILE_OPEN(NEW_COMM_RESTART_IO,&
          trim(path)//runfile_restart//rtag_v2(i_group_restart_io),&
          filemode,&
          finfo,&
          fhv,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_OPEN in cgyro_read_restart_v2 failed')
     return
  endif

  call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_COMPLEX16,&
          MPI_COMPLEX16,&
          'native',&
          finfo,&
          i_err)

  call MPI_FILE_READ_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_READ_AT in cgyro_read_restart_v2 failed')
     return
  endif

  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

end subroutine cgyro_read_restart_v2

! ---------------------------------------------------
! This routine is provided for backwards compatibility
! Eventually we should phase it out

subroutine cgyro_read_restart_v1

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
  !
  ! File chunking variables
  !
  integer :: i1,i2,j,n_loc,n_remain
  complex, dimension(:,:), allocatable :: h_chunk
  !---------------------------------------------------

  filemode = MPI_MODE_RDONLY
  disp     = 0

  n_loc = nc/n_chunk
  n_remain = nc-n_loc*n_chunk
  allocate(h_chunk(n_loc+n_remain,nv_loc))

  offset1 = size(h_chunk,kind=MPI_OFFSET_KIND)*i_proc
  
  if (offset1<0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_read_restart_v1')
     return
  endif

  do j=1,n_chunk

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart//rtag(j),&
          filemode,&
          finfo,&
          fhv,&
          i_err)
     if (i_err /= 0) then
        call cgyro_error('ERROR: (CGYRO) MPI_FILE_OPEN in cgyro_read_restart_v1 failed')
        return
     endif

     call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_COMPLEX16,&
          MPI_COMPLEX16,&
          'native',&
          finfo,&
          i_err)

     call MPI_FILE_READ_AT(fhv,&
          offset1,&
          h_chunk,&
          size(h_chunk),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)
     if (i_err /= 0) then
        call cgyro_error('ERROR: (CGYRO) MPI_FILE_READ_AT in cgyro_read_restart_v1 failed')
        return
     endif

     i1 = 1+(j-1)*n_loc
     i2 = j*n_loc
     if (j == n_chunk) then
        i2 = nc
        n_loc = n_loc+n_remain
     endif

     h_x(i1:i2,:) = h_chunk(1:n_loc,:)

     call MPI_FILE_CLOSE(fhv,i_err)
     call MPI_INFO_FREE(finfo,i_err)

  enddo

  deallocate(h_chunk)

end subroutine cgyro_read_restart_v1
