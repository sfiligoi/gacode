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

  ! Print this data on restart steps only; otherwise exit now
  if (mod(i_time,restart_step*print_step) /= 0) return

  if (restart_format==2) then
    if (mpiio_num_files==1) then
       call cgyro_write_restart_one
    else
       call cgyro_write_restart_v2
    endif
  elseif (restart_format==1) then
    if (n_chunk==1) then
       call cgyro_write_restart_one
    else
       call cgyro_write_restart_v1
    endif
  else
    call cgyro_error('ERROR: (CGYRO) Invalid restart_format found.')
  endif
end subroutine cgyro_write_restart

! Dump restart parameters to ASCII file.
subroutine cgyro_write_restart_tag
  use mpi
  use cgyro_globals
  use cgyro_io

  open(unit=io,file=trim(path)//runfile_restart_tag,status='replace')
  write(io,*) i_current
  write(io,fmtstr) t_current
  close(io)

  if (restart_format/=1) then
     open(unit=io,file=trim(path)//runfile_restart_tag_version,status='replace')
     write(io,*) restart_format
     close(io)
  endif

end subroutine cgyro_write_restart_tag

subroutine cgyro_write_restart_one

  use mpi
  use cgyro_globals
  use cgyro_io

  !----------------------------------------------
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
  !----------------------------------------------

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*i_proc
  if (offset1<0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_write_restart_one')
     return
  endif

  ! TODO Error handling
  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_INFO_SET(finfo,"striping_factor",mpiio_stripe_str,i_err)

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

  call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_WRITE_AT in cgyro_write_restart_one failed')
     return
  endif

  call MPI_FILE_SYNC(fhv,i_err)
  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

  !---------------------------------------------------------
  ! Dump restart parameters to ASCII file.
  !
  if (i_proc == 0) then
     call cgyro_write_restart_tag
  endif
  !---------------------------------------------------------

end subroutine cgyro_write_restart_one


subroutine cgyro_write_restart_v2

  use mpi
  use cgyro_globals
  use cgyro_io

  !----------------------------------------------
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
  !----------------------------------------------


  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*i_proc_restart_io

  if (offset1<0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_write_restart_v2')
     return
  endif

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_INFO_SET(finfo,"striping_factor",mpiio_stripe_str,i_err)

  call MPI_FILE_OPEN(NEW_COMM_RESTART_IO,&
          trim(path)//runfile_restart//rtag_v2(i_group_restart_io),&
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

  call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_WRITE_AT in cgyro_write_restart_v2 failed')
     return
  endif

  call MPI_FILE_SYNC(fhv,i_err)
  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

  !---------------------------------------------------------
  ! Wait for all groups to finish
  call MPI_Barrier(CGYRO_COMM_WORLD,i_err)

  !---------------------------------------------------------
  ! Dump restart parameters to ASCII file.
  !
  if (i_proc == 0) then
     call cgyro_write_restart_tag

     ! remove ambiguity
     if (input_restart_format==1) then
        call cgyro_del_restart_v1
        input_restart_format=0 ! avoid doing it in the next iteration
     endif
  endif
  !---------------------------------------------------------

end subroutine cgyro_write_restart_v2

subroutine cgyro_del_restart_v2
  use mpi
  use cgyro_globals
  use cgyro_io

  integer :: j
  integer :: stat

  do j=1,mpiio_num_files
    open(unit=io,file=trim(path)//runfile_restart//rtag_v2(j),iostat=stat,status='old')
    if (stat == 0) close(io, status='delete')
  enddo
end subroutine cgyro_del_restart_v2

! --------------------------------------------------------------
! cgyro_write_restart_v1 is deprecated (As of Aug 2017)
! Leaving in the code for now for testing/debugging reasons
! Should be removed in the near future

subroutine cgyro_write_restart_v1

  use mpi
  use cgyro_globals
  use cgyro_io

  !----------------------------------------------
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
  !----------------------------------------------

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  n_loc = nc/n_chunk
  n_remain = nc-n_loc*n_chunk
  allocate(h_chunk(n_loc+n_remain,nv_loc))

  offset1 = size(h_chunk,kind=MPI_OFFSET_KIND)*i_proc

  if (offset1<0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_write_restart_v1')
     return
  endif

  do j=1,n_chunk

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_INFO_SET(finfo,"striping_factor",mpiio_stripe_str,i_err)

     call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart//rtag(j),&
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

     i1 = 1+(j-1)*n_loc
     i2 = j*n_loc
     if (j == n_chunk) then
        i2 = nc
        n_loc = n_loc+n_remain
     endif

     h_chunk(1:n_loc,:) = h_x(i1:i2,:)

     call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h_chunk,&
          size(h_chunk),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)
     if (i_err /= 0) then
       call cgyro_error('ERROR: (CGYRO) MPI_FILE_WRITE_AT in cgyro_write_restart_v1 failed')
        return
     endif

     call MPI_FILE_SYNC(fhv,i_err)
     call MPI_FILE_CLOSE(fhv,i_err)
     call MPI_INFO_FREE(finfo,i_err)

  enddo

  deallocate(h_chunk)

  !---------------------------------------------------------
  ! Dump restart parameters to ASCII file.
  !
  if (i_proc == 0) then
     call cgyro_write_restart_tag

     ! remove ambiguity
     if (input_restart_format==2) then
        call cgyro_del_restart_v2
        input_restart_format=0 ! avoid doing it in the next iteration
     endif
  endif
  !---------------------------------------------------------

end subroutine cgyro_write_restart_v1

subroutine cgyro_del_restart_v1
  use mpi
  use cgyro_globals
  use cgyro_io

  integer :: j
  integer :: stat

  do j=1,n_chunk
    open(unit=io,file=trim(path)//runfile_restart//rtag(j),iostat=stat,status='old')
    if (stat == 0) close(io, status='delete')
  enddo
end subroutine cgyro_del_restart_v1

