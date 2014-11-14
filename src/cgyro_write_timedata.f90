! Print fields

subroutine cgyro_write_timedata

  use mpi
  use cgyro_globals

  implicit none

  complex :: a_norm
  integer :: i_field
  logical :: lfe

  ! Print this data on print steps
  if (mod(i_time,prin_time) /= 0) return

  do i_field=1,n_field

     call write_distributed_complex(&
          trim(path)//runfile_field(i_field),&
          io_data,&
          size(field(:,:,i_field)),&
          field(:,:,i_field))

     if (n_toroidal == 1) then

        a_norm = field(n_radial/2+1,n_theta/2+1,1) 

        call write_balloon(&
             trim(path)//runfile_fieldb(i_field),&
             io_data,&
             field(:,:,i_field)/a_norm)
     endif

  enddo

  if (n_toroidal == 1) call write_freq(trim(path)//runfile_freq,io_data)

  call write_time(trim(path)//runfile_time,io_data)

  ! Check for manual halt signal
  if (i_proc == 0) then
     inquire(file=trim(path)//'halt',exist=lfe)
     if (lfe .eqv. .true.) then
        open(unit=io_data,file='halt',status='old')
        read(io_data,*) signal
        close(io_data)
     endif
  endif

  call MPI_BCAST(signal,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_write_timedata


!===============================================================================
! Individual parallel-safe I/O routines
!===============================================================================

!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex(datafile,io,n_fn,fn)

  use mpi
  use cgyro_globals

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

     do in=1,n_toroidal

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_toroidal

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   i_send,&
                   in,&
                   CGYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   0,&
                   in,&
                   CGYRO_COMM_WORLD,&
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
        close(io)

     endif

  end select

end subroutine write_distributed_complex


subroutine write_balloon(datafile,io,fn)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  complex, intent(in) :: fn(n_radial,n_theta)
  !
  integer :: ir,it
  !------------------------------------------------------

  if (i_proc > 0) return

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

     ! Construct ballooning-space form of field

     do ir=1,n_radial
        do it=1,n_theta
           f_balloon(ir,it) = fn(ir,it) &
                *exp(-2*pi*i_c*indx_r(ir)*k_theta*rmin)
        enddo
     enddo

     write(io,fmtstr) transpose(f_balloon(:,:))
     close(io)

     !-------------------------------------------------------

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     endfile(io)
     close(io)

  end select

end subroutine write_balloon

!=========================================================================================

subroutine write_time(datafile,io)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  !------------------------------------------------------

  if (i_proc > 0) return

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
     write(io,fmtstr) i_time*delta_t
     close(io)

     !-------------------------------------------------------

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     endfile(io)
     close(io)

  end select

end subroutine write_time

!====================================================================================

subroutine write_freq(datafile,io)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  !------------------------------------------------------

  ! Compute frequencies on all cores
  call cgyro_freq

  if (i_proc > 0) return

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

     ! Construct ballooning-space form of field

     write(io,fmtstr) freq
     print '(t2,1pe10.3,2x,2(1pe13.6,1x),2x,2(1pe10.3,1x))',i_time*delta_t,freq,freq_err
     close(io)

     !-------------------------------------------------------

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     endfile(io)
     close(io)

  end select

end subroutine write_freq

!====================================================================================

subroutine write_distribution(datafile,io)

  use mpi

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  complex, dimension(:,:), allocatable :: h_x_glob
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

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
     endif

     allocate(h_x_glob(nc,nv))

     ! Collect distribution onto process 0
     call MPI_GATHER(h_x(:,:),&
          size(h_x),&
          MPI_DOUBLE_COMPLEX,&
          h_x_glob(:,:),&
          size(h_x),&
          MPI_DOUBLE_COMPLEX,&
          0,&
          NEW_COMM_1,&
          i_err)

     if (i_proc == 0) then
        do iv=1,nv
           do ic=1,nc
              f_balloon(ir_c(ic),it_c(ic)) = h_x_glob(ic,iv) &
                   *exp(-2*pi*i_c*indx_r(ir_c(ic))*k_theta*rmin)
           enddo
           write(io_data,fmtstr) transpose(f_balloon(:,:))
        enddo
        close(io_data)
     endif

     deallocate(h_x_glob)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     endfile(io)
     close(io)

  end select

end subroutine write_distribution
