! Print fields

subroutine cgyro_write_timedata

  use cgyro_globals

  implicit none

  complex :: a_norm
  integer :: i_field

  do i_field=1,n_field

     call write_distributed_complex(&
          trim(path)//runfile_field(i_field),&
          myio,&
          size(field(:,:,i_field)),&
          field(:,:,i_field))

     if (n_toroidal == 1) then

        a_norm = field(n_radial/2+1,n_theta/2+1,1) 

        call write_balloon(&
             trim(path)//runfile_fieldb(i_field),&
             myio,&
             size(field(:,:,i_field)),&
             field(:,:,i_field)/a_norm)
     endif

  enddo

end subroutine cgyro_write_timedata

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


subroutine write_balloon(datafile,io,n_fn,fn)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
  integer :: ir,it
   !------------------------------------------------------


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

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_balloning_mode done]' 
  endif

end subroutine write_balloon

