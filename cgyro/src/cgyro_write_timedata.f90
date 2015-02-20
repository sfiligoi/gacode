!-----------------------------------------------------------------
! cgyro_write_timedata.f90
!
! PURPOSE:
!  Output of time-dependent data 
!-----------------------------------------------------------------

subroutine cgyro_write_timedata

  use mpi
  use cgyro_globals

  implicit none

  complex :: a_norm
  integer :: i_field
  logical :: lfe
  real :: vfreq(2)

  ! Print this data on print steps only; otherwise exit now
  if (mod(i_time,print_step) /= 0) return

  if (nonlinear_flag == 1) then

     ! Density flux
     call write_distributed_real(&
          trim(path)//runfile_flux(1),&
          size(flux(:,1)),&
          flux(:,1))
     ! Energy flux
     call write_distributed_real(&
          trim(path)//runfile_flux(2),&
          size(flux(:,2)),&
          flux(:,2))
  endif

  do i_field=1,n_field

     ! Complete field output 
     call write_distributed_complex(&
          trim(path)//runfile_field(i_field),&
          size(field(:,:,i_field)),&
          field(:,:,i_field))

     ! Field intensity
     call write_distributed_real(&
          trim(path)//runfile_power(i_field),&
          size(power(:,i_field)),&
          power(:,i_field))

     if (n_toroidal == 1 .and. n > 0) then

        ! Ballooning mode output for linear runs with a single mode

        a_norm = field(n_radial/2+1,n_theta/2+1,1) 

        call write_balloon(&
             trim(path)//runfile_fieldb(i_field),&
             field(:,:,i_field)/a_norm)
     endif

  enddo
 
  ! Linear frequency diagnostics for every value of n
  call cgyro_freq
  vfreq(1) = real(freq) 
  vfreq(2) = aimag(freq)
  call write_distributed_real(trim(path)//runfile_freq,size(vfreq),vfreq)
  if (n_toroidal == 1) call write_freq()

  call write_time(trim(path)//runfile_time)

  ! Check for manual halt signal
  if (i_proc == 0) then
     inquire(file=trim(path)//'halt',exist=lfe)
     if (lfe .eqv. .true.) then
        open(unit=io,file='halt',status='old')
        read(io,*) signal
        close(io)
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

subroutine write_distributed_complex(datafile,n_fn,fn)

  use mpi
  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
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
        i_group_send = in-1

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


!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of real distributed array.
!------------------------------------------------------

subroutine write_distributed_real(datafile,n_fn,fn)

  use mpi
  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
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
        i_group_send = in-1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   i_send,&
                   in,&
                   CGYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
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

end subroutine write_distributed_real

subroutine write_balloon(datafile,fn)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  complex, intent(in) :: fn(n_radial,n_theta)
  !
  integer :: ir,jr,it,np
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

     np = n_radial/2/box_size

     do ir=-np,np-1
        do it=1,n_theta
           jr = box_size*ir+n_radial/2+1
           f_balloon(ir+np+1,it) = fn(jr,it) &
                *exp(-2*pi*i_c*ir*k_theta*rmin)
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

subroutine write_time(datafile)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
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

     if (n_toroidal > 1 .and. i_proc == 0) then
         print '(a,1pe9.3,a,5(1pe9.3,1x))',&
                '[t = ',i_time*delta_t,&
                '] t_err: ',field_error
     endif

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

subroutine write_freq()

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !------------------------------------------------------
 
  select case (io_control)

  case(0,1,3)

     return

  case(2)

     ! Append

     if (i_proc == 0) then
        print '(a,1pe9.3,a,1pe10.3,1pe10.3,a,1pe9.3,a,5(1pe9.3,1x))',&
             '[t = ',i_time*delta_t,&
             '][w = ',freq,&
             '][dw = ',abs(freq_err),&
             '] t_err: ',field_error
     endif

  end select

end subroutine write_freq

!====================================================================================

subroutine write_distribution(datafile,indx)

  use mpi

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  integer, intent(in) :: indx
  character (len=*), intent(in) :: datafile
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
     if (indx == 1) then
        call MPI_GATHER(h_x(:,:),&
             size(h_x),&
             MPI_DOUBLE_COMPLEX,&
             h_x_glob(:,:),&
             size(h_x),&
             MPI_DOUBLE_COMPLEX,&
             0,&
             NEW_COMM_1,&
             i_err)
     else
        call MPI_GATHER(cap_h_c(:,:),&
             size(cap_h_c),&
             MPI_DOUBLE_COMPLEX,&
             h_x_glob(:,:),&
             size(h_x),&
             MPI_DOUBLE_COMPLEX,&
             0,&
             NEW_COMM_1,&
             i_err)
     endif

     if (i_proc == 0) then
        do iv=1,nv
           if (box_size == 1) then 
              do ic=1,nc
                 f_balloon(ir_c(ic),it_c(ic)) = h_x_glob(ic,iv) &
                      *exp(-2*pi*i_c*px(ir_c(ic))*k_theta*rmin)
              enddo
           else
              f_balloon(:,:) = 0.0
           endif
           write(io,fmtstr) transpose(f_balloon(:,:))
        enddo
        close(io)
     endif

     deallocate(h_x_glob)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     endfile(io)
     close(io)

  end select

end subroutine write_distribution
