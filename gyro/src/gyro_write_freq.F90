!------------------------------------------------------
! gyro_write_freq.f90 [caller gyro_write_timedata]
!
! PURPOSE:
!  This subroutine computes and prints the instantaeous 
!  eigenfrequency.
!
! NOTES:
!
!  Sign convention: exp(-i w t) = exp(-i w_r t) exp(w_i t)
!
!  - ITG modes have w_r < 0
!  - Trapped-electron modes have w_r > 0
!  - Unstable modes have w_i > 0
!------------------------------------------------------

subroutine gyro_write_freq(datafile,io)

  use mpi
  use gyro_globals
  use math_constants
#ifdef HAVE_HDF5
  use hdf5_api
#endif

  !--------------------------------------------------
  implicit none 
  !
  complex, dimension(n_x,n_blend) :: freq_loc
  real,    dimension(n_x,n_blend) :: mode_weight
  !
  complex, dimension(2) :: freq 
  complex, dimension(n_n,2) :: dummy
  !
  real :: df_r
  real :: df_i
  real :: total_weight
  !
  integer :: data_loop
  !  
  integer, intent(in) :: io  
  character (len=*), intent(in) :: datafile
  !--------------------------------------------------

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

     !-------------------------------------------------------
     ! First, calculate freq(1) (frequencies) and freq(2) 
     ! (deviations).
     !
     freq = (0.0,0.0)
     !
     if (abs(n_1(in_1)) > 99999) then
        call catch_error('|n| too large.')
     endif

     if (n_1(in_1) /= 0) then

        total_weight = 0.0

        if (minval(abs(field_blend_old)) == 0.0) then

           freq(1) = 0.0
           freq(2) = 1.0

        else

           do i=1,n_x
              do j=1,n_blend
                 mode_weight(i,j) = abs(field_blend(j,i,1)) 
                 total_weight = total_weight+mode_weight(i,j)
                 freq_loc(i,j) = &
                      (i_c/dt)*log(field_blend(j,i,1)/field_blend_old(j,i,1))
                 freq(1) = freq(1)+freq_loc(i,j)*mode_weight(i,j)
              enddo
           enddo

           freq(1) = freq(1)/total_weight

           df_r = 0.0
           df_i = 0.0

           do i=1,n_x
              do j=1,n_blend
                 df_r = df_r+ &
                      abs(real(freq_loc(i,j)-freq(1)))*mode_weight(i,j)
                 df_i = df_i+ &
                      abs(aimag(freq_loc(i,j)-freq(1)))*mode_weight(i,j)
              enddo
           enddo

           ! Want a fractional, not absolute, error
           freq(2) = (df_r+i_c*df_i)/total_weight/abs(freq(1)) 

        endif

     else

        freq(:) = 0.0

     endif
     !
     !---------------------------------------------------------------

     call collect_complex(freq,omega_linear,2)

! I/O output
     if(io_method < 3) then 
       if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
          if (output_flag == 1) then
          ! Output to file
            do in=1,n_n
             write(io,20) omega_linear(in,:)
            enddo ! in
            close(io)
          endif !output_flag
        endif !i_proc
      endif !io_method

#ifndef HAVE_HDF5
    if(io_method > 1) then
       if (i_proc == 0) then
        write(*,*) "This gyro was not built with HDF5, please use IO_METHOD=1."
        !**need the exit call here
       endif
    endif
#else
    if(io_method >1 ) then
      if (i_proc == 0 ) then
        write(*,*) "need to put things here to for add_h5 (blah)"
      endif
    endif
#endif
    
      


     ! Convergence check for single-n simulation:
     ! Halt will occur in gyro_fulladvance

     if (n_n == 1 .and. n_1(in_1) /= 0) then

        freq_err = abs(omega_linear(1,2))

        call MPI_BCAST(freq_err, &
             1,&
             MPI_DOUBLE_PRECISION,&
             0,&
             GYRO_COMM_WORLD,&
             i_err)

     endif

  
  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')

        do data_loop=0,data_step
           read(io,20) dummy
        enddo

        endfile(io)
        close(io)

     endif

  end select

  if (debug_flag == 1 .and. i_proc == 0) print *,'[gyro_write_freq called]'

20 format(4(es11.4,1x))

end subroutine gyro_write_freq
