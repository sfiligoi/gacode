!------------------------------------------------------
! write_freq.f90 [caller gyro_write_timedata]
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

subroutine write_freq(datafile,io)

  use mpi
  use gyro_globals
  use math_constants

  !--------------------------------------------------
  implicit none 
  !
  complex, dimension(n_x,n_blend) :: freq_loc
  real,    dimension(n_x,n_blend) :: mode_weight
  !
  complex, dimension(2) :: freq
  complex, dimension(n_n,2) :: freq_collect
  complex, dimension(2,n_n) :: dummy
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

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
     endif

     call collect_complex(freq,freq_collect,2)

     do in=1,n_n

        ! Output to screen

        freq_n(:) = freq_collect(in,:)

        if (i_proc == 0) then

           if (silent_flag == 0) print 10,'n =',n(in),&
                'freq =',freq_n(1), &
                'df =',freq_n(2)  

           if (output_flag == 1) then

              ! Output to file

              write(io,20) freq_n(1),freq_n(2)

           endif

        endif

     enddo ! in

     ! Convergence check for single-n simulation:
     ! Halt will occur in gyro_fulladvance

     if (n_n == 1 .and. n_1(in_1) /= 0) then

        freq_err = abs(freq_n(2))

        call MPI_BCAST(freq_err, &
             1,&
             MPI_DOUBLE_PRECISION,&
             0,&
             GYRO_COMM_WORLD,&
             i_err)

     endif

     if (i_proc == 0 .and. output_flag == 1) then
        close(io)
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

  if (debug_flag == 1 .and. i_proc == 0) print *,'[write_freq called]'

10 format(t2,a,i5,2(3x,a,2(es11.4,1x)))
20 format(4(es11.4,1x))

end subroutine write_freq
