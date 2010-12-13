!-----------------------------------------------------
! make_nl.essl.f90
!
! PURPOSE:
!  Allocate (and define) selected arrays for 
!  use with direct and FFT methods.
!
! NOTES:
! ** ESSL-specific routine!  
!-------------------------------------------------------

subroutine make_nl

  use gyro_globals
  use gyro_nl_private

  !--------------------------------------------
  implicit none
  !
  complex, dimension(:,:), allocatable :: x_fft
  real, dimension(:,:), allocatable :: y_fft
  !---------------------------------------------

  !----------------------------------------------------
  do nn=-n_max,n_max
     i_p(nn) = abs(nn)+1
     if (nn < 0) then 
        n_p(nn) = -n(i_p(nn))
     else
        n_p(nn) = n(i_p(nn))
     endif
  enddo
  !
  ! Nonlinear coupling coefficient:
  !
  c_nl_i(:) = -rhos_norm*q_s(:)/(r_s(:)*b_unit_s(:))
  !
  ! ... and add the nonuniform grid effect:
  !
  c_nl_i(:) = dr_eodr(:)*c_nl_i(:)
  !----------------------------------------------------

  !----------------------------------------------------
  ! FFT stuff
  !
  n_max_d = (3*n_max)/2+1
  n_fft   = 2*n_max_d+2
  !
  allocate(x_fft(0:n_fft/2,6*n_x))
  allocate(y_fft(0:n_fft-1,6*n_x))

  allocate(aux1_dcrft(22000))
  allocate(aux2_dcrft(20000))
  allocate(aux1_drcft(22000))
  allocate(aux2_drcft(20000))

  CALL DCRFT(1, &
       x_fft, &
       n_fft/2+1, &
       y_fft, &
       n_fft, &
       n_fft, &
       6*n_x, &
       -1, &
       1.0, &
       aux1_dcrft, &
       22000, &
       aux2_dcrft,&
       20000)

  CALL DRCFT(1, &
       y_fft, &
       n_fft, &
       x_fft, &
       n_fft/2+1, &
       n_fft, &
       6*n_x, &
       1, &
       1.0, &
       aux1_drcft, &
       22000, &
       aux2_drcft,&
       20000)

  deallocate(x_fft)
  deallocate(y_fft)
  !----------------------------------------------------

end subroutine make_nl
