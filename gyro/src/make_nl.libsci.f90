!-----------------------------------------------------
! make_nl.libsci.f90
!
! PURPOSE:
!  Allocate (and define) selected arrays for 
!  use with direct and FFT methods.
!
! NOTES:
! ** LIBSCI-specific routine!  
!
!  Swapped arguments in ZDFFTM_STRIDES as a fix 
!  for the X1E.
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
  c_nl_i(:) = dr_eodr(:)*c_NL_i(:)
  !----------------------------------------------------

  !----------------------------------------------------
  ! FFT stuff
  !
  n_max_d = (3*n_max)/2+1
  n_fft   = 2*n_max_d+2
  !
  nx_fft = n_fft/2+1
  ny_fft = n_fft

  ! Switched from ZDFFTM_STRIDES to DZFFTM_STRIDES and order of ny_fft,nx_fft
  call DZFFTM_STRIDES(n_fft,ny_fft,nx_fft)

  ! Original code, don't use until PrgEnv 55 is default
  ! call ZDFFTM_STRIDES(n_fft,nx_fft,ny_fft)

  allocate(x_fft(0:nx_fft-1,6*n_x))
  allocate(y_fft(0:ny_fft-1,6*n_x))

  allocate(table_cs(100+2*n_fft))
  allocate(work_cs((2*n_fft+4)*6*n_x))

  CALL ZDFFTM(0, &
       n_fft, &
       6*n_x, &
       1.0, &
       x_fft, &
       nx_fft, &
       y_fft, &
       ny_fft, &
       table_cs, &
       work_cs, &
       0)

  allocate(table_sc(100+2*n_fft))
  allocate(work_sc((2*n_fft+4)*6*n_x))

  CALL DZFFTM(0, &
       n_fft, &
       6*n_x, &
       1.0, &
       y_fft, &
       ny_fft, &
       x_fft, &
       nx_fft, &
       table_sc, &
       work_sc, &
       0)

  deallocate(x_fft)
  deallocate(y_fft)
  !----------------------------------------------------

end subroutine make_nl
