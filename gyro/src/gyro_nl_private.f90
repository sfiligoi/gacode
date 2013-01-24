!-----------------------------------------------------
! gyro_nl_private.f90
!
! PURPOSE:
!  Module for variables used in do_nl_x.f90.  This 
!  routine containts control parameters for 
!   1. Direct method
!   2. FFTW FFT
!   3. FFTW3 FFT
!   4. ESSL FFT
!-------------------------------------------------------

module gyro_nl_private

  integer :: nn
  integer :: n1
  integer :: n2

  integer :: i_split

  real, dimension(:), allocatable :: c_nl_i

  integer, dimension(:), allocatable :: n_p
  integer, dimension(:), allocatable :: i_p  


  !-----------------------------------------------
  ! Generic FFT dimensions
  !
  integer :: n_max_d
  integer :: n_fft
  !-----------------------------------------------

  !-----------------------------------------------
  ! FFTW parameters:
  !
  integer(8) :: plan_f 
  integer(8) :: plan_b
  !
  ! FFTW2
  real, allocatable :: v_fft(:,:)
  real, allocatable :: vt_fft(:,:)
  !
  ! FFTW3
  real, allocatable :: v_fft3(:,:,:)
  real, allocatable :: vt_fft3(:,:,:)
  !-----------------------------------------------

  !-----------------------------------------------
  ! ESSL-specific parameters:
  !
  real, dimension(:), allocatable :: aux1_dcrft
  real, dimension(:), allocatable :: aux2_dcrft
  !
  real, dimension(:), allocatable :: aux1_drcft
  real, dimension(:), allocatable :: aux2_drcft
  !-----------------------------------------------

end module gyro_nl_private
