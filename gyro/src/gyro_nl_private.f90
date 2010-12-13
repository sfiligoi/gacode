!-----------------------------------------------------
! gyro_nl_private.f90
!
! PURPOSE:
!  Module for variables used in do_nl_x.f90.  This 
!  routine containts control parameters for 
!   1. Direct method
!   2. FFTW FFT
!   3. ESSL FFT
!   4. LIBSC FFT
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
  ! Generic FFT parameters:
  integer :: n_max_d
  integer :: n_fft
  !-----------------------------------------------

  !-----------------------------------------------
  ! FFTW-specific parameters:
  !
  integer*8 :: plan_f 
  integer*8 :: plan_b
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

  !-----------------------------------------------
  ! LIBSCI-specific parameters:
  !
  integer :: nx_fft
  integer :: ny_fft
  !
  real, dimension(:), allocatable :: table_cs
  real, dimension(:), allocatable :: work_cs
  !
  real, dimension(:), allocatable :: table_sc
  real, dimension(:), allocatable :: work_sc
  !-----------------------------------------------

end module gyro_nl_private
