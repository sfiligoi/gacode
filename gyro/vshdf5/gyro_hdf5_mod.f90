!--------------------------------------------------
! A module with common data to make communication 
! with the distributed arrays from write_hdf5.F90
! easier.
!------------------------------------------------
 module hdf5_mod
  real, dimension(:,:,:), allocatable :: alpha_phi
  real, dimension(:,:,:), allocatable :: alpha_phi_fine

  !---------------------------------------------------------
  ! These should go into gyro_globals.f90, but hard code them
  ! for now.
  !---------------------------------------------------------
  ! For HDF5 plotting, we can specify the number of toroidal
  ! angles used in the alpha expansion.
!  integer :: n_alpha_plot=20
!  integer :: n_alpha_fine=1
  !---------------------------------------------------------

 end module 
