!--------------------------------------------------
! A module with common data to make communication 
! with the distributed arrays in some of the vshdf5
! easier.
!------------------------------------------------
 module gyro_vshdf5_mod
  real, dimension(:,:,:), allocatable :: alpha_phi
  real, dimension(:,:,:), allocatable :: alpha_phi_wedge

  !---------------------------------------------------------
  ! These should go into gyro_globals.f90, but hard code them
  ! for now.
  !---------------------------------------------------------
  ! For HDF5 plotting, we can specify the number of toroidal
  ! angles used in the alpha expansion.
!  integer :: n_torangle_3d=20
!  integer :: n_torangle_wedge=1
  !---------------------------------------------------------

 end module 
