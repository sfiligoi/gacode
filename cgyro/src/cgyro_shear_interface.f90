module cgyro_shear_interface

  real, dimension(:,:), allocatable :: fcoef0
  real, dimension(:,:), allocatable :: gcoef0
  complex, dimension(:,:,:), allocatable :: omega_s0
  real, dimension(:,:,:), allocatable :: jvec_c0
  complex, dimension(:,:), allocatable :: omega_cap_h0

end module cgyro_shear_interface
