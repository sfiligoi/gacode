module cgyro_shear_interface

  real, dimension(:,:), allocatable :: fcoef0
  real, dimension(:,:), allocatable :: gcoef0
  real, dimension(:), allocatable :: sum_den_x0
  complex, dimension(:,:,:), allocatable :: omega_s0
  complex, dimension(:,:,:), allocatable :: jvec_c0
  complex, dimension(:,:), allocatable :: omega_cap_h0

end module cgyro_shear_interface
