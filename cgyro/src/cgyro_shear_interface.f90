module cgyro_shear_interface

  complex, dimension(:,:,:), allocatable :: omega_s0
  complex, dimension(:,:), allocatable :: omega_cap_h0

  real, dimension(:,:,:), allocatable :: jvec_c0
  real, dimension(:), allocatable :: k_perp0    

end module cgyro_shear_interface
