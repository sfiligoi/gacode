module input_data

  integer :: n_r
  integer :: n_pass
  integer :: n_trap
  integer :: n_energy
  integer :: n_n
  integer :: n_field
  integer :: n_ion
  integer :: n_kinetic
  integer :: n_spec

  real :: rho 
  real, dimension(:), allocatable :: ky
  real, dimension(:), allocatable :: kx
  real, dimension(:), allocatable :: r
 
end module input_data
