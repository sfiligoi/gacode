module gyro2gkv_globals

  integer :: n_x
  integer :: n_theta_section
  integer :: n_pass
  integer :: n_trap
  integer :: n_energy
  integer :: n_theta_plot
  integer :: n0
  integer :: n_n
  integer :: d_n
  integer :: n_explicit_damp
  integer :: nonlinear_flag
  integer :: electron_method
  integer :: n_field
  integer :: n_ion
  integer :: n_kinetic
  integer :: n_spec
  integer :: field_r0_flag
  integer :: field_r0_grid
  integer :: n_grid_exp
  integer :: boundary_method

  real, dimension(:), allocatable :: r
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: theta_plot
  real, dimension(:), allocatable :: theta_r0_plot
  integer, dimension(:), allocatable :: n

  integer :: n_time
  integer :: i_time
  integer :: i_n
  integer :: exists_u
  integer :: exists_n
  integer :: exists_e
  integer :: exists_diff_i
  integer :: f_flag


  real :: pi
  real, dimension(:), allocatable :: t  

  real*8, dimension(:,:,:), allocatable :: diff0
  real, dimension(:,:,:), allocatable :: diff_density
  real, dimension(:,:,:), allocatable :: diff_energy
  real, dimension(:,:,:), allocatable :: diff_density_em
  real, dimension(:,:,:), allocatable :: diff_energy_em
  
  real*8, dimension(:,:,:,:), allocatable :: phi_in
  real, dimension(:,:,:,:), allocatable :: phi0

  real*8, dimension(:,:,:,:), allocatable :: phi_in_full
  real, dimension(:,:,:,:,:), allocatable :: phi0_full

  real*8, dimension(:,:,:,:), allocatable :: a_parallel_in
  real, dimension(:,:,:,:), allocatable :: a_parallel0

  real*8, dimension(:,:,:,:), allocatable :: a_parallel_in_full
  real, dimension(:,:,:,:,:), allocatable :: a_parallel0_full


  real*8, dimension(:,:,:), allocatable :: phi_r0_in
  real, dimension(:,:,:,:), allocatable :: phi_r0

  real*8, dimension(:,:,:,:), allocatable :: moment_n_in
  real, dimension(:,:,:,:), allocatable :: moment_n

  real*8, dimension(:,:,:,:), allocatable :: moment_n_full_in
  real, dimension(:,:,:,:,:), allocatable :: moment_n_full

  real*8, dimension(:,:,:,:), allocatable :: moment_e_in
  real, dimension(:,:,:,:), allocatable :: moment_e

  real*8, dimension(:,:,:,:), allocatable :: moment_e_full_in
  real, dimension(:,:,:,:,:), allocatable :: moment_e_full

  character (len=70) :: tag(5)
  integer :: lt(5)

end module gyro2gkv_globals
