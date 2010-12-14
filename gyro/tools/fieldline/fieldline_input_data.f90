module fieldline_input_data

  real, parameter :: pi=3.141592653589793238
  complex, parameter :: i_c=(0.0,1.0)

  integer :: n_eq=3

  ! GYRO profile data

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
  real, dimension(:), allocatable :: r_s
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: aspect_s
  real, dimension(:), allocatable :: shat_s
  real :: rhos_norm

  real, dimension(:,:), allocatable :: g_theta ! WG

  ! Derived data

  real :: r_length
  real :: y_length
  real :: r_natural
  real :: dr
  real :: dx
  real :: dtheta
  real :: lambda
  real :: eps_lambda

  integer, dimension(:), allocatable :: n
  complex, dimension(:,:), allocatable :: phase
  real, dimension(:), allocatable :: theta

  real, dimension(:), allocatable :: ky
  complex, dimension(:,:,:), allocatable :: a
  complex, dimension(:,:,:), allocatable :: anp

  real, dimension(:), allocatable :: rp
  real, dimension(:), allocatable :: xp
  real, dimension(:), allocatable :: qp

  real :: r0
  real :: q0
  real :: s0
  real :: rmaj0

  integer :: n_box
  integer :: n_x_decimate
  integer :: n_n_decimate
  integer :: n_turn
  integer :: n_sub
  integer :: n_ic
  integer :: n_time_slice

  integer :: ierr
 
end module fieldline_input_data
