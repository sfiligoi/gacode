module gyro_collision_private

  integer :: ij
  integer :: ijp
  integer :: k_amp
  integer :: n_x_max
  integer :: n_zero
  integer :: n_fini
  integer :: n_t2
  integer :: n_x2

  complex :: val

  ! collision_method:
  !  1 -> original RBF collisions with A_dot update
  !  2 -> new RBF scheme
  !  3 -> ebelli collisions (Connor)
  !  4 -> ebelli collisions (Connor-like with symmetric nu_d)

  !-------------------------------------------------
  ! Energy-dependent collision coefficients.  These
  ! go directly into the collision solve.
  !
  real, allocatable, dimension(:,:,:) ::   nu_total
  real, allocatable, dimension(:,:,:,:) :: nu_coll_d

  real, dimension(:,:,:), allocatable :: rs_coll_const
  real, dimension(:,:,:,:), allocatable :: rs_nunu_const
  integer, dimension(:,:), allocatable :: indx_coll

  !-------------------------------------------------

  real, allocatable, dimension(:,:,:) :: xi

end module gyro_collision_private
