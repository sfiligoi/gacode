module expromake_globals

  integer :: mode
  integer :: nx

  ! Input variables

  real, dimension(3) :: exm_z
  real :: exm_b_ref
  real :: exm_arho
  real :: exm_kappa
  real :: exm_delta
  real :: exm_te_axis
  real :: exm_alte
  real :: exm_ti_axis
  real :: exm_alti

  integer :: exm_set_b_ref
  integer :: exm_set_arho
  integer :: exm_set_kappa
  integer :: exm_set_delta
  integer :: exm_set_te_axis
  integer :: exm_set_alte
  integer :: exm_set_ti_axis
  integer :: exm_set_alti

end module expromake_globals
