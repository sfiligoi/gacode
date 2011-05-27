module peek_globals

  integer :: n_time
  integer :: i_spec
  integer :: i_field
  real, dimension(:), allocatable :: t
  real, dimension(:,:,:), allocatable :: gbflux_1
  real, dimension(:,:,:), allocatable :: gbflux_2
  real, dimension(:,:,:,:), allocatable :: gbflux_t
  character (len=10) :: tag(3)=(/'Phi ','Apar','Bpar'/)
  real :: window

end module peek_globals
