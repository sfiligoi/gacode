module pneo_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: ntot
  integer :: ni,nj
  integer, dimension(:), allocatable :: ic,jc
  real, dimension(:,:), allocatable :: data_vec
  real, dimension(:,:), allocatable :: data_tot

end module pneo_globals
