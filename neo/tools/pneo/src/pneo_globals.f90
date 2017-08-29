module pneo_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: ntot
  integer, dimension(:), allocatable :: ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8
  real, dimension(:,:), allocatable :: indata_vec
  real, dimension(:,:), allocatable :: indata_tot
  real, dimension(:,:), allocatable :: outdata_vec
  real, dimension(:,:), allocatable :: outdata_tot


end module pneo_globals
