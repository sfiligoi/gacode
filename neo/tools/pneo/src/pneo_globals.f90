module pneo_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: ntot
  integer, dimension(:), allocatable :: ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8,ic9
  real, dimension(:,:), allocatable :: indata_loc
  real, dimension(:,:), allocatable :: indata
  real, dimension(:,:), allocatable :: outdata_j_loc
  real, dimension(:,:), allocatable :: outdata_j
  real, dimension(:,:), allocatable :: outdata_u_loc
  real, dimension(:,:), allocatable :: outdata_u
  real, dimension(:,:), allocatable :: outdata_g_loc
  real, dimension(:,:), allocatable :: outdata_g
  real, dimension(:,:), allocatable :: outdata_q_loc
  real, dimension(:,:), allocatable :: outdata_q
  
end module pneo_globals
