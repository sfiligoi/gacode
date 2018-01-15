module pneo_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: ntot
  integer, dimension(:), allocatable :: ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8,ic9
  real, dimension(:,:), allocatable :: indata_loc
  real, dimension(:,:), allocatable :: indata
  real, dimension(:,:), allocatable :: ingeodata_loc
  real, dimension(:,:), allocatable :: ingeodata
  real, dimension(:,:), allocatable :: outdata_j_loc
  real, dimension(:,:), allocatable :: outdata_j
  real, dimension(:,:), allocatable :: outdata_ke_loc
  real, dimension(:,:), allocatable :: outdata_ke
  real, dimension(:,:), allocatable :: outdata_ki_loc
  real, dimension(:,:), allocatable :: outdata_ki
  real, dimension(:,:), allocatable :: outdata_kc_loc
  real, dimension(:,:), allocatable :: outdata_kc
  
  

end module pneo_globals
