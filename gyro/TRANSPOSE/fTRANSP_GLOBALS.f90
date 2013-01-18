module fTRANSP_GLOBALS

  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: s0 
  integer :: s_dim

  integer :: i
  integer :: j
  integer :: k

  integer :: n_i
  integer :: n_j
  integer :: n_k
  integer :: n_m

  integer :: n_ij
  integer :: n_jk
  integer :: p_ij
  integer :: p_jk
  integer :: n_ij_loc
  integer :: n_jk_loc
  integer :: p_ij_loc
  integer :: p_jk_loc

  integer :: n_ijk

  integer :: i_recv
  integer :: i_from

  complex, dimension(:,:,:), allocatable :: q_send
  complex, dimension(:,:,:), allocatable :: q_recv

  integer, dimension(:), allocatable :: s

  integer, dimension(:,:), allocatable :: i_rc
  integer, dimension(:,:), allocatable :: jk
  integer, dimension(:,:), allocatable :: ij
  integer, dimension(:), allocatable :: i_ij,j_ij
  integer, dimension(:), allocatable :: j_jk,k_jk

  integer, dimension(:,:), allocatable :: i_map
  integer, dimension(:,:), allocatable :: p_jk_loc_map
  integer, dimension(:,:), allocatable :: i_map2
  integer, dimension(:,:), allocatable :: p_jk_loc_map2

  integer :: TRANSP_COMM

end module fTRANSP_GLOBALS
