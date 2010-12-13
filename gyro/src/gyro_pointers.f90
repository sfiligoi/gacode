module gyro_pointers

  !----------------------------------------------------
  !  Pointer dimensions
  !
  ! n_nek = n_n*n_energy*n_lambda  
  !
  integer :: p_nek  
  !
  integer :: p_nek_loc
  !
  integer :: n_nek_1
  integer :: n_nek_loc_1
  integer :: n_nek_loc_true
  !
  ! n_ine = n_r*n_n*n_energy  
  !
  integer :: p_ine  
  !
  integer :: p_ine_loc
  !
  integer :: n_ine_1
  integer :: n_ine_loc_1
  !
  ! n_ki = n_lambda*n_x  
  !
  integer :: p_ki  
  !
  integer :: p_ki_loc
  !
  integer :: n_ki_1
  integer :: n_ki_loc_1
  !
  !----------------------------------------------------

  !----------------------------------------------------
  ! Pointers used for 1-D loop parallelization
  !
  integer, dimension(:), allocatable :: nek_n
  integer, dimension(:), allocatable :: nek_e
  integer, dimension(:), allocatable :: nek_k
  !
  integer, dimension(:), allocatable :: ine_i
  integer, dimension(:), allocatable :: ine_n
  integer, dimension(:), allocatable :: ine_e
  !
  integer, dimension(:), allocatable :: ki_k
  integer, dimension(:), allocatable :: ki_i
  !
  !----------------------------------------------------

end module gyro_pointers
