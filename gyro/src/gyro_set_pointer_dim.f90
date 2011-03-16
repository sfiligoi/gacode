!-------------------------------------------------------------
! gyro_set_pointer_dim.f90
!
! PURPOSE:
!  Compute dimenions and pointers used for array distribution.
!-------------------------------------------------------------

subroutine gyro_set_pointer_dim

  use gyro_globals
  use gyro_pointers

  !--------------------------------
  implicit none
  !
  integer, external :: parallel_dim
  !--------------------------------

  ! Subgroup dimensions:

  n_nek_1     = n_n_1*n_energy*n_lambda
  n_nek_loc_1 = parallel_dim(n_nek_1,n_proc_1) 

  n_ine_1     = n_x*n_n_1*n_energy

  If (linsolve_method /= 2)   &
    n_ine_loc_1 = parallel_dim(n_ine_1,n_proc_1)

  n_ki_1     = n_lambda*n_x
  n_ki_loc_1 = parallel_dim(n_ki_1,n_proc_1) 

end subroutine gyro_set_pointer_dim
