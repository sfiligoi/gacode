!-----------------------------------------------------
! make_ampere_matrix.f90
!
! PURPOSE:
!  Define sparse form of Ampere matrix (also used 
!  for A_par advance in the collision step) and 
!  factorize using UMFPACK.
!-----------------------------------------------------

subroutine make_ampere_matrix

  use gyro_globals
  use gyro_collision_private
  use math_constants

  !------------------------
  implicit none
  !-----------------------

  ! Row dimension of sparse Ampere matrix
  n_ampere_row = n_x*n_blend

  ! radial band width
  n_x2 = 2*mg_dx-ig_dx+1

  if (boundary_method == 1) then

     ! nonzero elements in n=0 Ampere matrix:
     n_zero = (n_x-1)*n_x2*n_blend**2+n_x*n_blend

     ! nonzero elements in n>0 Ampere matrix:
     n_fini = n_x*n_x2*n_blend**2

  else

     ! nonzero elements in nonperiodic Ampere matrices:
     n_zero = (n_x*n_x2-mg_dx*(mg_dx+1))*n_blend**2

     ! same number of nonzeros in n>0 nonperiodic matrices:
     n_fini = n_zero

  endif

  if (n_1(in_1) == 0) then
     n_ampere = n_zero
  else
     n_ampere = n_fini
  endif

  lindx(2)  = 2*n_ampere
  lvalue(2) = n_ampere

  allocate(m_ampere(lvalue(2)))
  allocate(indx_ampere(lindx(2)))

  if (n_1(in_1) == 0 .and. boundary_method == 1) then
     n_x_max  = n_x-1
  else
     n_x_max  = n_x
  endif

  allocate(aa_mm(n_x,-mg_dx:mg_dx-ig_dx,n_blend,n_blend))
  call make_ampere_blend

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_ampere,n_ampere_row,2,0)
  else
     call sparse_solve_mumps(n_ampere,n_ampere_row,2,0)
  endif

  !---------------------------------------------
  ! These are large matrices and deallocation is
  ! important:
  !
  deallocate(aa_mm)
  !---------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_ampere_matrix done]'
  endif

end subroutine make_ampere_matrix
