!-----------------------------------------------------
! make_maxwell_matrix.f90
!
! PURPOSE:
!  Define sparse form of Poisson-Ampere (Maxwell) 
!  matrix and factorize using UMFPACK.
! 
! REVISIONS
! 10 Sept 01: jc
!  Adapted from field_ampere.f90.
! 16 Jan 03: jc
!  Upgrade for IMEX scheme, renamed from 
!  field_ampere_sparse.
! 29 July 03: mf
!  Changed initial size of lvalue/lindx to be 
!  larger, make them grow faster.
!-----------------------------------------------------

subroutine make_maxwell_matrix

  use gyro_globals
  use maxwell_private

  !------------------------
  implicit none
  !------------------------

  !----------------------------------------------------------------
  ! Begin by doing a hand-count of the nonzero elements in 
  ! the sparse Maxwell matrix:
  !
  if (n_field == 1) then

     ! ELECTROSTATIC

     ! Row dimension of sparse Poisson matrix
     n_maxwell_row = n_x*n_blend

     ! radial band width
     n_gyro = 2*m_gyro-i_gyro+1

     ! theta/blending band width
     ! n_t2 = 2*blend_fit_order-1
     ! IMEX: no sparse assumption
     n_t2 = n_blend

     if (boundary_method == 1) then

        ! nonzero elements in n=0 Poisson matrix:
        n_zero = (n_x-1)*n_gyro*n_blend*n_t2+n_x*n_blend

        ! nonzero elements in n>0 Poisson matrix:
        n_fini = n_x*n_gyro*n_blend*n_t2

     else

        ! nonzero elements in n=0 nonperiodic Poisson matrices:
        n_zero = (n_x*n_gyro-m_gyro*(m_gyro+1))*n_blend*n_t2

        ! same number of nonzeros in n>0 nonperiodic matrices:
        n_fini = n_zero

     endif

  else if(n_field == 2) then
     
     ! ELECTROMAGNETIC -- A_parallel only

     ! Row dimension of sparse Maxwell submatrix
     n_maxwell_row = 2*n_x*n_blend

     n_dx = 2*mg_dx-ig_dx+1

     ! radial band width
     n_gyro = 2*m_gyro-i_gyro+1

     ! theta/blending band width
     ! n_t2 = 2*blend_fit_order-1
     ! IMEX: no sparse assumption
     n_t2 = n_blend

     if (boundary_method == 1) then

        ! nonzero elements in n=0 Ampere matrix:
        n_zero = (n_x-1)*n_gyro*n_blend*n_t2+n_x*n_blend & 
             + (n_x-1)*n_dx*n_blend*n_t2+n_x*n_blend & 
             + 2*(n_x-1)*n_blend*n_blend 

        ! nonzero elements in n>0 Ampere matrix:
        n_fini = n_x*n_gyro*n_blend*n_t2 + &
             n_x*n_dx*n_blend*n_t2 + & 
             2*n_x*n_blend*n_blend

     else

        ! nonzero elements in nonperiodic Ampere matrices:
        n_zero = (n_x*n_dx-mg_dx*(mg_dx+1))*n_blend*n_t2 & 
             + (n_x*n_gyro-m_gyro*(m_gyro+1))*n_blend*n_t2 &
             + 2*n_x*n_blend*n_blend

        n_fini = n_zero

     endif

  else

     ! ELECTROMAGNETIC -- A_parallel and B_parallel

     ! Row dimension of sparse Maxwell submatrix
     n_maxwell_row = 3*n_x*n_blend

     n_dx = 2*mg_dx-ig_dx+1

     ! radial band width
     n_gyro = 2*m_gyro-i_gyro+1

     ! theta/blending band width
     ! n_t2 = 2*blend_fit_order-1
     ! IMEX: no sparse assumption
     n_t2 = n_blend

     if (boundary_method == 1) then

        ! nonzero elements in n=0 Ampere matrix:
        n_zero = (n_x-1)*n_gyro*n_blend*n_t2+n_x*n_blend & 
             + (n_x-1)*n_dx*n_blend*n_t2+n_x*n_blend &
             + (n_x-1)*n_gyro*n_blend*n_t2+n_x*n_blend &
             + 2*(n_x-1)*n_gyro*n_blend*n_t2 &
             + 4*(n_x-1)*n_blend*n_blend 

        ! nonzero elements in n>0 Ampere matrix:
        n_fini = 4*n_x*n_gyro*n_blend*n_t2 + &
             n_x*n_dx*n_blend*n_t2 + & 
             4*n_x*n_blend*n_blend

     else

        ! nonzero elements in nonperiodic Ampere matrices:
        n_zero = (n_x*n_dx-mg_dx*(mg_dx+1))*n_blend*n_t2 & 
             + 4*(n_x*n_gyro-m_gyro*(m_gyro+1))*n_blend*n_t2 &
             + 4*n_x*n_blend*n_blend

        n_fini = n_zero

     endif


  endif
  !----------------------------------------------------------------

  if (n_1(in_1) == 0) then
     n_maxwell = n_zero
  else
     n_maxwell = n_fini
  endif

  lindx(3)  = 2*n_maxwell
  lvalue(3) = n_maxwell

  allocate(m_maxwell(lvalue(3)))
  allocate(indx_maxwell(lindx(3)))

  if (n_1(in_1) == 0 .and. boundary_method == 1) then
     n_x_max  = n_x-1
  else
     n_x_max  = n_x
  endif

  allocate(ap_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
  call make_poisson_blend(1)

  if (n_field > 1) then
     allocate(aa_mm(n_x,-mg_dx:mg_dx-ig_dx,n_blend,n_blend))
     call make_ampere_blend
  endif

  if (n_field == 3) then
     allocate(ab_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
     allocate(abp_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
     call make_ampereperp_blend
  endif

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_maxwell,n_maxwell_row,3,0)
  else
     call sparse_solve_mumps(n_maxwell,n_maxwell_row,3,0)
  endif

  !---------------------------------------------
  ! These are large matrices and deallocation is
  ! important:
  !
  deallocate(ap_mm)
  if (n_field > 1) deallocate(aa_mm)
  if (n_field == 3) then
     deallocate(ab_mm); deallocate(abp_mm)
  endif
  !---------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_maxwell_matrix done]'
  endif

end subroutine make_maxwell_matrix
