module neo_sparse_solve

  implicit none

  public :: SOLVE_factor, SOLVE_do

  real,    dimension(10), private ::  cntl
  integer, dimension(20), private ::  icntl
  real,    dimension(20), private ::  rinfo
  integer, dimension(40), private ::  uinfo
  integer, dimension(20), private ::  keep
  
contains

  subroutine SOLVE_factor(n_elem)
    use neo_globals
    implicit none
    integer, intent (in) :: n_elem
    
    call UMD21I(keep,cntl,icntl)
 
    !--------------------------------------------
    ! Error and diagnostic messages from UMFPACK,
    ! to unit=5.
    !
    icntl(1) = 6
    icntl(2) = 6
    !
    ! Use icntl(3) = 2 for more verbosity.
    icntl(3) = 0
    !
    ! Use icntl(8)=n for n levels of iterative 
    ! refinement (set job=1).
    icntl(8) = 0
    !--------------------------------------------

    ! Factor Maxwell matrices
    call UMD2FA(n_row,&
         n_elem,&
         0,&
         .false.,&
         n_max,&
         2*n_max,&
         a,&
         a_indx,&
         keep,&
         cntl,&
         icntl,&
         uinfo,&
         rinfo)

    if(uinfo(1) < 0) then
       call neo_error('ERROR: Matrix factorization failed')
       return
    endif


  end subroutine SOLVE_factor

  subroutine SOLVE_do
    use neo_globals
    implicit none
    real, dimension(:), allocatable :: w, gs

    allocate(gs(n_row))
    allocate(w(2*n_row))

    call UMD2SO(n_row,&
         0,&
         .false.,&
         n_max,&
         2*n_max,&
         a,&
         a_indx,&
         keep,&
         g,&
         gs,&
         w,&
         cntl,&
         icntl,&
         uinfo,&
         rinfo)
    
    g = gs

    deallocate(gs)
    deallocate(w)

  end subroutine SOLVE_do

end module neo_sparse_solve
