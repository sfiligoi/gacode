!------------------------------------------------
! parallel_dim.f90
!
! PURPOSE:
!  Simple routine to get local array dimension
!  for array distributed along n_proc processors.
!
! REVISIONS:
! 26 Sept 01: jc 
!  Documented
!------------------------------------------------

integer function parallel_dim(n_grid,n_proc) 

  !-------------------------------
  implicit none
  !
  integer, intent(in) :: n_grid  
  integer, intent(in) :: n_proc
  integer :: i
  !-------------------------------

  if (n_grid < n_proc) then
     print '(a)','ERROR: (parallel_dim) More processors than gridpoints.'
     stop
  endif

  i = n_grid/n_proc
  if (modulo(n_grid,n_proc) > 0) then
     i = i+1
  endif

  parallel_dim = i

end function parallel_dim
