module BLEND_private

  integer :: n
  integer :: n_fit

  real :: d

  complex :: phase

  complex, dimension(:,:), allocatable :: cs

  complex, dimension(:), allocatable :: c0

  !------------------------------------------
  ! For LAPACK:
  !
  integer, dimension(:), allocatable :: i_piv
  integer :: info
  !------------------------------------------

end module BLEND_private
