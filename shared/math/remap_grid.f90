!-----------------------------------------------------
! remap_grid.f90
!
! PURPOSE:
!  Generate nonuniform grid based on a weight 
!  function.
!
! NOTES:
!  Given n_p "weight functions" f_p at points r_p 
! (possibly nonuniform), generate a nonuniform grid 
! r and uniform grid x (taken from r on input) 
! such that:
!
!        /      /
!    L = | dx = | dr = x(n)-x(1) = r(n)-r(1)
!        /      /
!
!  The grid spacings satisfy dx = f dr, and as 
!  a consequence:
!
!             /
!             | f dr = L
!             /
!
!  These formulae have the implication that 
!
!             dx = constant = f dr
! 
!  Also output is the array dxdr:
!
!                dxdr = dx / dr
!
! REVISIONS
! 11 Mar 02: jc
!  Created.
!---------------------------------------------

subroutine remap_grid(r_p,f_p,n_p,r,dxdr,x,n)

  !----------------------------------
  implicit none
  !
  integer, intent(in) :: n_p
  integer, intent(in) :: n
  !
  real, dimension(n_p), intent(in) :: r_p
  real, dimension(n_p), intent(in) :: f_p
  !
  real, dimension(n), intent(inout) :: r
  real, dimension(n), intent(inout) :: dxdr
  real, dimension(n), intent(inout) :: x
  !
  real, dimension(n_p) :: f_h
  !
  integer :: p
  integer :: i
  integer :: j
  !
  real :: sum0,sum1
  real :: length
  real :: dx
  real :: r0,r1
  real :: f0
  !
  !---------------------------------
  ! Integration intervals:
  !
  integer, parameter :: n_fine = 256
  !----------------------------------

  length = r(n)-r(1)
  dx = length/n_fine

  !                         /
  !  First, compute sum0 -> | f dr 
  !                         /
  
  sum0 = 0.0
  do p=1,n_fine

     r0 = r(1)+p*dx

     do j=1,n_p
        if (r_p(j) >= r0) exit
     enddo

     f0 = (f_p(j+1)*(r0-r_p(j))+f_p(j)*(r_p(j+1)-r0))/ & 
          (r_p(j+1)-r_p(j))

     sum0 = sum0+f0*dx

  enddo

  !                                 /
  ! Define an average f_h such that | dr f_h = 1
  !                                 /

  f_h(:) = length*f_p(:)/sum0

  x(:) = r(:)

  !                 r_i
  !                  /
  ! Solve (i-1) dx = | dr f_h(r)  for r_i. 
  !                  /  
  !                 r_1

  do i=1,n

     sum1 = x(1)-x(i)

     do p=1,n_fine

        r1 = x(1)+p*dx
        r0 = r1-dx

        do j=1,n_p
           if (r_p(j) > r1) exit
        enddo

        f0 = (f_h(j+1)*(r0-r_p(j))+f_h(j)*(r_p(j+1)-r0))/ & 
             (r_p(j+1)-r_p(j))

        sum0 = sum1
        sum1 = sum1+f0*dx

        if (sum1 > 0.0) then

           ! Root is bracketed inside (r0,r1)
           
           ! Get r(i) by linear interpolation 
           ! (i.e, secant method).

           r(i) = r0-sum0*(r1-r0)/(sum1-sum0)

           dxdr(i) = (f_h(j+1)*(r(i)-r_p(j))+f_h(j)*(r_p(j+1)-r(i)))/ & 
                (r_p(j+1)-r_p(j))

           goto 10

        endif

     enddo

10   continue

  enddo ! i

end subroutine remap_grid
