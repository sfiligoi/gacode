!--------------------------------------------------------------------
! tgyro_scaleint.f90
!
! PURPOSE:
!  Given the linear (z=f') or negative logarithmic (z=-f'/f) gradient 
!  of a function f, find the backward integral:
!
! type=log
!                    x*
!                    /
!  f(x) = f(x*) exp[ | z(s) ds ]
!                    /
!                    x
!
! type=lin
!                 x*
!                 /
!  f(x) = f(x*) - | z(s) ds
!                 /
!                 x
!
!
!  In other words, We obtain f by analytic integration of the piecewise
!  linear gradient function z(x).
!
! NOTES:
!        z(n) : (input) gradients
!        x(n) : (input) radii
!           n : (input) dimension of input arrays
!           f : (inout) f defined at n 
!---------------------------------------------------------------------

subroutine math_scaleintv(z,x,n,f,type)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: z(n),x(n)
  real, intent(inout), dimension(n) :: f
  character(len=3), intent(in) :: type

  integer :: i

  if (type == 'log') then
     do i=n,2,-1
        f(i-1) = f(i)*exp(0.5*(z(i)+z(i-1))*(x(i)-x(i-1)))
     enddo
  else
     do i=n,2,-1
        f(i-1) = f(i)-0.5*(z(i)+z(i-1))*(x(i)-x(i-1))
     enddo
  endif

end subroutine math_scaleintv
