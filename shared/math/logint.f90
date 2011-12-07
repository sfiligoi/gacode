!------------------------------------------------------
! logint.f90
!
! PURPOSE:
!  Given the (negative) logarithmic gradient f of a 
!  function fi:
!
!                     x*
!                     /
!  fi(x) = fi(x*)exp[ | f(y) dy ]
!                     /
!                     x
!  
!  obtain fi by trapezoidal integration.  It is assumed 
!  that fi(x*) is known.
!
! NOTES:
!      f(1:n) : (input) gradients
!       fi(n) : (input) boundary condition
!   fi(1:n-1) : (output) interior values
!------------------------------------------------------

subroutine logint(fi,f,r,n,i0)

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: i0

  real, intent(inout), dimension(n) :: fi
  real, intent(in), dimension(n) :: f
  real, intent(in),dimension(n) :: r

  integer :: i


  if (i0 > 1) then
     do i=i0,2,-1
        fi(i-1) = fi(i)*exp(0.5*(f(i)+f(i-1))*(r(i)-r(i-1)))
     enddo
  endif

  if (i0 < n) then
     do i=i0+1,n
        fi(i) = fi(i-1)*exp(-0.5*(f(i)+f(i-1))*(r(i)-r(i-1)))
     enddo
  endif

end subroutine logint
