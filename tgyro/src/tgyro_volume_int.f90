!---------------------------------------------------
! tgyro_volume_int.f90
!
! Integrate source to obtain integrated power:
!
! s    : source (input)
! volp : derivative of volume, V' (input)
! p    : power (output; erg/s if source in erg/cm^3/s)
!
!        r
!        /
! p(r) = | dx V'(x) s(x)
!        /
!        0
!
!---------------------------------------------------

subroutine tgyro_volume_int(s,p)

  use tgyro_globals, only : n_r,r,volp

  implicit none

  integer :: i

  real :: x1,x2,x3
  real :: f1,f2,f3

  real, dimension(n_r), intent(in) :: s 
  real, dimension(n_r), intent(inout) :: p 
 
  integer, parameter :: itrap=0

  ! Integrated power in erg/s

  if (n_r <= 2 .or. itrap == 1) then

     ! Integration using trapezoidal rule

     p(1) = 0.0
     do i=2,n_r
        p(i) = p(i-1)+0.5*(s(i-1)*volp(i-1)+s(i)*volp(i))*(r(i)-r(i-1))
     enddo

  else

     ! Integration using quadratic polynomial interpolation

     ! NOTE: this was introduced to significantly improve 
     !       accuracy for sparse (n_r < 10) grids.

     p(1) = 0.0
     do i=2,n_r

        if (i == 2) then

           ! Use points 1-3 for integral from 1-2

           x1 = r(i-1)
           x2 = r(i)
           x3 = r(i+1)

           f1 = s(i-1)*volp(i-1)
           f2 = s(i)*volp(i)
           f3 = s(i+1)*volp(i+1)

           p(i) = p(i-1)+(x2-x1)*(3*x3-x2-2*x1)*f1/6/(x3-x1) &
                + (x2-x1)*(3*x3-2*x2-x1)*f2/6/(x3-x2) &
                - (x2-x1)**3*f3/6/(x3-x1)/(x3-x2)
        else

           ! Use points 1-3 for integral from 2-3

           x1 = r(i-2)
           x2 = r(i-1)
           x3 = r(i)

           f1 = s(i-2)*volp(i-2) 
           f2 = s(i-1)*volp(i-1)
           f3 = s(i)*volp(i)

           p(i) = p(i-1)+(x3-x2)*(2*x3+x2-3*x1)*f3/6/(x3-x1) &
                + (x3-x2)*(x3+2*x2-3*x1)*f2/6/(x2-x1) &
                - (x3-x2)**3*f1/6/(x2-x1)/(x3-x1)
        endif

     enddo

  endif

end subroutine tgyro_volume_int

subroutine tgyro_volume_ave(f,r,volp,fave,n)

  implicit none

  integer :: i
  integer, intent(in) :: n
  real, intent(in), dimension(n) :: f,r,volp
  real, intent(inout) :: fave
  real :: vol

  fave = 0.0
  vol  = 0.0
  do i=2,n
     fave = fave+0.5*(f(i-1)*volp(i-1)+f(i)*volp(i))*(r(i)-r(i-1))
     vol  = vol+0.5*(volp(i-1)+volp(i))*(r(i)-r(i-1))
  enddo
  fave = fave/vol

end subroutine tgyro_volume_ave
