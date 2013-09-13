subroutine tgyro_residual(f,g,res,n,method)

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: method
  real, intent(in), dimension(n) :: f
  real, intent(in), dimension(n) :: g
  real, intent(inout), dimension(n) :: res

  select case (method)

  case (1) 

     ! ORIGINAL METHOD:
     res = (f-g)**2/(f**2+g**2)

  case (2)

     ! SIMPLE NORM:
     res = abs(f-g)

  case (3)

     ! SQUARE RESIDUAL
     res = 0.5*(f-g)**2

  case (4)

     ! BALANCED
     res = (f-g)**2/MAX((f**2+g**2),1.0)

  case default

     print *,'Error in tgyro_residual'
     stop

  end select

end subroutine tgyro_residual
