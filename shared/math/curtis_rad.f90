! Based on J. Candy's rotuine make_energy_grid_curtis
! Makes a radial grid based on the open
! Curtis-Clenshaw quadrature, also known as Fejer's second quadrature.
! This routine computes only the derivative weights 
! (i.e. the integration weights are not returned)

subroutine curtis_rad(n_radial,ra,rb,r,wr_deriv)

  !------------------------------------------
  implicit none
  !
  integer, intent(in) :: n_radial
  real, intent(in)    :: ra, rb
  real, intent(inout) :: r(n_radial)
  real, intent(inout) :: wr_deriv(n_radial,n_radial)
  !
  integer :: k, kp, n, nm
  real :: pi
  real :: xa, xb
  real :: thetak
  real :: thetakp
  real, dimension(:), allocatable :: x
  real, dimension(:,:), allocatable :: wd
  !------------------------------------------

  pi = 3.1415926535897932

  xa = 0.5*(rb-ra)
  xb = 0.5*(rb+ra)

  nm = n_radial-1
  allocate(x(0:nm))
  allocate(wd(0:nm,0:nm))

  ! Grid points
  do k=0,nm
     thetak = (2*k+1)*pi/(2.0*(nm+1))
     x(k)     = -cos(thetak)
  enddo

  ! Unscaled derivative matrix: wd(k,kp)
  do k=0,nm
     do kp=0,nm

        thetak  = (2*k+1)*pi/(2.0*(nm+1))
        thetakp = (2*kp+1)*pi/(2.0*(nm+1))

        ! First term is zero:
        wd(k,kp) = 0.0 
        do n=1,nm
           wd(k,kp) = wd(k,kp)-(2.0*n)/(nm+1)*cos(n*thetakp)* &
                sin(n*thetak)/sin(thetak)
        enddo

     enddo
  enddo

  ! Scaled quantities:

  r(:) = xa*x(:) + xb
  wr_deriv(:,:) = wd(:,:)/xa

end subroutine curtis_rad
