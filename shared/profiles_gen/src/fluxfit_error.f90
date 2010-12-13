subroutine fluxfit_error(err)

  use prgen_fluxfit_globals

  implicit none

  ! Output
  real, intent(inout) :: err

  ! Number of points in search arc
  integer, parameter :: nt=16

  ! Iteration tolerance
  real, parameter :: tol=1e-9

  ! Internal variables
  integer :: i,p,k
  real, dimension(:), allocatable :: d
  real, dimension(:), allocatable :: tt
  real :: dmin,dmin0,q
  real :: r,z,t

  allocate(d(nt))
  allocate(tt(nt))

  err = 0.0
  do i=1,nd-1
     if (i < nd/2) then
        do k=1,nt 
           tt(k) = -0.5*pi+(k-1)*2.0*pi/(nt-1)
           call fluxfit_f_model(tt(k),r,z)
           d(k) = sqrt((r-rd(i))**2+(z-zd(i))**2)
        enddo
     else
        do k=1,nt 
           tt(k) = 0.5*pi+(k-1)*2.0*pi/(nt-1)
           call fluxfit_f_model(tt(k),r,z)
           d(k) = sqrt((r-rd(i))**2+(z-zd(i))**2)
        enddo
     endif

     call fluxfit_minmax(d,tt,nt,dmin,t,"min")

     dmin0 = 0.0
     p = 0

     do while (abs(dmin-dmin0) > tol)
        dmin0 = dmin
        p = p+1
        do k=1,nt 
           q = (nt*1.0)**p
           tt(k) = t-pi/q+(k-1.0)*2*pi/(nt-1)/q
           call fluxfit_f_model(tt(k),r,z)
           d(k) = sqrt((r-rd(i))**2+(z-zd(i))**2)
        enddo
        call fluxfit_minmax(d,tt,nt,dmin,t,"min")
     enddo
     err = err+dmin/(nd-1)/rmin
  enddo

  deallocate(d)
  deallocate(tt)

end subroutine fluxfit_error
