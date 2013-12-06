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

  real, dimension(n_r), intent(in) :: s 
  real, dimension(n_r), intent(inout) :: p 

  ! Integrated power in erg/s

  p(1) = 0.0
  do i=2,n_r
     p(i) = p(i-1)+0.5*(s(i-1)*volp(i-1)+s(i)*volp(i))*(r(i)-r(i-1))
  enddo

end subroutine tgyro_volume_int
