!---------------------------------------------------------
! cgyro_shear_advect.f90
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!---------------------------------------------------------

subroutine cgyro_shear_advect

  use cgyro_globals

  implicit none

  integer :: j,ir
  real :: wdt
  complex, dimension(n_radial,n_theta) :: kp
  complex, dimension(n_radial,n_theta) :: kpp

  wdt = omega_eb*delta_t

  do j=1,n_theta

     kp(:,:) = h_x(ic_c(:,j),:)

     do ir=1,n_radial-1
        kpp(ir,:) = kp(ir,:)+wdt*(kp(ir+1,:)-kp(ir,:))
     enddo
     kpp(n_radial,:) = 0.0

     h_x(ic_c(:,j),:) = kpp(:,:)

  enddo

  call cgyro_field_c

end subroutine cgyro_shear_advect
