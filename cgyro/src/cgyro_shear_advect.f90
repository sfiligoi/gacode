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
  complex, dimension(0:n_radial+2,nv_loc) :: kp
  complex, dimension(n_radial,nv_loc) :: kpp

  wdt = omega_eb*delta_t

  do j=1,n_theta

     kp(:,:) = 0.0
     kp(1:n_radial,:) = h_x(ic_c(:,j),:)
     
     do ir=1,n_radial
        kpp(ir,:) = kp(ir,:)+wdt*(-kp(ir+2,:)+6*kp(ir+1,:)-3*kp(ir,:)-2*kp(ir-1,:))/6
     enddo

     h_x(ic_c(:,j),:) = kpp(1:n_radial,:)

  enddo

  call cgyro_field_c

end subroutine cgyro_shear_advect
