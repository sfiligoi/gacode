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
  complex, dimension(n_radial,nv_loc) :: kp
  complex, dimension(n_radial,nv_loc) :: kpp

  wdt = omega_eb*delta_t

  do j=1,n_theta

     kp(:,:) = h_x(ic_c(:,j),:)

     do ir=1,n_radial-2
        kpp(ir,:) = kp(ir,:)+wdt*0.5*(-kp(ir+2,:)+4*kp(ir+1,:)-3*kp(ir,:))
     enddo
     ir = n_radial-1
     kpp(ir,:) = kp(ir,:)+wdt*0.5*(4*kp(ir+1,:)-3*kp(ir,:))
     ir = n_radial
     kpp(ir,:) = 0.0

     h_x(ic_c(:,j),:) = kpp(:,:)

  enddo

  call cgyro_field_c

end subroutine cgyro_shear_advect
