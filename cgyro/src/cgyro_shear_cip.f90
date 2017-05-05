!-------------------------------------------------------------------
! cgyro_shear_cip.f90
!
! PURPOSE:
!  CIP advection solve.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!--------------------------------------------------------------------

subroutine cgyro_shear_cip

  use cgyro_globals

  implicit none

  real :: z0,dx
  integer :: ir,j,iup
  complex, dimension(n_radial) :: h0,fp
  complex, dimension(0:n_radial+1) :: fc,gc
  complex :: a,b

  z0 = -omega_eb*delta_t
  dx = 1.0
  
  do iv_loc=1,nv_loc
     do j=1,n_theta

        h0(:) = h_x(ic_c(:,j),iv_loc)

        fc(0)          = h0(1)
        fc(1:n_radial) = h0(1:n_radial)
        fc(n_radial+1) = h0(n_radial)

        gc(:) = 0.0
        do ir=1,n_radial
           gc(ir) = (fc(ir+1)-fc(ir-1))/(2*dx)
        enddo

        do ir=1,n_radial
           iup = ir-1
           a   = (gc(ir)+gc(iup))/dx**2-2*(fc(ir)-fc(iup))/dx**3
           b   = 3*(fc(iup)-fc(ir))/dx**2+(2*gc(ir)+gc(iup))/dx
           fp(ir) = a*z0**3+b*z0**2+gc(ir)*z0+fc(ir)
        enddo

        h_x(ic_c(:,j),iv_loc) = fp(:)

     enddo
  enddo

  call cgyro_field_c

end subroutine cgyro_shear_cip
