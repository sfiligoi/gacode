!-----------------------------------------------------------------
! cgyro_zftest_em.f90
!
! PURPOSE:
!  initializes h and evaluates relevant integrals 
!  for an electromagnetic zonal flow test 
!-----------------------------------------------------------------

subroutine cgyro_zftest_em

  use mpi
  use cgyro_globals
  use cgyro_io
!  use cgyro_field_coefficients

  implicit none

  integer :: ir,it,is,ie,ix
  real :: arg, ang

  do iv=nv1,nv2

     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)       
     ie = ie_v(iv)
 
     do ic=1,nc

         ir = ir_c(ic) 
         it = it_c(ic)

         if (px(ir) /= 0) then
            arg = k_perp(ic)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
                 *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)           
            h_x(ic,iv_loc) = (z(is)/temp(is)) * 2 * (1-bessel_j0(abs(arg)))

         endif
     enddo
  enddo

end subroutine cgyro_zftest_em
