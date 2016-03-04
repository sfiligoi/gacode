!---------------------------------------------------------
! cgyro_shear_dft.f90
!
! PURPOSE:
!  DFT shear algorithm.  DFT coefficients are shifted
!  continuously to the left or right depending upon 
!  sign of omega_eb.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!---------------------------------------------------------

subroutine cgyro_shear_dft

  use cgyro_globals
  use timer_lib

  implicit none

  integer :: ix,p
  integer :: ir,it
  real :: a

  include 'fftw3.f03'

  a = omega_eb*delta_t

  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do it=1,n_theta

        do ir=1,n_radial
           ix = ir-1-n_radial/2
           if (ix < 0) ix = ix+n_radial
           fp(ix) = h_x(ic_c(ir,it),iv_loc)
        enddo

        call fftw_execute_dft(plan_p2j,fp,fj)
        do ir=1,n_radial
           p = ir-1-n_radial/2
           ix = p
           if (ix < 0) ix=ix+n_radial
           fj(ix) = exp(2*pi*i_c*p*a/n_radial)*fj(ix)
        enddo
        call fftw_execute_dft(plan_j2p,fj,fp)

        do ir=1,n_radial
           ix = ir-1-n_radial/2
           if (ix < 0) ix = ix+n_radial
           h_x(ic_c(ir,it),iv_loc) = fp(ix)/n_radial
        enddo

     enddo
  enddo

  call cgyro_field_c

end subroutine cgyro_shear_dft
