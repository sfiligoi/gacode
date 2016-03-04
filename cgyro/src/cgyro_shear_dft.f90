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

  integer :: j,jx,it

  include 'fftw3.f03'

  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do it=1,n_theta
        
        fp(:) = h_x(ic_c(:,it),iv_loc)

        call fftw_execute_dft(plan_p2j,fp,fj)
        do j=1,n_radial
           jx = j-1
           if (jx > n_radial/2-1) jx=jx-n_radial
           fj(j) = exp(2*i_c*pi*jx*omega_eb*delta_t/n_radial)*fj(j)
        enddo
        call fftw_execute_dft(plan_j2p,fj,fp)

        h_x(ic_c(:,it),iv_loc) = fp(:)

     enddo
  enddo

  call cgyro_field_c

end subroutine cgyro_shear_dft
