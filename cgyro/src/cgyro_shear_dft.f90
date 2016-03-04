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

  integer :: ir
  real :: dp
  complex, dimension(n_theta,nv_loc) :: a1


  dp = omega_eb*delta_t
  do i=1,n_radial
     p = i-n_radial/2-1
     do ip=1,n_radial
        m = ip-n_radial/2-1
        cr(i,ip) = exp(i_c*m*p*(pi_2/n_x))
     enddo
  enddo

  ! Forward shearing
  if (gtime > 0.5) then

     gtime = gtime-1.0

     a1 = h_x(ic_c(1,:),:)

     do ir=2,n_radial
        h_x(ic_c(ir-1,:),:) = h_x(ic_c(ir,:),:)
     enddo

     h_x(ic_c(n_radial,:),:) = a1*gamma_e_decay

     call cgyro_field_c

  endif

  ! Backward shearing
  if (gtime < -0.5) then

     gtime = gtime+1.0

     a1 = h_x(ic_c(n_radial,:),:) 

     do ir=n_radial-1,1,-1
        h_x(ic_c(ir+1,:),:) = h_x(ic_c(ir,:),:)
     enddo

     h_x(ic_c(1,:),:) = a1*gamma_e_decay

     call cgyro_field_c

  endif

end subroutine cgyro_shear
