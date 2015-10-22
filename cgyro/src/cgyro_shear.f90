!---------------------------------------------------------
! cgyro_shear.f90
!
! PURPOSE:
!  Spectral shear algorithm.  Wavenumbers are shifted
!  to the left or right depending upon sign of omega_eb.
!---------------------------------------------------------

subroutine cgyro_shear

  use cgyro_globals
  use timer_lib

  implicit none

  integer :: ir
  complex, dimension(n_theta,nv_loc) :: a1

  
  gtime = gtime+omega_eb*delta_t

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
