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

  integer :: ir,it

  call timer_lib_in('stream')

  gtime = gtime+omega_eb*delta_t

  ! Forward shearing
  if (gtime > 0.5) then

     gtime = gtime-1.0

     do ir=2,n_radial
        h_x(ic_c(ir-1,:),:)     = h_x(ic_c(ir,:),:)
        cap_h_c(ic_c(ir-1,:),:) = cap_h_c(ic_c(ir,:),:)
        psi(ic_c(ir-1,:),:)     = psi(ic_c(ir,:),:)
        field(ir-1,:,:)         = field(ir,:,:)
     enddo
     h_x(ic_c(n_radial,:),:)     = 0.0
     cap_h_c(ic_c(n_radial,:),:) = 0.0
     psi(ic_c(n_radial,:),:)     = 0.0
     field(n_radial,:,:)         = 0.0

  endif

  ! Backward shearing
  if (gtime < -0.5) then

     gtime = gtime+1.0

     do ir=n_radial-1,1,-1
        h_x(ic_c(ir+1,:),:)     = h_x(ic_c(ir,:),:)
        cap_h_c(ic_c(ir+1,:),:) = cap_h_c(ic_c(ir,:),:)
        psi(ic_c(ir+1,:),:)     = psi(ic_c(ir,:),:)
        field(ir+1,:,:)         = field(ir,:,:)
     enddo
     h_x(ic_c(1,:),:)     = 0.0
     cap_h_c(ic_c(1,:),:) = 0.0
     psi(ic_c(1,:),:)     = 0.0
     field(1,:,:)         = 0.0

  endif

  call timer_lib_out('stream')

end subroutine cgyro_shear
