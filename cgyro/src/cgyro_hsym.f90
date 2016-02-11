!-----------------------------------------------------------------
! cgyro_hsym.f90
!
! PURPOSE:
!  Compute moment used in conservative advection scheme:
!
!              /
!     |xi| g - | dxi |xi| g
!              /
!-----------------------------------------------------------------

subroutine cgyro_hsym

  use parallel_lib

  use cgyro_globals

  implicit none

  integer :: is,ie,ix
  complex, dimension(nc_loc) :: tmp
  real, dimension(nc_loc) :: jtmp

  call parallel_lib_rtrans(g_x,cap_h_v)
  cap_h_v_prime(:,:) = (0.0,0.0)

  do is=1,n_species
     tmp(:) = 0.0
     jtmp(:) = 0.0
     do ie=1,n_energy
        do ix=1,n_xi
           iv = iv_v(ie,ix,is)
           tmp(:) = tmp(:)+cap_h_v(:,iv)*w_xi(ix)*w_e(ie)*abs(xi(ix))*jvec_v(1,:,iv)*sqrt(energy(ie))
           jtmp(:) = jtmp(:)+w_xi(ix)*w_e(ie)*jvec_v(1,:,iv)**2
        enddo
     enddo
     do ie=1,n_energy
        do ix=1,n_xi
           iv = iv_v(ie,ix,is)
           cap_h_v_prime(:,iv) = sqrt(energy(ie))*abs(xi(ix))*cap_h_v(:,iv)-tmp(:)*jvec_v(1,:,iv)
        enddo
     enddo
  enddo

  call parallel_lib_f(cap_h_v_prime,cap_h_ct)
  g_x = transpose(cap_h_ct)

end subroutine cgyro_hsym
