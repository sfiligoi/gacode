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

  call parallel_lib_rtrans(g_x,cap_h_v)
!$omp workshare
  cap_h_v_prime(:,:) = (0.0,0.0)
!$omp end workshare


!$omp parallel do collapse(2) &
!$omp& private(is,ie,ix,tmp)
  do is=1,n_species
     do ie=1,n_energy
        tmp(:) = 0.0
        do ix=1,n_xi
           tmp(:) = tmp(:)+&
                cap_h_v(:,iv_v(ie,ix,is))*(w_xi(ix)*abs(xi(ix)))
        enddo
        do ix=1,n_xi
           cap_h_v_prime(:,iv_v(ie,ix,is)) = abs(xi(ix))*cap_h_v(:,iv_v(ie,ix,is))-tmp(:)
        enddo
     enddo
  enddo

  call parallel_lib_f(cap_h_v_prime,cap_h_ct)
  g_x = transpose(cap_h_ct)

end subroutine cgyro_hsym
