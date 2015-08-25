!-----------------------------------------------------------------
! cgyro_hsym.f90
!
! PURPOSE:
!  Compute moment used in conservative advection scheme
!-----------------------------------------------------------------

subroutine cgyro_hsym

  use parallel_lib

  use cgyro_globals

  implicit none

  integer :: is,ie,ix
  complex, dimension(nc_loc) :: moment

  call parallel_lib_r(transpose(h_x),cap_h_v)
  cap_h_v_prime(:,:) = (0.0,0.0)

  do is=1,n_species
     do ie=1,n_energy
        moment(:) = 0.0
        do ix=1,n_xi
           moment(:) = moment(:)+&
                cap_h_v(:,iv_v(ie,ix,is))*0.5*w_xi(ix)*abs(xi(ix))
        enddo
        do ix=1,n_xi
           cap_h_v_prime(:,iv_v(ie,ix,is)) = moment
        enddo
     enddo
  enddo

  call parallel_lib_f(cap_h_v_prime,cap_h_ct)
  h_xs = transpose(cap_h_ct)

end subroutine cgyro_hsym
