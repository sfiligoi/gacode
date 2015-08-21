subroutine cgyro_hsym

  use parallel_lib

  use cgyro_globals

  implicit none

  integer :: is,ir,it,ie,ix,jx
  integer :: ixp
  
  call parallel_lib_r(transpose(h_x),cap_h_v)
  cap_h_v_prime(:,:) = (0.0,0.0)

  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1
     it = it_c(ic)
     ir = ir_c(ic)

     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        !cap_h_v_prime(ic_loc,iv) = cap_h_v(ic_loc,iv_v(ie,ix,is)) &
        !      +cap_h_v(ic_loc,iv_v(ie,n_xi-ix+1,is))
        do ixp=1,n_xi
           cap_h_v_prime(ic_loc,iv) = cap_h_v_prime(ic_loc,iv)+&
                cap_h_v(ic_loc,iv_v(ie,ixp,is))*0.5*w_xi(ixp)*abs(xi(ixp))
        enddo
     enddo
  enddo

  call parallel_lib_f(cap_h_v_prime,cap_h_ct)
  h_xs = transpose(cap_h_ct)

end subroutine cgyro_hsym
