subroutine cgyro_hsym

  use parallel_lib

  use cgyro_globals

  implicit none

  integer :: is,ir,it,ie,ix,jx
  integer :: icp,icm

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

        cap_h_v_prime(ic_loc,iv) = cap_h_v(ic_loc,iv_v(ie,ix,is)) &
             +cap_h_v(ic_loc,iv_v(ie,n_xi-ix+1,is))

     enddo
  enddo

  call parallel_lib_f(cap_h_v_prime,cap_h_ct)
  h_xs = transpose(cap_h_ct)*0.25

  cap_h_ct = h_xs

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ir=1,n_radial
        do it=1,n_theta

           icp = ic_c(ir,it)
           icm = ic_c(ir,n_theta-it+1)
           h_xs(ic,iv_loc) = cap_h_ct(icp,iv_loc)+cap_h_ct(icm,iv_loc)

        enddo
     enddo
  enddo

end subroutine cgyro_hsym
