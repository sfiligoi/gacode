subroutine cgyro_step_collision

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is,ir,it,ie,ix
  integer :: ivp

  ! compute new collisional cap_H: H = h + ze/T G phi
  ! assumes have cap_h_x

  call timer_lib_in('coll_comm')
  call parallel_lib_r(transpose(cap_h_c),cap_h_v)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

!$omp parallel private(ic,ic_loc,cvec,bvec,ivp,iv)
do ic=ibeg,iend
   ic_loc = ic_locv(ic)
  !ic_loc = 0
  !do ic=nc1,nc2
   !  ic_loc = ic_loc+1

     ! Set-up the RHS: H = f + ze/T G phi

     cvec(:) = cap_h_v(ic_loc,:)

     ! This is a key loop for performance
     bvec = (0.0,0.0)
     do ivp=1,nv
        do iv=1,nv
           bvec(iv) = bvec(iv)+cmat(iv,ivp,ic_loc)*cvec(ivp)
        enddo
     enddo

     cap_h_v(ic_loc,:) = bvec(:)

  enddo
!$omp end parallel

  call timer_lib_out('coll')

  ! Compute the new phi
  if (collision_field_model == 1) then
     call cgyro_field_v
  endif

  call timer_lib_in('coll_comm')
  call parallel_lib_f(cap_h_v,cap_h_ct)
  cap_h_c = transpose(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic)
        it = it_c(ic)

        psi(ic,iv_loc) = j0_c(ic,iv_loc)*efac(iv_loc,1)*field(ir,it,1)
        if (n_field > 1) then
           psi(ic,iv_loc) = psi(ic,iv_loc)+j0_c(ic,iv_loc)*efac(iv_loc,2)*field(ir,it,2)
           if (n_field > 2) then
              psi(ic,iv_loc) = psi(ic,iv_loc) &
                   +j0perp_c(ic,iv_loc)*efac(iv_loc,3)/Bmag(it)*field(ir,it,3)
           endif
        endif

        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*z(is)/temp(is)

     enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision
