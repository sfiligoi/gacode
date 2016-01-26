subroutine cgyro_step_collision

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is,ie,ix
  integer :: ivp
  complex, dimension(size(cap_h_v,2),nc1:nc2) :: cvecm
  complex, dimension(size(cap_h_v,2),nc1:nc2) :: bvecm

  ! compute new collisional cap_H: H = h + ze/T G phi
  ! assumes have cap_h_x


  call timer_lib_in('coll_comm')
  call parallel_lib_rtrans(cap_h_c,cap_h_v)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

#ifdef _OPENACC
!$acc  data present(cmat) &
!$acc& pcreate(bvecm,cvecm)  pcopy(cap_h_v)

!$acc  parallel 
!$acc  loop gang private(ic_loc,ivp,iv)
#else
!$omp parallel private(ic_loc,ivp,iv)
!$omp do
#endif
  do ic=nc1,nc2
     ! ic_loc = ic_locv(ic)
     ic_loc = ic-nc1+1

     ! Set-up the RHS: H = f + ze/T G phi

     cvecm(:,ic) = cap_h_v(ic_loc,:)

     ! This is a key loop for performance
     bvecm(:,ic) = (0.0,0.0)
     do ivp=1,nv
!$acc   loop vector
        do iv=1,nv
           bvecm(iv,ic) = bvecm(iv,ic)+ &
                               cmat(iv,ivp,ic_loc)*cvecm(ivp,ic)
        enddo
     enddo

     cap_h_v(ic_loc,:) = bvecm(:,ic)

  enddo
#ifdef _OPENACC
!$acc end parallel
!$acc end data
#else
!$omp end do
!$omp end parallel
#endif

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
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*z(is)/temp(is)
     enddo
  enddo


  call timer_lib_out('coll')

end subroutine cgyro_step_collision
