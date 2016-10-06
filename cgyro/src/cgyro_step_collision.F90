!------------------------------------------------------------------------------
! cgyro_step_collision.F90
!
! PURPOSE:
!  Take an implicit collision step using the pre-computed collision 
!  matrix.  Effectively, we compute the new collisional cap_H: 
!
!                       H = h + ze/T G phi
!------------------------------------------------------------------------------

subroutine cgyro_step_collision

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is,ie,ix
  integer :: ivp
  complex, dimension(size(cap_h_v,2),nc1:nc2) :: cvec
  complex, dimension(size(cap_h_v,2),nc1:nc2) :: bvec
  real :: cvec_re,cvec_im

  if (collision_model == 5) then
     call cgyro_step_collision_simple
     return
  endif

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_comm')
  call parallel_lib_rtrans(cap_h_c,cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------

  call timer_lib_in('coll')

#ifdef _OPENACC
!$acc  data present(cmat) &
!$acc& pcreate(bvec,cvec)  pcopy(cap_h_v)

!$acc  parallel 
!$acc  loop gang private(ic_loc,ivp,iv,cvec_re,cvec_im)
#else
!$omp parallel private(ic_loc,ivp,iv,cvec_re,cvec_im)
!$omp do
#endif
  do ic=nc1,nc2
     ic_loc = ic-nc1+1

     ! Set-up the RHS: H = f + ze/T G phi

!$acc loop vector
     do iv=1,nv
       cvec(iv,ic) = cap_h_v(ic_loc,iv)
     enddo

!$acc loop vector
     do iv=1,nv
       bvec(iv,ic) = (0.0,0.0)
     enddo

     ! This is a key loop for performance
     do ivp=1,nv
        cvec_re = real(cvec(ivp,ic),kind=kind(cmat))
        cvec_im = aimag(cvec(ivp,ic))
!$acc   loop vector
        do iv=1,nv
           bvec(iv,ic) = bvec(iv,ic)+ &
                 cmplx(cmat(iv,ivp,ic_loc)*cvec_re, &
                       cmat(iv,ivp,ic_loc)*cvec_im)
        enddo
     enddo

!$acc loop vector
     do iv=1,nv
       cap_h_v(ic_loc,iv) = bvec(iv,ic)
     enddo

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

!$omp parallel do &
!$omp& private(iv,iv_loc,is,ix,ie,ic)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  if (collision_field_model == 0) then
     call cgyro_field_c
  endif
  
  call timer_lib_out('coll')

end subroutine cgyro_step_collision
