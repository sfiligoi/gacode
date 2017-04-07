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

  integer :: is,ivp
  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_comm')
  call parallel_lib_rtrans(cap_h_c,cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------

  call timer_lib_in('coll')

!$omp parallel do private(iv,ivp,cvec,bvec,cvec_re,cvec_im)
  do ic_loc=1,nc2-nc1+1

     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(iv) = cap_h_v(ic_loc,iv)
     enddo

     bvec(:) = (0.0,0.0)

     ! This is a key loop for performance
     do ivp=1,nv
        cvec_re = real(cvec(ivp))
        cvec_im = aimag(cvec(ivp))
        do iv=1,nv
           bvec(iv) = bvec(iv)+ &
                cmplx(cmat(iv,ivp,ic_loc)*cvec_re, &
                cmat(iv,ivp,ic_loc)*cvec_im)
        enddo
     enddo

     do iv=1,nv
        cap_h_v(ic_loc,iv) = bvec(iv)
     enddo

  enddo

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

!$omp parallel do private(iv_loc,is,ic)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        chi(ic,iv_loc) = sum(jxvec_c(:,ic,iv_loc)*field(:,ic))
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  call timer_lib_out('coll')

  if (collision_field_model == 0) then
     call cgyro_field_c
  endif

end subroutine cgyro_step_collision
