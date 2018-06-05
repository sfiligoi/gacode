!------------------------------------------------------------------------------
! cgyro_step_collision.gpu.F90
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

  integer :: is,ivp,it,j,k

#ifdef _OPENACC
  integer, parameter :: NVMAX = 758
  complex, dimension(NVMAX) :: bvec,cvec
#else
  complex, dimension(nv) :: bvec,cvec
#endif
  real :: cvec_re,cvec_im
  real :: cval

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_mem')
  call parallel_lib_rtrans_pack(cap_h_c)
  call timer_lib_out('coll_mem')
  call timer_lib_in('coll_comm')
  call parallel_lib_r_do(cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------

  call timer_lib_in('coll')

  call parallel_lib_nj_loc(nj_loc)

!$acc parallel loop gang private(bvec, cvec) &
!$acc& present(cmat) 
  do ic=nc1,nc2
!$acc cache(bvec,cvec)

     ic_loc = ic-nc1+1

     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

!$acc loop vector
     do iv=1,nv
        bvec(iv) = (0.0,0.0)
        cvec(iv) = cap_h_v(ic_loc,iv)
     enddo

     ! This is a key loop for performance
!$acc loop seq
     do ivp=1,nv
        cvec_re = real(cvec(ivp))
        cvec_im = aimag(cvec(ivp))
!$acc loop vector
        do iv=1,nv
           cval = cmat(iv,ivp,ic_loc)
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re,cval*cvec_im)
        enddo
     enddo

    ! pack communication array while bvec still in cache
    do k=1,nproc
!$acc loop vector
       do j=1,nj_loc
          fsendf(j,ic_loc,k) = bvec(j+(k-1)*nj_loc)
       enddo
    enddo

    if (collision_field_model == 1) then
       ! cap_h_v not re-used else
!$acc loop vector
       do iv=1,nv
          cap_h_v(ic_loc,iv) = bvec(iv)
       enddo
     endif
  enddo
!$acc end parallel

  call timer_lib_out('coll')

  ! Compute the new phi
  if (collision_field_model == 1) then
     call cgyro_field_v
  endif

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

!$omp parallel do private(iv_loc,is,ic,iv)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        chi(ic,iv_loc) = sum(jxvec_c(:,ic,iv_loc)*field(:,ic))
     enddo
     ! this should be coll_mem timer, but not easy with OMP
     do ic=1,nc
        cap_h_c(ic,iv_loc) = cap_h_ct(iv_loc,ic)
     enddo
     is = is_v(iv)
     do ic=1,nc
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  call timer_lib_out('coll')

  if (collision_field_model == 0) then
     call cgyro_field_c
  endif

end subroutine cgyro_step_collision
