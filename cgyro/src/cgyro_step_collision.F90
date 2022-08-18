!------------------------------------------------------------------------------
! cgyro_step_collision.F90
!
! PURPOSE:
!  Take an implicit collision step using the pre-computed collision 
!  matrix.  Effectively, we compute the new collisional cap_H: 
!
!                       H = h + ze/T G phi
!------------------------------------------------------------------------------

#ifndef _OPENACC
subroutine cgyro_step_collision_cpu

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is,ivp,it,j,k,nj_loc
  complex, dimension(nv) :: bvec,cvec
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

!$omp parallel do private(ic,ic_loc,it,iv,ivp,cvec,bvec,cvec_re,cvec_im,cval) firstprivate(collision_model)
  do ic=nc1,nc2

     ic_loc = ic-nc1+1
     ! Set-up the RHS: H = f + ze/T G phi

     it = it_c(ic)

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
           cval = cmat(iv,ivp,ic_loc)
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
     enddo

    do k=1,nproc
       do j=1,nj_loc
          fsendf(j,ic_loc,k) = bvec(j+(k-1)*nj_loc)
       enddo
    enddo

    if (collision_field_model == 1) then
       if (.not.(n == 0 .and. ae_flag == 1)) then
          ! cap_h_v not re-used else
          do iv=1,nv
             cap_h_v(ic_loc,iv) = bvec(iv)
          enddo
       endif
    endif
 enddo

  call timer_lib_out('coll')

  ! Compute the new phi
  if (collision_field_model == 1) then
     if (.not.(n == 0 .and. ae_flag == 1)) then
        call cgyro_field_v
     endif
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
     enddo
     ! this should be coll_mem timer , but not easy with OMP
     do ic=1,nc
        cap_h_c(ic,iv_loc) = cap_h_ct(iv_loc,ic)
     enddo
     is = is_v(iv)
     do ic=1,nc
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision_cpu

  ! else _OPENACC
#else

subroutine cgyro_step_collision_gpu

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: j,k
  integer :: is,ivp,nj_loc

  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im
  real :: cval

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_mem')
  call parallel_lib_rtrans_pack_gpu(cap_h_c)
  call timer_lib_out('coll_mem')
  call timer_lib_in('coll_comm')
  call parallel_lib_r_do_gpu(cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------


  call parallel_lib_nj_loc(nj_loc)

  ! Note that if GPU is absent, gpu_bibmem_flag will be reset to 0

  if (gpu_bigmem_flag == 1) then
     call timer_lib_in('coll')

!$acc parallel loop gang private(bvec,cvec) &
!$acc& present(cmat,cap_h_v,fsendf) 
     do ic=nc1,nc2

        ic_loc = ic-nc1+1

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
           if (.not.(n == 0 .and. ae_flag == 1)) then
!$acc loop vector
              do iv=1,nv
                 cap_h_v(ic_loc,iv) = bvec(iv)
              enddo
           endif
        endif
     enddo

     call timer_lib_out('coll')

  else

     call timer_lib_in('coll_mem')
!$acc update host(cap_h_v)
     call timer_lib_out('coll_mem')

     call timer_lib_in('coll')

!$omp parallel do private(ic_loc,iv,ivp,cvec,bvec,cvec_re,cvec_im,cval)
     do ic=nc1,nc2

        ic_loc = ic-nc1+1

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
              cval = cmat(iv,ivp,ic_loc)
              bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re,cval*cvec_im)
           enddo
        enddo

        ! Pack communication array while bvec still in cache
        do k=1,nproc
           do j=1,nj_loc
              fsendf(j,ic_loc,k) = bvec(j+(k-1)*nj_loc)
           enddo
        enddo

        if (collision_field_model == 1) then
           if (.not.(n == 0 .and. ae_flag == 1)) then
              ! cap_h_v not re-used else
              do iv=1,nv
                 cap_h_v(ic_loc,iv) = bvec(iv)
              enddo
           endif
        endif
     enddo

    call timer_lib_out('coll')

    call timer_lib_in('coll_mem')
!$acc update device(fsendf)
    if (collision_field_model == 1) then
       if (.not.(n == 0 .and. ae_flag == 1)) then
!$acc update device (cap_h_v)
       endif
    endif
    call timer_lib_out('coll_mem')

  endif

  ! Compute the new phi
  if (collision_field_model == 1) then
     if (.not.(n == 0 .and. ae_flag == 1)) then
        call cgyro_field_v_gpu
     endif
  endif

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do_gpu(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

!$acc parallel loop collapse(2) gang vector private(iv_loc,is) &
!$acc&         present(is_v,psi,cap_h_c,cap_h_ct,cap_h_c,jvec_c,field,z,temp,h_x) &
!$acc&         default(none)
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        cap_h_c(ic,iv_loc) = cap_h_ct(iv_loc,ic)
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision_gpu

  ! endif infdef _OPENACC
#endif

subroutine cgyro_step_collision

  use timer_lib
  use cgyro_globals

  implicit none

#ifdef _OPENACC
  call cgyro_step_collision_gpu
#else
  call cgyro_step_collision_cpu
#endif

  if (collision_field_model == 0 .or. (n == 0 .and. ae_flag == 1)) then
     call cgyro_field_c
  endif

end subroutine cgyro_step_collision

