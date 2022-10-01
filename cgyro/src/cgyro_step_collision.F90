!------------------------------------------------------------------------------
! cgyro_step_collision.F90
!
! PURPOSE:
!  Take an implicit collision step using the pre-computed collision 
!  matrix.  Effectively, we compute the new collisional cap_H: 
!
!                       H = h + ze/T G phi
!------------------------------------------------------------------------------

subroutine cgyro_calc_collision_cpu_fp64(nj_loc,update_chv)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  !
  integer :: ivp,j,k
  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im
  real :: cval
  ! --------------------------------------------------

!$omp parallel do private(ic,ic_loc,iv,ivp,cvec,bvec,cvec_re,cvec_im,cval,j,k) &
!$omp&            shared(cap_h_v,fsendf,cmat)
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
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
     enddo

    do k=1,nproc
       do j=1,nj_loc
          fsendf(j,ic_loc,k) = bvec(j+(k-1)*nj_loc)
       enddo
    enddo

    if (update_chv) then
          do iv=1,nv
             cap_h_v(ic_loc,iv) = bvec(iv)
          enddo
    endif
 enddo

end subroutine cgyro_calc_collision_cpu_fp64

subroutine cgyro_calc_collision_cpu_fp32(nj_loc,update_chv)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  !
  integer :: ivp,j,k,v1,v2
  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im
  real :: cval
  ! --------------------------------------------------

!$omp parallel do private(ic,ic_loc,iv,ivp,cvec,bvec,cvec_re,cvec_im,cval,j,k,v1,v2) &
!$omp&            shared(cap_h_v,fsendf,cmat_fp32,cmat_stripes)
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
        v1=max(1,ivp-collision_full_stripes)
        v2=min(ivp+collision_full_stripes,nv)
        do iv=1,v1-1
           cval = cmat_fp32(iv,ivp,ic_loc)
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
        do iv=v2+1,nv
           cval = cmat_fp32(iv,ivp,ic_loc)
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
        do iv=v1,v2
           cval = cmat_stripes(iv-ivp,ivp,ic_loc)
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
     enddo

    do k=1,nproc
       do j=1,nj_loc
          fsendf(j,ic_loc,k) = bvec(j+(k-1)*nj_loc)
       enddo
    enddo

    if (update_chv) then
          do iv=1,nv
             cap_h_v(ic_loc,iv) = bvec(iv)
          enddo
    endif
 enddo
end subroutine cgyro_calc_collision_cpu_fp32

subroutine cgyro_calc_collision_cpu(nj_loc,update_chv)

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  ! --------------------------------------------------

  if (collision_precision_mode == 0) then
     call cgyro_calc_collision_cpu_fp64(nj_loc,update_chv)
  else
     call cgyro_calc_collision_cpu_fp32(nj_loc,update_chv)
  endif

end subroutine cgyro_calc_collision_cpu

#ifndef _OPENACC

subroutine cgyro_calc_collision_simple_cpu(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !

  integer :: is,ie,ix,jx,it,ir,j,k
  integer :: ivp
  complex, dimension(:,:,:),allocatable :: bvec,cvec
  complex :: bvec_flat(nv)
  real :: cvec_re,cvec_im
  ! --------------------------------------------------

  allocate(bvec(n_xi,n_energy,n_species))
  allocate(cvec(n_xi,n_energy,n_species))

!$omp parallel do private(ic_loc,ivp,iv,is,ix,jx,ie,ir,it,cvec_re,cvec_im,bvec,cvec,bvec_flat,k,j)
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
     ir = ir_c(ic)
     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(ix_v(iv),ie_v(iv),is_v(iv)) = cap_h_v(ic_loc,iv)
     enddo

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. my_toroidal == 0) then
        bvec = cvec
     else

        bvec = 0.0

        do is=1,n_species
           do ie=1,n_energy              
              do jx=1,n_xi

                 cvec_re = real(cvec(jx,ie,is))
                 cvec_im = aimag(cvec(jx,ie,is))

                 do ix=1,n_xi
                    bvec(ix,ie,is) = bvec(ix,ie,is)+ &
                         cmplx(cmat_simple(ix,jx,ie,is,it)*cvec_re, &
                         cmat_simple(ix,jx,ie,is,it)*cvec_im)
                 enddo
              enddo
           enddo
        enddo
     endif

     do iv=1,nv
        bvec_flat(iv) = bvec(ix_v(iv),ie_v(iv),is_v(iv))
     enddo

     do k=1,nproc
        do j=1,nj_loc
           fsendf(j,ic_loc,k) = bvec_flat(j+(k-1)*nj_loc)
        enddo
     enddo

  enddo

  deallocate(bvec,cvec)
end subroutine cgyro_calc_collision_simple_cpu

  ! ==================================================

subroutine cgyro_step_collision_cpu(use_simple)

  use parallel_lib
  use timer_lib

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  logical, intent(in) :: use_simple
  !

  integer :: is,nj_loc
  logical :: update_chv
  complex :: my_psi,my_ch
  ! --------------------------------------------------

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
  if (use_simple) then
     call cgyro_calc_collision_simple_cpu(nj_loc)
  else
     update_chv = .FALSE.
     if (collision_field_model == 1) then
       if (.not.(my_toroidal == 0 .and. ae_flag == 1)) then
          ! cap_h_v not re-used else
          update_chv = .TRUE.
       endif
     endif

     call cgyro_calc_collision_cpu(nj_loc, update_chv)
  endif

  call timer_lib_out('coll')

  if (.not. use_simple) then
    ! Compute the new phi
    if (collision_field_model == 1) then
      if (.not.(my_toroidal == 0 .and. ae_flag == 1)) then
        call cgyro_field_v
      endif
    endif
  endif

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

!$omp parallel do private(iv_loc,is,ic,iv,my_psi,my_ch) firstprivate(nc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ! this should be coll_mem timer , but not easy with OMP
     do ic=1,nc
        my_ch = cap_h_ct(iv_loc,ic)
        my_psi = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        h_x(ic,iv_loc) = my_ch-my_psi*(z(is)/temp(is))
        cap_h_c(ic,iv_loc) = my_ch
     enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision_cpu

  ! else _OPENACC
#else

  ! ==================================================

subroutine cgyro_calc_collision_gpu_fp64(nj_loc,update_chv)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  !

  integer :: j,k,ivp
  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

!$acc parallel loop gang firstprivate(nproc,nj_loc,nv,update_chv) &
!$acc& present(cmat,cap_h_v,fsendf)  private(k,ic,j,ic_loc)
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
!$acc loop vector collapse(2) private(b_re,b_im,cval,ivp,iv)
     do k=1,nproc
        do j=1,nj_loc
           iv = j+(k-1)*nj_loc
           b_re = 0.0
           b_im = 0.0
!$acc loop seq private(cval)
           do ivp=1,nv
              cval = cmat(iv,ivp,ic_loc)
              b_re = b_re + cval*real(cap_h_v(ic_loc,ivp))
              b_im = b_im + cval*aimag(cap_h_v(ic_loc,ivp))
           enddo

           fsendf(j,ic_loc,k) = cmplx(b_re,b_im)
        enddo
     enddo

     if (update_chv) then
!$acc loop collapse(2) vector private(iv)
        do k=1,nproc
           do j=1,nj_loc
              iv = j+(k-1)*nj_loc
              cap_h_v(ic_loc,iv) = fsendf(j,ic_loc,k)
           enddo
        enddo
     endif
  enddo
end subroutine cgyro_calc_collision_gpu_fp64

subroutine cgyro_calc_collision_gpu_b2_fp64(nj_loc,update_chv)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  !

  integer :: bsplit
  integer :: j,k,ivp,b
  integer :: n_ic_loc,d_ic_loc
  integer, dimension((gpu_bigmem_flag*2)+1) :: bic
  integer :: bs,be,bb
  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

  bsplit = gpu_bigmem_flag*2

  n_ic_loc = nc2-nc1+1
  d_ic_loc = n_ic_loc/bsplit

  bic(1) = nc1
  do b=2,bsplit
     bic(b) = bic(b-1) + d_ic_loc
  enddo
  bic(bsplit+1) = nc2+1

  do b=1,bsplit
    bs = bic(b)-nc1+1
    be = bic(b+1)-nc1
    ! by keeping only 2 alive at any time, we limit GPU memory use
    bb = modulo(b,2)+2
    ! ensure there is not another even/odd already runnning
!$acc wait(bb)
    ! now launch myself
!$acc parallel loop gang firstprivate(nproc,nj_loc,nv,update_chv) &
!$acc& copyin(cmat(:,:,bs:be)) present(cap_h_v,fsendf)  private(k,j,ic_loc) async(bb)
    do ic_loc=bs,be
!$acc loop vector collapse(2) private(b_re,b_im,cval,ivp,iv)
      do k=1,nproc
        do j=1,nj_loc
           iv = j+(k-1)*nj_loc
           b_re = 0.0
           b_im = 0.0
!$acc loop seq private(cval)
           do ivp=1,nv
              cval = cmat(iv,ivp,ic_loc)
              b_re = b_re + cval*real(cap_h_v(ic_loc,ivp))
              b_im = b_im + cval*aimag(cap_h_v(ic_loc,ivp))
           enddo

           fsendf(j,ic_loc,k) = cmplx(b_re,b_im)
        enddo
      enddo

      if (update_chv) then
!$acc loop collapse(2) vector private(iv)
        do k=1,nproc
           do j=1,nj_loc
              iv = j+(k-1)*nj_loc
              cap_h_v(ic_loc,iv) = fsendf(j,ic_loc,k)
           enddo
        enddo
      endif
    enddo ! ic
  enddo ! b

  ! wait for all the async kernels to terminate
!$acc wait(2)
!$acc wait(3)

end subroutine cgyro_calc_collision_gpu_b2_fp64

subroutine cgyro_calc_collision_gpu_fp32(nj_loc,update_chv)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  !

  integer :: j,k,ivp,dv
  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

!$acc parallel loop gang firstprivate(nproc,nj_loc,nv,collision_full_stripes,update_chv) &
!$acc& present(cmat_fp32,cmat_stripes,cap_h_v,fsendf)  private(k,ic,j,ic_loc)
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
!$acc loop vector collapse(2) private(b_re,b_im,cval,ivp,iv)
     do k=1,nproc
        do j=1,nj_loc
           iv = j+(k-1)*nj_loc
           b_re = 0.0
           b_im = 0.0
!$acc loop seq private(cval,dv)
           do ivp=1,nv
              dv = iv-ivp
              if (abs(dv) .GT. collision_full_stripes) then
                 cval = cmat_fp32(iv,ivp,ic_loc)
              else
                 cval = cmat_stripes(dv,ivp,ic_loc)
              endif
              b_re = b_re + cval*real(cap_h_v(ic_loc,ivp))
              b_im = b_im + cval*aimag(cap_h_v(ic_loc,ivp))
           enddo

           fsendf(j,ic_loc,k) = cmplx(b_re,b_im)
        enddo
     enddo

     if (update_chv) then
!$acc loop collapse(2) vector private(iv)
        do k=1,nproc
           do j=1,nj_loc
              iv = j+(k-1)*nj_loc
              cap_h_v(ic_loc,iv) = fsendf(j,ic_loc,k)
           enddo
        enddo
     endif
  enddo
end subroutine cgyro_calc_collision_gpu_fp32

subroutine cgyro_calc_collision_gpu_b2_fp32(nj_loc,update_chv)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  !

  integer :: bsplit
  integer :: j,k,ivp,b,dv
  integer :: n_ic_loc,d_ic_loc
  integer, dimension((gpu_bigmem_flag*2)+1) :: bic
  integer :: bs,be,bb
  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

  bsplit = gpu_bigmem_flag*2

  n_ic_loc = nc2-nc1+1
  d_ic_loc = n_ic_loc/bsplit

  bic(1) = nc1
  do b=2,bsplit
     bic(b) = bic(b-1) + d_ic_loc
  enddo
  bic(bsplit+1) = nc2+1

  do b=1,bsplit
    bs = bic(b)-nc1+1
    be = bic(b+1)-nc1
    ! by keeping only 2 alive at any time, we limit GPU memory use
    bb = modulo(b,2)+2
    ! ensure there is not another even/odd already runnning
!$acc wait(bb)
    ! now launch myself
!$acc parallel loop gang firstprivate(nproc,nj_loc,nv,collision_full_stripes,update_chv) &
!$acc& copyin(cmat_fp32(:,:,bs:be),cmat_stripes(:,:,bs:be)) present(cap_h_v,fsendf)  private(k,j,ic_loc) async(bb)
    do ic_loc=bs,be
!$acc loop vector collapse(2) private(b_re,b_im,cval,ivp,iv)
       do k=1,nproc
         do j=1,nj_loc
           iv = j+(k-1)*nj_loc
           b_re = 0.0
           b_im = 0.0
!$acc loop seq private(cval,dv)
           do ivp=1,nv
              dv = iv-ivp
              if (abs(dv) .GT. collision_full_stripes) then
                 cval = cmat_fp32(iv,ivp,ic_loc)
              else
                 cval = cmat_stripes(dv,ivp,ic_loc)
              endif
              b_re = b_re + cval*real(cap_h_v(ic_loc,ivp))
              b_im = b_im + cval*aimag(cap_h_v(ic_loc,ivp))
           enddo

           fsendf(j,ic_loc,k) = cmplx(b_re,b_im)
         enddo
       enddo

       if (update_chv) then
!$acc loop collapse(2) vector private(iv)
         do k=1,nproc
           do j=1,nj_loc
              iv = j+(k-1)*nj_loc
              cap_h_v(ic_loc,iv) = fsendf(j,ic_loc,k)
           enddo
         enddo
      endif
    enddo ! ic
  enddo ! b

  ! wait for all the async kernels to terminate
!$acc wait(2)
!$acc wait(3)

end subroutine cgyro_calc_collision_gpu_b2_fp32

subroutine cgyro_calc_collision_gpu(nj_loc,update_chv)

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  ! --------------------------------------------------

  if (collision_precision_mode == 0) then
     call cgyro_calc_collision_gpu_fp64(nj_loc,update_chv)
  else
     call cgyro_calc_collision_gpu_fp32(nj_loc,update_chv)
  endif

end subroutine cgyro_calc_collision_gpu

subroutine cgyro_calc_collision_gpu_b2(nj_loc,update_chv)

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  logical, intent(in) :: update_chv
  ! --------------------------------------------------

  if (collision_precision_mode == 0) then
     call cgyro_calc_collision_gpu_b2_fp64(nj_loc,update_chv)
  else
     call cgyro_calc_collision_gpu_b2_fp32(nj_loc,update_chv)
  endif

end subroutine cgyro_calc_collision_gpu_b2

subroutine cgyro_calc_collision_simple_gpu(nj_loc)

  use parallel_lib

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !

  integer :: is,ie,ix,jx,it,ir,j,k
  integer :: ivp

  real, dimension(n_xi,n_energy,n_species) :: bvec_re,cvec_re
  real, dimension(n_xi,n_energy,n_species) :: bvec_im,cvec_im
  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

!$acc parallel loop gang &
!$acc&         present(ix_v,ie_v,is_v,ir_c,it_c,px,cap_h_v,cmat_simple,fsendf) &
!$acc&         private(bvec_re,bvec_im,cvec_re,cvec_im) &
!$acc&         private(ic_loc,ir,it,b_re,b_im,cval) &
!$acc&         private(is,ie,ix,jx,iv,k,j)
  do ic=nc1,nc2

     ic_loc = ic-nc1+1
     ir = ir_c(ic)
     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. my_toroidal == 0) then

        ! shortcut all the logic, just fill fsenf
!$acc loop collapse(2) vector private(iv)
        do k=1,nproc
           do j=1,nj_loc
              iv=j+(k-1)*nj_loc
              fsendf(j,ic_loc,k) = cap_h_v(ic_loc,iv)
           enddo
        enddo
     else
!$acc loop vector
        do iv=1,nv
           cvec_re(ix_v(iv),ie_v(iv),is_v(iv)) = real(cap_h_v(ic_loc,iv))
           cvec_im(ix_v(iv),ie_v(iv),is_v(iv)) = aimag(cap_h_v(ic_loc,iv))
        enddo

!$acc loop vector collapse(3) private(b_re,b_im,cval,jx)
        do is=1,n_species
           do ie=1,n_energy              
              do ix=1,n_xi
                 b_re = 0.0
                 b_im = 0.0
!$acc loop seq private(cval)
                 do jx=1,n_xi
                     cval = cmat_simple(ix,jx,ie,is,it)
                     b_re = b_re + cval*cvec_re(jx,ie,is)
                     b_im = b_im + cval*cvec_im(jx,ie,is)
                 enddo
                 bvec_re(ix,ie,is) = b_re
                 bvec_im(ix,ie,is) = b_im
              enddo
           enddo
        enddo

!$acc loop collapse(2) vector private(iv)
        do k=1,nproc
           do j=1,nj_loc
              iv=j+(k-1)*nj_loc
              fsendf(j,ic_loc,k) = cmplx(bvec_re(ix_v(iv),ie_v(iv),is_v(iv)),bvec_im(ix_v(iv),ie_v(iv),is_v(iv)))
           enddo
        enddo
     endif

  enddo

end subroutine cgyro_calc_collision_simple_gpu

  ! ==================================================

subroutine cgyro_step_collision_gpu(use_simple)

  use parallel_lib
  use timer_lib

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  logical, intent(in) :: use_simple
  !

  integer :: is,nj_loc
  logical :: update_chv
  complex :: my_psi,my_ch
  ! --------------------------------------------------

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

  if (use_simple) then
      call timer_lib_in('coll')
      call cgyro_calc_collision_simple_gpu(nj_loc)
      call timer_lib_out('coll')

  else

    update_chv = .FALSE.
    if (collision_field_model == 1) then
       if (.not.(my_toroidal == 0 .and. ae_flag == 1)) then
          ! cap_h_v not re-used else
          update_chv = .TRUE.
       endif
    endif

    ! Note that if GPU is absent, gpu_bibmem_flag will be reset to 0
    if (gpu_bigmem_flag == 1) then
      call timer_lib_in('coll')
      call cgyro_calc_collision_gpu(nj_loc, update_chv)
      call timer_lib_out('coll')

    elseif (gpu_bigmem_flag .GT. 1) then
      call timer_lib_in('coll')
      call cgyro_calc_collision_gpu_b2(nj_loc, update_chv)
      call timer_lib_out('coll')

    else

      call timer_lib_in('coll_mem')
!$acc update host(cap_h_v)
      call timer_lib_out('coll_mem')

      call timer_lib_in('coll')
      call cgyro_calc_collision_cpu(nj_loc, update_chv)
      call timer_lib_out('coll')

      call timer_lib_in('coll_mem')
!$acc update device(fsendf)
      if (collision_field_model == 1) then
        if (.not.(my_toroidal == 0 .and. ae_flag == 1)) then
!$acc update device (cap_h_v)
        endif
      endif
      call timer_lib_out('coll_mem')

    endif !bigmem

    ! Compute the new phi
    if (collision_field_model == 1) then
      if (.not.(my_toroidal == 0 .and. ae_flag == 1)) then
        call cgyro_field_v_gpu
      endif
    endif

  endif ! use_simple

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do_gpu(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

!$acc parallel loop collapse(2) gang vector private(iv_loc,is,my_psi,my_ch) &
!$acc&         present(is_v,cap_h_c,cap_h_ct,cap_h_c,jvec_c,field,z,temp,h_x) &
!$acc&         default(none)
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        my_psi = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        my_ch = cap_h_ct(iv_loc,ic)
        h_x(ic,iv_loc) = my_ch-my_psi*(z(is)/temp(is))
        cap_h_c(ic,iv_loc) = my_ch
     enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision_gpu

  ! endif infdef _OPENACC
#endif

  ! ==================================================

subroutine cgyro_step_collision

  use timer_lib
  use cgyro_globals

  implicit none

#ifdef _OPENACC
  call cgyro_step_collision_gpu(.FALSE.)
#else
  call cgyro_step_collision_cpu(.FALSE.)
#endif

  if (collision_field_model == 0 .or. (my_toroidal == 0 .and. ae_flag == 1)) then
     call cgyro_field_c
  endif

end subroutine cgyro_step_collision

subroutine cgyro_step_collision_simple

  use timer_lib
  use cgyro_globals

  implicit none

#ifdef _OPENACC
  call cgyro_step_collision_gpu(.TRUE.)
#else
  call cgyro_step_collision_cpu(.TRUE.)
#endif

  call cgyro_field_c

end subroutine cgyro_step_collision_simple

