!------------------------------------------------------------------------------
! cgyro_step_collision.F90
!
! PURPOSE:
!  Take an implicit collision step using the pre-computed collision 
!  matrix.  Effectively, we compute the new collisional cap_H: 
!
!                       H = h + ze/T G phi
!------------------------------------------------------------------------------

#if defined(_OPENACC) || defined(OMPGPU)
#define CGYRO_GPU_ROUTINES
#endif

subroutine cgyro_calc_collision_cpu_fp32(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !
  integer :: ivp,j,k,itor,ism
  integer :: vcount
  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im
  real :: cval
  ! --------------------------------------------------

  vcount = nv/nv_loc
!$omp parallel do collapse(3) &
!$omp&            private(ic,ic_loc,iv,ivp,cvec,bvec,cvec_re,cvec_im,cval,j,k) &
!$omp&            shared(cap_h_v,fsendf,cmat_fp32) firstprivate(vcount,nv,nj_loc)
  do itor=nt1,nt2
   do ic=nc_cl1,nc_cl2
    do ism=1,n_sim
     ic_loc = ic-nc_cl1+1
     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(iv) = cap_h_v(ic_loc,itor,iv,ism)
     enddo

     bvec(:) = (0.0,0.0)

     ! This is a key loop for performance
     do ivp=1,nv
        cvec_re = real(cvec(ivp))
        cvec_im = aimag(cvec(ivp))
        do iv=1,nv
           cval = cmat_fp32(iv,ivp,ic_loc,itor)
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
     enddo

    do k=1,vcount
       do j=1,nj_loc ! == nv_loc
          fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = bvec(j+(k-1)*nj_loc)
       enddo
    enddo

   enddo
  enddo
 enddo

end subroutine cgyro_calc_collision_cpu_fp32

subroutine cgyro_calc_collision_cpu_fp64(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !
  integer :: ivp,j,k,itor,ism
  integer :: vcount
  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im
  real :: cval
  ! --------------------------------------------------

  vcount = nv/nv_loc
!$omp parallel do collapse(3) &
!$omp&            private(ic,ic_loc,iv,ivp,cvec,bvec,cvec_re,cvec_im,cval,j,k) &
!$omp&            shared(cap_h_v,fsendf,cmat) firstprivate(vcount,nv,nj_loc)
  do itor=nt1,nt2
   do ic=nc_cl1,nc_cl2 ! == nc_loc_coll
    do ism=1,n_sim

     ic_loc = ic-nc_cl1+1
     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(iv) = cap_h_v(ic_loc,itor,iv,ism)
     enddo

     bvec(:) = (0.0,0.0)

     ! This is a key loop for performance
     do ivp=1,nv
        cvec_re = real(cvec(ivp))
        cvec_im = aimag(cvec(ivp))
        do iv=1,nv
           cval = cmat(iv,ivp,ic_loc,itor)
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
     enddo

    do k=1,vcount
       do j=1,nj_loc ! == nv_loc
          fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = bvec(j+(k-1)*nj_loc)
       enddo
    enddo

   enddo
  enddo
 enddo

end subroutine cgyro_calc_collision_cpu_fp64

subroutine cgyro_calc_collision_cpu_m1(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !
  integer :: ivp,j,k,itor,ism
  integer :: ie,is,ix,iep,isp,ixp
  integer :: vcount
  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im
  real :: cval
  ! --------------------------------------------------

  vcount = nv/nv_loc
!$omp parallel do collapse(3) &
!$omp&            private(ic,ic_loc,iv,ivp,cvec,bvec,cvec_re,cvec_im,cval,j,k) &
!$omp&            private(ie,is,ix,iep,isp,ixp) firstprivate(vcount,nv,nj_loc) &
!$omp&            shared(cap_h_v,fsendf,cmat_fp32,cmat_stripes,cmat_e1,ie_v,is_v,ix_v)
  do itor=nt1,nt2
   do ic=nc_cl1,nc_cl2
    do ism=1,n_sim

     ic_loc = ic-nc_cl1+1
     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(iv) = cap_h_v(ic_loc,itor,iv,ism)
     enddo

     bvec(:) = (0.0,0.0)

     ! This is a key loop for performance
     do ivp=1,nv
        cvec_re = real(cvec(ivp))
        cvec_im = aimag(cvec(ivp))
        iep = ie_v(ivp)
        isp = is_v(ivp)
        ixp = ix_v(ivp)
        do iv=1,nv
           cval = cmat_fp32(iv,ivp,ic_loc,itor)
           ie = ie_v(iv)
           is = is_v(iv)
           ix = ix_v(iv)
           if (ie<=n_low_energy) then
             cval = cval + cmat_e1(ix,is,ie,ivp,ic_loc,itor)
           else
             if ((iep==ie) .AND. (isp==is)) then
               cval = cval + cmat_stripes(ix,is,ie,ixp,ic_loc,itor)
             endif
           endif
           bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
        enddo
     enddo

    do k=1,vcount
       do j=1,nj_loc
          fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = bvec(j+(k-1)*nj_loc)
       enddo
    enddo

   enddo
  enddo
 enddo
end subroutine cgyro_calc_collision_cpu_m1

subroutine cgyro_calc_collision_cpu(nj_loc)

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  ! --------------------------------------------------

  if (collision_precision_mode == 1) then
     call cgyro_calc_collision_cpu_m1(nj_loc)
  else if (collision_precision_mode == 32) then
     call cgyro_calc_collision_cpu_fp32(nj_loc)
  else
     call cgyro_calc_collision_cpu_fp64(nj_loc)
  endif

end subroutine cgyro_calc_collision_cpu

#ifndef CGYRO_GPU_ROUTINES

subroutine cgyro_calc_collision_simple_cpu(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !

  integer :: is,ie,ix,jx,it,ir,j,k
  integer :: ivp,itor,ism
  integer :: vcount
  complex, dimension(:,:,:),allocatable :: bvec,cvec
  complex :: bvec_flat(nv)
  real :: cvec_re,cvec_im
  ! --------------------------------------------------

  allocate(bvec(n_xi,n_energy,n_species))
  allocate(cvec(n_xi,n_energy,n_species))

  vcount = nv/nv_loc
!$omp parallel do collapse(3) &
!$omp&            firstprivate(vcount,nv,nj_loc) &
!$omp&            private(ic_loc,ivp,iv,is,ix,jx,ie,ir,it,cvec_re,cvec_im,bvec,cvec,bvec_flat,k,j)
  do itor=nt1,nt2
   do ic=nc_cl1,nc_cl2
    do ism=1,n_sim
     ic_loc = ic-nc_cl1+1
     ir = ir_c(ic)
     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(ix_v(iv),ie_v(iv),is_v(iv)) = cap_h_v(ic_loc,itor,iv,ism)
     enddo

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. itor == 0) then
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
                         cmplx(cmat_simple(ix,jx,ie,is,it,itor)*cvec_re, &
                         cmat_simple(ix,jx,ie,is,it,itor)*cvec_im)
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
           fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = bvec_flat(j+(k-1)*nj_loc)
        enddo
     enddo

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

  integer :: is,nj_loc,itor
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
     call cgyro_calc_collision_cpu(nj_loc)
  endif

  call timer_lib_out('coll')

  if (.not. use_simple) then
    ! Compute the new phi
    if (collision_field_model == 1) then
        ! noop if (my_toroidal == 0 .and. ae_flag == 1))
        call cgyro_field_v_notae
    endif
  endif

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

!$omp parallel do collapse(2) &
!$omp&            private(iv_loc,is,ic,iv,my_psi,my_ch) firstprivate(nc)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        my_ch = cap_h_ct(iv_loc,itor,ic)
        my_psi = sum(jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        h_x(ic,iv_loc,itor) = my_ch-my_psi*(z(is)/temp(is))
        cap_h_c(ic,iv_loc,itor) = my_ch
     enddo
   enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision_cpu

  ! else CGYRO_GPU_ROUTINES
#else

  ! ==================================================

subroutine cgyro_calc_collision_gpu_fp32(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !

  integer :: j,k,ivp,itor,ism
  integer :: vcount
  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

  vcount = nv/nv_loc
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(b_re,b_im,cval,ivp,iv) firstprivate(nproc,nj_loc,nv,n_sim,vcount) &
!$omp&         private(k,ic,j,ic_loc,ism)
#else
!$acc parallel loop collapse(4) gang vector &
!$acc& private(b_re,b_im,cval,ivp,iv) firstprivate(nproc,nj_loc,nv,n_sim,vcount) &
!$acc& present(cmat_fp32,cap_h_v,fsendf)  private(k,ic,j,ic_loc,ism)
#endif
  do itor=nt1,nt2
    do ic=nc_cl1,nc_cl2  ! ==nc_loc_coll
      do k=1,vcount
        do j=1,nj_loc ! == nv_loc
           ic_loc = ic-nc_cl1+1
           iv = j+(k-1)*nj_loc
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(b_re,b_im,ivp)
#endif
           do ism=1,n_sim  ! keep ism as inner loop for cmat_fp32 locality
            b_re = 0.0
            b_im = 0.0
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(cval)
#endif
            do ivp=1,nv
              cval = cmat_fp32(iv,ivp,ic_loc,itor)
              b_re = b_re + cval*real(cap_h_v(ic_loc,itor,ivp,ism))
              b_im = b_im + cval*aimag(cap_h_v(ic_loc,itor,ivp,ism))
            enddo

            fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = cmplx(b_re,b_im)
           enddo
        enddo
      enddo
    enddo
  enddo
end subroutine cgyro_calc_collision_gpu_fp32

subroutine cgyro_calc_collision_gpu_fp64(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !

  integer :: j,k,ivp,itor,ism
  integer :: vcount
  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

  vcount = nv/nv_loc
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(b_re,b_im,cval,ivp,iv) firstprivate(nproc,nj_loc,nv,n_sim,vcount) &
!$omp&         private(k,ic,j,ic_loc,ism)
#else
!$acc parallel loop collapse(4) gang vector &
!$acc& private(b_re,b_im,cval,ivp,iv) firstprivate(nproc,nj_loc,nv,n_sim,vcount) &
!$acc& present(cmat,cap_h_v,fsendf)  private(k,ic,j,ic_loc,ism)
#endif
  do itor=nt1,nt2
    do ic=nc_cl1,nc_cl2 ! == nc_loc_coll
      do k=1,vcount
        do j=1,nj_loc ! == nv_loc
           ic_loc = ic-nc_cl1+1
           iv = j+(k-1)*nj_loc
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(b_re,b_im,ivp)
#endif
           do ism=1,n_sim  ! keep ism as inner loop for cmat locality
             b_re = 0.0
             b_im = 0.0
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(cval)
#endif
             do ivp=1,nv
              cval = cmat(iv,ivp,ic_loc,itor)
              b_re = b_re + cval*real(cap_h_v(ic_loc,itor,ivp,ism))
              b_im = b_im + cval*aimag(cap_h_v(ic_loc,itor,ivp,ism))
             enddo

             fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = cmplx(b_re,b_im)
           enddo
        enddo
      enddo
    enddo
  enddo
end subroutine cgyro_calc_collision_gpu_fp64

subroutine cgyro_calc_collision_gpu_m1(nj_loc)

  use parallel_lib
  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !

  integer :: j,k,ivp,itor,ism
  integer :: ie,is,ix,iep,isp,ixp
  integer :: vcount
  real :: b_re,b_im
  real :: h_re,h_im
  real :: cval
  ! --------------------------------------------------

  vcount = nv/nv_loc
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(b_re,b_im,cval,ivp,iv) firstprivate(nproc,nj_loc,nv,n_sim,vcount) &
!$omp&         private(k,ic,j,ic_loc,ie,is,ix,ism) &
!$omp&         private(iep,isp,ixp,h_re,h_im)
#else
!$acc parallel loop collapse(4) gang vector &
!$acc& private(b_re,b_im,cval,ivp,iv) firstprivate(nproc,nj_loc,nv,n_sim,vcount) &
!$acc& present(cmat_fp32,cmat_stripes,cmat_e1,cap_h_v,fsendf,ie_v,is_v,ix_v) &
!$acc& private(k,ic,j,ic_loc,ie,is,ix,ism) &
!$acc& private(iep,isp,ixp,h_re,h_im)
#endif
  do itor=nt1,nt2
   do ic=nc_cl1,nc_cl2 ! == nc_loc_coll
     do k=1,vcount
        do j=1,nj_loc ! == nv_loc
          ic_loc = ic-nc_cl1+1
          iv = j+(k-1)*nj_loc
          ie = ie_v(iv)
          is = is_v(iv)
          ix = ix_v(iv)
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(b_re,b_im,ivp)
#endif
          do ism=1,n_sim  ! keep ism as inner loop for cmat locality
            b_re = 0.0
            b_im = 0.0
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(cval,h_re,h_im,iep,isp,ixp)
#endif
            do ivp=1,nv
              cval = cmat_fp32(iv,ivp,ic_loc,itor)
              h_re = real(cap_h_v(ic_loc,itor,ivp,ism))
              h_im = aimag(cap_h_v(ic_loc,itor,ivp,ism))
              if (ie<=n_low_energy) then
                 cval = cval + cmat_e1(ix,is,ie,ivp,ic_loc,itor)
              else
                 iep = ie_v(ivp)
                 isp = is_v(ivp)
                 ixp = ix_v(ivp)
                 if ((ie==iep) .AND. (is==isp)) then
                    cval = cval + cmat_stripes(ix,is,ie,ixp,ic_loc,itor)
                 endif
              endif
              b_re = b_re + cval*h_re
              b_im = b_im + cval*h_im
            enddo

            fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = cmplx(b_re,b_im)
          enddo
        enddo
     enddo

   enddo
  enddo
end subroutine cgyro_calc_collision_gpu_m1

subroutine cgyro_calc_collision_gpu(nj_loc)

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  ! --------------------------------------------------

  if (collision_precision_mode == 1) then
     call cgyro_calc_collision_gpu_m1(nj_loc)
  else if (collision_precision_mode == 32) then
     call cgyro_calc_collision_gpu_fp32(nj_loc)
  else
     call cgyro_calc_collision_gpu_fp64(nj_loc)
  endif

end subroutine cgyro_calc_collision_gpu

subroutine cgyro_calc_collision_simple_gpu(nj_loc)

  use parallel_lib

  use cgyro_globals

  ! --------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nj_loc
  !

  integer :: is,ie,ix,jx,it,ir,j,k,jv
  integer :: ivp,itor, ism
  integer :: vcount

  real :: b_re,b_im
  real :: cval
  ! --------------------------------------------------

  vcount = nv/nv_loc
#if defined(OMPGPU)
!$omp target teams distribute collapse(2) &
!$omp&         private(ic_loc,ir,it,b_re,b_im,cval,ism) &
!$omp&         private(is,ie,ix,jx,iv,k,j,jv) firstprivate(vcount,nj_loc)
#else
!$acc parallel loop collapse(2) gang &
!$acc&         present(ix_v,ie_v,is_v,iv_v,ir_c,it_c,px,cap_h_v,cmat_simple,fsendf) &
!$acc&         private(ic_loc,ir,it,b_re,b_im,cval,ism) &
!$acc&         private(is,ie,ix,jx,iv,k,j,jv) firstprivate(vcount,nj_loc)
#endif
  do itor=nt1,nt2
   do ic=nc_cl1,nc_cl2 ! == nc_loc_coll

     ic_loc = ic-nc_cl1+1
     ir = ir_c(ic)
     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. itor == 0) then

        ! shortcut all the logic, just fill fsenf
#if defined(OMPGPU)
!$omp parallel do simd collapse(3) private(iv)
#else
!$acc loop collapse(3) vector private(iv)
#endif
        do k=1,vcount
          do ism=1,n_sim
           do j=1,nj_loc ! == nv_loc
              iv=j+(k-1)*nj_loc
              fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = cap_h_v(ic_loc,itor,iv,ism)
           enddo
         enddo
        enddo
     else
#if defined(OMPGPU)
!$omp parallel do simd collapse(3) private(iv,jv,is,ie,ix,jx,b_re,b_im,cval)
#else
!$acc loop collapse(3) vector private(iv,jv,is,ie,ix,jx,b_re,b_im,cval)
#endif
        do k=1,vcount
         do ism=1,n_sim
           do j=1,nj_loc ! == nv_loc
              iv = j+(k-1)*nj_loc
              is = is_v(iv)
              ie = ie_v(iv)
              ix = ix_v(iv)
              b_re = 0.0
              b_im = 0.0
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(cval,jv)
#endif
              do jx=1,n_xi
                 jv = iv_v(ie,jx,is)
                 cval = cmat_simple(ix,jx,ie,is,it,itor)
                 ! cvec_re(ix_v(jv),ie_v(jv),is_v(jv)) = real(cap_h_v(ic_loc,itor,jv))
                 ! b_re = b_re + cval*cvec_re(jx,ie,is)
                 b_re = b_re + cval*real(cap_h_v(ic_loc,itor,jv,ism))
                 b_im = b_im + cval*aimag(cap_h_v(ic_loc,itor,jv,ism))
              enddo
              fsendf(j,itor,ic_loc,k+(ism-1)*vcount) = cmplx(b_re,b_im)
           enddo
         enddo
        enddo
     endif

   enddo
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

  integer :: is,nj_loc,itor
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

    ! Note that if GPU is absent, gpu_bigmem_flag will be reset to 0
    if (gpu_bigmem_flag .GT. 0) then
      call timer_lib_in('coll')
      call cgyro_calc_collision_gpu(nj_loc)
      call timer_lib_out('coll')
    else

      call timer_lib_in('coll_mem')
#if defined(OMPGPU)
!$omp target update from(cap_h_v)
#else
!$acc update host(cap_h_v)
#endif
      call timer_lib_out('coll_mem')

      call timer_lib_in('coll')
      call cgyro_calc_collision_cpu(nj_loc)
      call timer_lib_out('coll')

      call timer_lib_in('coll_mem')
#if defined(OMPGPU)
!$omp target update to(fsendf)
#else
!$acc update device(fsendf)
#endif
      call timer_lib_out('coll_mem')

    endif !bigmem

    ! Compute the new phi
    if (collision_field_model == 1) then
        ! noop if (my_toroidal == 0 .and. ae_flag == 1))
        call cgyro_field_v_notae_gpu
    endif

  endif ! use_simple

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do_gpu(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&         private(iv_loc,is,my_psi,my_ch)
#else
!$acc parallel loop collapse(3) gang vector private(iv_loc,is,my_psi,my_ch) &
!$acc&         present(is_v,cap_h_c,cap_h_ct,cap_h_c,jvec_c,field,z,temp,h_x) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) default(none)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        my_psi = sum(jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        my_ch = cap_h_ct(iv_loc,itor,ic)
        h_x(ic,iv_loc,itor) = my_ch-my_psi*(z(is)/temp(is))
        cap_h_c(ic,iv_loc,itor) = my_ch
     enddo
   enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision_gpu

  ! endif infdef CGYRO_GPU_ROUTINES
#endif

  ! ==================================================

subroutine cgyro_step_collision

  use timer_lib
  use cgyro_globals

  implicit none

#ifdef CGYRO_GPU_ROUTINES
  call cgyro_step_collision_gpu(.FALSE.)
#else
  call cgyro_step_collision_cpu(.FALSE.)
#endif

  if (collision_field_model == 0) then
     call cgyro_field_c(.TRUE.)
  else if (nt1 == 0 .and. ae_flag == 1) then
     call cgyro_field_c_ae
  endif

end subroutine cgyro_step_collision

subroutine cgyro_step_collision_simple

  use timer_lib
  use cgyro_globals

  implicit none

#ifdef CGYRO_GPU_ROUTINES
  call cgyro_step_collision_gpu(.TRUE.)
#else
  call cgyro_step_collision_cpu(.TRUE.)
#endif

  call cgyro_field_c(.TRUE.)

end subroutine cgyro_step_collision_simple

