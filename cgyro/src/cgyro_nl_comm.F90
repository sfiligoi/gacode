!-----------------------------------------------------------------
! cgyro_nl_comm.F90
!
! PURPOSE:
!  Nonlinear communication routines
!-----------------------------------------------------------------

module cgyro_nl_comm

  implicit none

contains

subroutine cgyro_nl_dealias_init

  use cgyro_globals

  implicit none

  integer :: itor,l0,l,p,m
  ! max_l0 = box_size*nt2*sign_qs
  integer :: pvec_count(box_size*nt2*sign_qs)

  ! Wavenumber M from CGYRO paper
  m = n_radial/2

  !
  ! max_pvec_count is global and we are defining it here
  !
  ! We will have arrays that depend on this later on
  ! so, precompute early
  !
  max_pvec_count = 1 ! make it >0, to always have valid arrays

  do itor=nt1,nt2
    if (itor/=0) then ! pvec not used else
      ! Total number of ballooning angles for finite-n ballooning mode
      l0 = box_size*itor*sign_qs

      pvec_count(1:l0) = 0

      ! Sort p indices by ballooning angle index l
      do p=-m,m-1
         l = mod(mod(p,l0)+l0,l0)+1
         pvec_count(l) = pvec_count(l)+1
         if (pvec_count(l) > max_pvec_count) then
           max_pvec_count = pvec_count(l)
         endif
      enddo
    endif
  enddo
end subroutine cgyro_nl_dealias_init


! === Internal ===
! Extended-angle filter for nonlinear dealiasing
! fraw (input unfiltered field)
! f    (output filtered field)

! itor=0 special case
! Pure periodicity 
pure recursive subroutine impfilter5_i0(&
                    n_theta,n_radial,a0,a1,a2,a3,fraw,f)
  implicit none

  integer, intent(in) :: n_theta,n_radial
  real, intent(in) :: a0,a1,a2,a3
  complex, intent(in) :: fraw(n_radial,n_theta)
  complex, intent(out):: f(n_radial,n_theta)
  !--------
  integer :: it, ir
  integer :: jm1,jm2,jm3,jp1,jp2,jp3

  do it=1,n_theta
     do ir=1,n_radial
        jm1 = it-1; if (jm1 < 1) jm1 = n_theta
        jm2 = it-2; if (jm2 < 1) jm2 = jm2+n_theta
        jm3 = it-3; if (jm3 < 1) jm3 = jm3+n_theta
        jp1 = it+1; if (jp1 > n_theta) jp1 = 1
        jp2 = it+2; if (jp2 > n_theta) jp2 = jp2-n_theta
        jp3 = it+3; if (jp3 > n_theta) jp3 = jp3-n_theta

        ! Pure periodicity 
        f(ir,it) = &
                a0*fraw(ir,it) + &
                a1*(fraw(ir,jm1) + fraw(ir,jp1)) + &
                a2*(fraw(ir,jm2) + fraw(ir,jp2)) + &
                a3*(fraw(ir,jm3) + fraw(ir,jp3))
     enddo
  enddo
end subroutine impfilter5_i0

! guaranteed that itor/=0
pure recursive subroutine impfilter5_n0(&
                    n_theta,n_radial,max_pvec_count,&
                    a0,a1,a2,a3,m,phase,l0,&
                    fraw,f,itor)

  implicit none

  integer, intent(in) :: n_theta,n_radial,max_pvec_count
  real, intent(in) :: a0,a1,a2,a3
  integer, intent(in) :: m
  complex, intent(in) :: phase
  integer, intent(in) :: l0
  complex, intent(in) :: fraw(n_radial,n_theta)
  complex, intent(out):: f(n_radial,n_theta)
  integer, intent(in) :: itor
  integer :: jm1,jm2,jm3,jp1,jp2,jp3
  integer :: ir,it,l,p,nex,iex,panel
  integer :: pvec(l0,max_pvec_count),pvec_count(l0)
  ! pre-allocate max possible size
  ! max_nex = n_theta*max_pvec_count
  complex :: fex(n_theta*max_pvec_count)
  integer :: ir_ex(n_theta*max_pvec_count)
  integer :: it_ex(n_theta*max_pvec_count)
  complex, parameter :: i_c  = (0.0,1.0)

  pvec_count(:) = 0

  ! Sort p indices by ballooning angle index l
  do p=-m,m-1
     l = mod(mod(p,l0)+l0,l0)+1
     pvec_count(l) = pvec_count(l)+1
     pvec(l,pvec_count(l)) = p
  enddo

  ! Construct ballooning modes
  ! all _ex quantities refer to extended angle
  ! nex = total number of points along extended angle
  ! iex = extended angle index
  do l=1,l0
     nex = n_theta*pvec_count(l)
     iex = 0
     do panel=1,pvec_count(l)
        p = pvec(l,panel)
        ir = p+m+1
        do it=1,n_theta
           iex = iex+1
           ir_ex(iex) = ir
           it_ex(iex) = it
           ! add phase 
           fex(iex) = fraw(ir,it)*exp(-i_c*p*phase)
        enddo
     enddo
     do iex=1,nex
        
        jm1 = iex-1; if (jm1 < 1) jm1 = nex
        jm2 = iex-2; if (jm2 < 1) jm2 = jm2+nex
        jm3 = iex-3; if (jm3 < 1) jm3 = jm3+nex
        jp1 = iex+1; if (jp1 > nex) jp1 = 1
        jp2 = iex+2; if (jp2 > nex) jp2 = jp2-nex
        jp3 = iex+3; if (jp3 > nex) jp3 = jp3-nex

        ir = ir_ex(iex) 
        it = it_ex(iex) 

        ! filter
        f(ir,it) = &
             a0*fex(iex) + &
             a1*(fex(jm1)+fex(jp1)) + &           
             a2*(fex(jm2)+fex(jp2)) + &           
             a3*(fex(jm3)+fex(jp3))            

        ! dephase
        p = ir-m-1
        f(ir,it) = f(ir,it)*exp(i_c*p*phase)

     enddo
  enddo

end subroutine impfilter5_n0

pure recursive subroutine impfilter5(&
                    n_theta,n_radial,max_pvec_count,&
                    dealias_order,dealias,&
                    box_size,q,sign_qs,&
                    fraw,f,itor)

  implicit none

  integer, intent(in) :: n_theta,n_radial,max_pvec_count
  integer, intent(in) :: dealias_order
  real, intent(in) :: dealias
  integer, intent(in) :: box_size
  real, intent(in) :: q
  integer, intent(in) :: sign_qs
  complex, intent(in) :: fraw(n_radial,n_theta)
  complex, intent(out):: f(n_radial,n_theta)
  integer, intent(in) :: itor
  real :: a0,a1,a2,a3
  integer :: m,l0
  complex :: phase
  real, parameter    :: pi   = 3.1415926535897932

  ! Maximally flat filters (3, 5, or 7-point):
  !
  ! H(0) = 1
  ! H'(0) = H''(0) = ... = 0
  ! H(pi) = 1-dealias
  !
  ! dealias = 0 (no filter)
  ! dealias = 1 (max filter)
  ! dealias_order = 1, 2, or 3 (for 3,5,7 pt)

  if (dealias_order == 1) then
     ! 3-point filter
     a0 = (2-dealias)/2.0
     a1 = dealias/4.0
     a2 = 0.0
     a3 = 0.0
  elseif (dealias_order == 2) then
     ! 5-point filter
     a0 = (8-3*dealias)/8.0
     a1 = dealias/4.0
     a2 = -dealias/16.0
     a3 = 0.0
  else
     ! 7-point filter (default)
     a0 = (16-5*dealias)/16.0
     a1 = (15/64.0)*dealias
     a2 = -(6/64.0)*dealias
     a3 = (1/64.0)*dealias
  endif

  if (itor == 0) then
     call impfilter5_i0(n_theta,n_radial,a0,a1,a2,a3,fraw,f)
  else
     ! Wavenumber M from CGYRO paper
     m = n_radial/2

     ! Phase factor
     phase = 2.0*pi*q/box_size

     ! Total number of ballooning angles for finite-n ballooning mode
     l0 = box_size*itor*sign_qs
     call impfilter5_n0(n_theta,n_radial,max_pvec_count,&
                        a0,a1,a2,a3,m,phase,l0,&
                        fraw,f,itor)
  endif

end subroutine impfilter5

subroutine hx_dealias_one(iv_loc_m,itor,hfil)
  use cgyro_globals

  implicit none
  
  integer, intent(in) :: iv_loc_m,itor
  complex, dimension(n_radial,n_theta), intent(out) :: hfil
  ! --------
  complex, dimension(n_radial,n_theta) :: hraw
  integer :: ir,it
  
  ! Construct h_x_dealias array
  do it=1,n_theta
    do ir=1,n_radial
      hraw(ir,it) = h_x(ic_c(ir,it),iv_loc_m,itor)
    enddo
  enddo
  ! Extended-angle dealiasing filter
  call impfilter5(n_theta,n_radial,max_pvec_count,&
                  dealias_order,dealias,&
                  box_size,q,sign_qs,&
                  hraw,hfil,itor)
end subroutine hx_dealias_one

subroutine field_dealias_one(i_field,itor,hfil)
  use cgyro_globals

  implicit none
  
  integer, intent(in) :: i_field,itor
  complex, dimension(n_radial,n_theta), intent(out) :: hfil
  ! --------
  complex, dimension(n_radial,n_theta) :: hraw
  integer :: ir,it
  
  ! Construct field_dealias array
  do it=1,n_theta
    do ir=1,n_radial
      hraw(ir,it) = field(i_field,ic_c(ir,it),itor)
    enddo
  enddo
  ! Extended angle dealiasing filter
  call impfilter5(n_theta,n_radial,max_pvec_count,&
                  dealias_order,dealias,&
                  box_size,q,sign_qs,&
                  hraw,hfil,itor)
end subroutine field_dealias_one

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm_test
  use parallel_lib
  use cgyro_globals

  implicit none

  if (fA_req_valid) call parallel_slib_test(fA_req)
  if (g_req_valid) call parallel_slib_test(g_req)
  if (fB_req_valid) call parallel_slib_test(fB_req)

end subroutine cgyro_nl_fftw_comm_test

!
! Comm is a transposea
! Reminder: nc ~= n_radial*n_theta
! First half of the transpose is done locally
!  from (theta,radial,nv_loc,nt_loc) -> (radial, nt_loc, theta, nv_loc)
! Then AlltoAll finishes the transpose
!  from (radial, nt_loc, theta, nv_loc_1, nv_loc_2) x toroidal_procs -> (radial, nt_loc, theta, nv_loc_1 , toroidal_procs) x nv_loc_2
! Implies nv_loc_2 == toroidal_procs
!

! NOTE: call cgyro_nl_fftw_comm1/2_async before cgyro_nl_fftw
subroutine cgyro_nl_fftw_comm1_f64_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex :: h_loc
  complex, dimension(n_radial,n_theta) :: hfil

  call timer_lib_in('nl_mem')

  if (nsplitB > 0) then

  ! TODO: GPU
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA,fpackB) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA,nsplitB) default(none)
#else
!$omp parallel do collapse(2) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc,hfil)
#endif
  do iv_loc_m=1,nv_loc
   do itor=nt1,nt2
    call hx_dealias_one(iv_loc_m,itor,hfil)
    do it=1,n_theta
     do ir=1,n_radial
       h_loc = hfil(ir,it)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
       else
          iexch_base = 1+itor0*nsplitB
          fpackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA)) = h_loc
       endif
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA,fpackB,nsplit,nsplitA,nsplitB)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
       else
          iexch_base = 1+itor0*nsplitB
          fpackB(1:n_radial,1:nt_loc,iexch_base+(isplit0-nsplitA)) = (0.0,0.0)
       endif
    enddo
  endif

  else ! nsplitB==0

  ! TODO: GPU
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA) default(none)
#else
!$omp parallel do collapse(2) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc,hfil)
#endif
  do iv_loc_m=1,nv_loc
   do itor=nt1,nt2
    call hx_dealias_one(iv_loc_m,itor,hfil)
    do it=1,n_theta
     do ir=1,n_radial
       h_loc = hfil(ir,it)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA,nsplit,nsplitA)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
    enddo
  endif

  endif ! if nspliB>0

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! split the comm in two, so we can start working on first as soon as it is ready
  call parallel_slib_f_nc_async(nsplitA,fpackA,fA_nl,fA_req)
  fA_req_valid = .TRUE.
  ! send only the first half, use comm3 to send the other half
  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm1_f64_async

subroutine cgyro_nl_fftw_comm1_f32_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex(KIND=REAL32) :: h_loc
  complex, dimension(n_radial,n_theta) :: hfil

  call timer_lib_in('nl_mem')

  if (nsplitB > 0) then

  ! TODO: GPU
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA32,fpackB32) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA,nsplitB) default(none)
#else
!$omp parallel do collapse(2) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc,hfil)
#endif
  do iv_loc_m=1,nv_loc
   do itor=nt1,nt2
    call hx_dealias_one(iv_loc_m,itor,hfil)
    do it=1,n_theta
     do ir=1,n_radial
       h_loc = hfil(ir,it)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA32(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
       else
          iexch_base = 1+itor0*nsplitB
          fpackB32(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA)) = h_loc
       endif
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA32,fpackB32,nsplit,nsplitA,nsplitB)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       if (isplit0 < nsplitA) then
          iexch_base = 1+itor0*nsplitA
          fpackA32(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
       else
          iexch_base = 1+itor0*nsplitB
          fpackB32(1:n_radial,1:nt_loc,iexch_base+(isplit0-nsplitA)) = (0.0,0.0)
       endif
    enddo
  endif

  else ! nsplitB==0

  ! TODO: GPU
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent &
!$acc&         private(iexch0,itor0,isplit0,iexch_base,h_loc) &
!$acc&         present(ic_c,h_x,fpackA32) &
!$acc&         present(n_theta,nv_loc,nt1,nt2,n_radial,nsplit,nsplitA) default(none)
#else
!$omp parallel do collapse(2) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base,h_loc,hfil)
#endif
  do iv_loc_m=1,nv_loc
   do itor=nt1,nt2
    call hx_dealias_one(iv_loc_m,itor,hfil)
    do it=1,n_theta
     do ir=1,n_radial
       h_loc = hfil(ir,it)
       iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA32(ir,itor-nt1+1,iexch_base+isplit0) = h_loc
     enddo
    enddo
   enddo
  enddo

  if ( (nv_loc*n_theta) < (nsplit*n_toroidal_procs) ) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(iexch0,itor0,isplit0,iexch_base)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) &
!$acc&         present(fpackA32,nsplit,nsplitA)
#endif
    do iexch0=nv_loc*n_theta,nsplit*n_toroidal_procs-1
       itor0 = iexch0/nsplit
       isplit0 = modulo(iexch0,nsplit)
       iexch_base = 1+itor0*nsplitA
       fpackA32(1:n_radial,1:nt_loc,iexch_base+isplit0) = (0.0,0.0)
    enddo
  endif

  endif ! if nspliB>0

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! split the comm in two, so we can start working on first as soon as it is ready
  call parallel_slib_f_nc32_async(nsplitA,fpackA32,fA_nl32,fA_req)
  fA_req_valid = .TRUE.
  ! send only the first half, use comm3 to send the other half
  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm1_f32_async

subroutine cgyro_nl_fftw_comm1_async

  use cgyro_globals

  implicit none
  !-----------------------------------

  if (nl_single_flag > 1) then
    call cgyro_nl_fftw_comm1_f32_async
  else
    call cgyro_nl_fftw_comm1_f64_async
  endif
end subroutine cgyro_nl_fftw_comm1_async

subroutine cgyro_nl_fftw_comm3_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  if (nsplitB > 0) then
    call timer_lib_in('nl_comm')

    if (nl_single_flag > 1) then
      call parallel_slib_f_nc32_async(nsplitB,fpackB32,fB_nl32,fB_req)
    else
      call parallel_slib_f_nc_async(nsplitB,fpackB,fB_nl,fB_req)
    endif
    fB_req_valid = .TRUE.

    call timer_lib_out('nl_comm')
  endif
end subroutine cgyro_nl_fftw_comm3_async

subroutine cgyro_nl_fftw_comm1_r64(ij)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  integer :: ir,it,iv_loc_m,ic_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex :: my_psi
  real :: psi_mul

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc_wait(nsplitA,fA_nl,fpackA,fA_req)
  fA_req_valid = .FALSE.
  if (nsplitB > 0) then
    ! no major compute to overlap
    call parallel_slib_r_nc_wait(nsplitB,fB_nl,fpackB,fB_req)
    fB_req_valid = .FALSE.
  endif
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  psi_mul = (q*rho/rmin)*(2*pi/length)

  if (nsplitB > 0) then

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi) firstprivate(px_zero)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero) &
!$acc&         present(ic_c,rhs,fpackA,fpackB) copyin(psi_mul,nl_min) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA,nsplitB) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  ((ir == 1) .or. (ir == px_zero)) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              if (isplit0 < nsplitA) then
                 iexch_base = 1+itor0*nsplitA
                 my_psi = fpackA(ir,itor-nt1+1,iexch_base+isplit0)
              else
                 iexch_base = 1+itor0*nsplitB
                 my_psi = fpackB(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))
              endif
           endif           
           if (itor < nl_min) then
              my_psi = (0.0,0.0)
           endif
           
           ! RHS -> -[f,g] = [f,g]_{r,-alpha}
           rhs(ic_loc_m,iv_loc_m,itor,ij) = rhs(ic_loc_m,iv_loc_m,itor,ij)+psi_mul*my_psi
        enddo
      enddo
    enddo
  enddo

  else ! nsplitB==0

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi) firstprivate(px_zero)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero) &
!$acc&         present(ic_c,rhs,fpackA) copyin(psi_mul,nl_min) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  ((ir == 1) .or. (ir == px_zero)) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              iexch_base = 1+itor0*nsplitA
              my_psi = fpackA(ir,itor-nt1+1,iexch_base+isplit0)
           endif           
           if (itor < nl_min) then
              my_psi = (0.0,0.0)
           endif
           
           ! RHS -> -[f,g] = [f,g]_{r,-alpha}
           rhs(ic_loc_m,iv_loc_m,itor,ij) = rhs(ic_loc_m,iv_loc_m,itor,ij)+psi_mul*my_psi
        enddo
      enddo
    enddo
  enddo

  endif ! if nsplitB>0

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw_comm1_r64


subroutine cgyro_nl_fftw_comm1_r32(ij)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  integer :: ir,it,iv_loc_m,ic_loc_m,itor
  integer :: iexch0,itor0,isplit0,iexch_base
  complex(KIND=REAL32) :: my_psi
  real(KIND=REAL32) :: psi_mul

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc32_wait(nsplitA,fA_nl32,fpackA32,fA_req)
  fA_req_valid = .FALSE.
  if (nsplitB > 0) then
    ! no major compute to overlap
    call parallel_slib_r_nc32_wait(nsplitB,fB_nl32,fpackB32,fB_req)
    fB_req_valid = .FALSE.
  endif
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  psi_mul = (q*rho/rmin)*(2*pi/length)

  if (nsplitB > 0) then

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi) firstprivate(px_zero)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero) &
!$acc&         present(ic_c,rhs,fpackA32,fpackB32) copyin(psi_mul,nl_min) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA,nsplitB) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  ((ir == 1) .or. (ir == px_zero)) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              if (isplit0 < nsplitA) then
                 iexch_base = 1+itor0*nsplitA
                 my_psi = fpackA32(ir,itor-nt1+1,iexch_base+isplit0)
              else
                 iexch_base = 1+itor0*nsplitB
                 my_psi = fpackB32(ir,itor-nt1+1,iexch_base+(isplit0-nsplitA))
              endif
           endif           
           if (itor < nl_min) then
              my_psi = (0.0,0.0)
           endif
           
           ! RHS -> -[f,g] = [f,g]_{r,-alpha}
           rhs(ic_loc_m,iv_loc_m,itor,ij) = rhs(ic_loc_m,iv_loc_m,itor,ij)+psi_mul*my_psi
        enddo
      enddo
    enddo
  enddo

  else ! nsplitB==0

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) &
!$omp&         private(ic_loc_m,my_psi) firstprivate(px_zero)
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang vector independent private(ic_loc_m,my_psi) &
!$acc&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero) &
!$acc&         present(ic_c,rhs,fpackA32) copyin(psi_mul,nl_min) &
!$acc&         present(nt1,nt2,nv_loc,n_theta,n_radial,nsplit,nsplitA) copyin(ij) default(none)
#else
!$omp parallel do collapse(2) private(ic_loc_m,my_psi) &
!$omp&         private(iexch0,itor0,isplit0,iexch_base) firstprivate(px_zero)
#endif
  do itor=nt1,nt2
    do iv_loc_m=1,nv_loc
      do it=1,n_theta
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (itor == 0) .and.  ((ir == 1) .or. (ir == px_zero)) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              iexch0 = (iv_loc_m-1) + (it-1)*nv_loc
              itor0 = iexch0/nsplit
              isplit0 = modulo(iexch0,nsplit)
              iexch_base = 1+itor0*nsplitA
              my_psi = fpackA32(ir,itor-nt1+1,iexch_base+isplit0)
           endif           
           if (itor < nl_min) then
              my_psi = (0.0,0.0)
           endif
           
           ! RHS -> -[f,g] = [f,g]_{r,-alpha}
           rhs(ic_loc_m,iv_loc_m,itor,ij) = rhs(ic_loc_m,iv_loc_m,itor,ij)+psi_mul*my_psi
        enddo
      enddo
    enddo
  enddo

  endif ! if nsplitB>0

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw_comm1_r32

subroutine cgyro_nl_fftw_comm1_r(ij)
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  if (nl_single_flag .EQ. 0) then
    call cgyro_nl_fftw_comm1_r64(ij)
  else
    call cgyro_nl_fftw_comm1_r32(ij)
  endif

end subroutine cgyro_nl_fftw_comm1_r

!
! Comm2 is a transpose
! Reminder: nc ~= n_radial*n_theta
! First half of the transpose is done locally with sub-sampling
!  from (n_field,n_theta,n_radial,nt_loc) -> (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! Then AlltoAll finishes the transpose
!  (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_proc)xn_toroidal_proc -> (n_field,n_radial,n_jtheta,nt_loc,n_toroida_procl)xn_toroidal_proc
! 

subroutine cgyro_nl_fftw_comm2_f64_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,it_loc,itm,itl,itf
  integer :: itor,mytor
  integer :: iltheta_min
  complex :: gval
  complex, dimension(n_radial,n_theta,n_field,nt1:nt2) :: field_dealias

  call timer_lib_in('nl_mem')

  ! TODO:GPU
  ! field-> gval processing is not uniform, create whole field_dealias
  !          It is small, anyway
!$omp parallel do collapse(2)
  do itor=nt1,nt2
    do itf=1,n_field
      call field_dealias_one(itf,itor,field_dealias(:,:,itf,itor))
    enddo
  enddo

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(5) &
!$omp&         private(itor,it,iltheta_min,mytor,gval)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(5) independent &
!$acc&         private(itor,it,iltheta_min,mytor,gval) &
!$acc&         present(field_dealias,gpack) &
!$acc&         present(n_toroidal_procs,nt_loc,n_jtheta,nv_loc,nt1) &
!$acc&         present(n_theta,n_radial,n_field,nsplit) &
!$acc&         default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(it_loc,itor,mytor,it,ir,iltheta_min,gval)
#endif
  do itm=1,n_toroidal_procs
   do itl=1,nt_loc
    do it_loc=1,n_jtheta
     do ir=1,n_radial
      do itf=1,n_field
       iltheta_min = 1+((itm-1)*nsplit)/nv_loc
       it = it_loc+iltheta_min-1
       itor = itl+(itm-1)*nt_loc
       gval = (0.0,0.0)
       if (it <= n_theta) then
         mytor = nt1+itl-1
         gval = field_dealias(ir,it,itf,mytor)
       endif
       ! else just padding
       gpack(itf,ir,it_loc,itor) = gval
      enddo
     enddo
    enddo
   enddo
  enddo

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd_async(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  g_req_valid = .TRUE.

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm2_f64_async

subroutine cgyro_nl_fftw_comm2_f32_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,it_loc,itm,itl,itf
  integer :: itor,mytor
  integer :: iltheta_min
  complex(KIND=REAL32) :: gval
  complex, dimension(n_radial,n_theta,n_field,nt1:nt2) :: field_dealias

  call timer_lib_in('nl_mem')

  ! TODO:GPU
  ! field-> gval processing is not uniform, create whole field_dealias
  !          It is small, anyway
!$omp parallel do collapse(2)
  do itor=nt1,nt2
    do itf=1,n_field
      call field_dealias_one(itf,itor,field_dealias(:,:,itf,itor))
    enddo
  enddo

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(5) &
!$omp&         private(itor,it,iltheta_min,mytor,gval)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(5) independent &
!$acc&         private(itor,it,iltheta_min,mytor,gval) &
!$acc&         present(field_dealias,gpack32) &
!$acc&         present(n_toroidal_procs,nt_loc,n_jtheta,nv_loc,nt1) &
!$acc&         present(n_theta,n_radial,n_field,nsplit) &
!$acc&         default(none)
#else
!$omp parallel do collapse(3) &
!$omp&         private(it_loc,itor,mytor,it,ir,iltheta_min,gval)
#endif
  do itm=1,n_toroidal_procs
   do itl=1,nt_loc
    do it_loc=1,n_jtheta
     do ir=1,n_radial
      do itf=1,n_field
       iltheta_min = 1+((itm-1)*nsplit)/nv_loc
       it = it_loc+iltheta_min-1
       itor = itl+(itm-1)*nt_loc
       gval = (0.0,0.0)
       if (it <= n_theta) then
         mytor = nt1+itl-1
         gval = field_dealias(ir,it,itf,mytor)
       endif
       ! else just padding
       gpack32(itf,ir,it_loc,itor) = gval
      enddo
     enddo
    enddo
   enddo
  enddo

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd32_async(n_field,n_radial,n_jtheta,gpack32,g_nl32,g_req)
  g_req_valid = .TRUE.

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm2_f32_async

subroutine cgyro_nl_fftw_comm2_async

  use cgyro_globals

  implicit none
  !-----------------------------------

  if (nl_single_flag > 1) then
    call cgyro_nl_fftw_comm2_f32_async
  else
    call cgyro_nl_fftw_comm2_f64_async
  endif
end subroutine cgyro_nl_fftw_comm2_async

end module cgyro_nl_comm

