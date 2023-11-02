subroutine cgyro_rhs_comm_async
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
       call cgyro_nl_fftw_comm2_async
       call cgyro_nl_fftw_comm1_async
  endif

end subroutine cgyro_rhs_comm_async

subroutine cgyro_rhs_r_comm_async(ij)
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij

  if (nonlinear_flag == 1) then
       call cgyro_nl_fftw_comm1_r(ij)
  endif

end subroutine cgyro_rhs_r_comm_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_rhs_comm_test
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
     call cgyro_nl_fftw_comm_test
  endif

end subroutine cgyro_rhs_comm_test


subroutine cgyro_rhs(ij,update_cap)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  logical, intent(in) :: update_cap
  !--------------------------------
  integer :: is,itor
  integer :: id,jc
  real :: rval,rval2
  complex :: rhs_stream,h_el,cap_el,my_psi,rhs_el

  ! both h_x and field are ready by now
  call cgyro_rhs_comm_async

  call timer_lib_in('str')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(h_x,g_x,rhs,field,cap_h_c)
#endif


#if (!defined(OMPGPU)) && defined(_OPENACC)

!$acc data  &
!$acc& present(rhs) &
!$acc& present(cap_h_c,z,temp,jvec_c) &
!$acc& present(is_v,ix_v,ie_v,it_c) &
!$acc& present(omega_cap_h,omega_h,omega_s) &
!$acc& present(omega_stream,xi,vel) &
!$acc& present(dtheta,dtheta_up,icd_c)

#endif

  ! get the initial rhs initialization
  ! that depends on h_x,cap_h_c and field only
#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv,ic,iv_loc,rhs_el,h_el,cap_el,my_psi)
#elif defined(_OPENACC)
!$acc  parallel loop gang vector collapse(3) & 
!$acc& private(iv,ic,iv_loc,rhs_el,h_el,cap_el,my_psi) async(1)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        h_el = h_x(ic,iv_loc,itor)
        if (update_cap) then
           is = is_v(iv)
           my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
           cap_el = h_el+my_psi*z(is)/temp(is)
           cap_h_c(ic,iv_loc,itor) = cap_el
        else
           cap_el = cap_h_c(ic,iv_loc,itor)
        endif
        ! Diagonal terms
        rhs_el = &
             omega_cap_h(ic,iv_loc,itor)*cap_el+&
             omega_h(ic,iv_loc,itor)*h_el

        rhs(ic,iv_loc,itor,ij) = rhs_el + &
             sum(omega_s(:,ic,iv_loc,itor)*field(:,ic,itor))
     enddo
   enddo
  enddo

  call cgyro_rhs_comm_test

  ! Prepare suitable distribution (g, not h) for conservative upwind method
  if (n_field > 1) then

#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv_loc,is)
#elif defined(_OPENACC)
!$acc parallel loop  collapse(3) independent gang vector &
!$acc&         private(iv_loc,is) &
!$acc&         present(is_v,z,temp,jvec_c) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) &
!$acc&         default(none) async(1)
#endif
     do itor=nt1,nt2
      do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           is = is_v(iv)

           g_x(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+ & 
                (z(is)/temp(is))*jvec_c(2,ic,iv_loc,itor)*field(2,ic,itor)
        enddo
      enddo
     enddo

  else

#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop  collapse(3) gang vector &
!$acc&         independent private(iv_loc) &
!$acc&         present(nt1,nt2,nv1,nv2,nc)  &
!$acc&         default(none) async(1)
#endif
     do itor=nt1,nt2
      do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           g_x(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)
        enddo
      enddo
     enddo

  endif

#if (!defined(OMPGPU)) && defined(_OPENACC)
  call cgyro_rhs_comm_test
!$acc wait(1)
#endif

  call cgyro_rhs_comm_test

  call timer_lib_out('str')

  ! Wavenumber advection shear terms
  ! depends on h_x and field only, updates rhs
  call cgyro_advect_wavenumber(ij)

  call cgyro_rhs_comm_test

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then
     ! assumes someone already started the input comm
     ! and will finish the output comm
     call cgyro_nl_fftw()
  endif

  call cgyro_upwind

  call cgyro_rhs_comm_test

  call timer_lib_in('str')

  ! add stream to rhs
#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc) &
!$omp&  map(to:cap_h_c)
#elif defined(_OPENACC)
!$acc  parallel loop gang vector collapse(3) & 
!$acc& private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc) async(1)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ! Parallel streaming with upwind dissipation 
        rval  = omega_stream(it_c(ic),is,itor)*vel(ie_v(iv))*xi(ix_v(iv))
        rval2 = abs(omega_stream(it_c(ic),is,itor))

        rhs_stream = 0.0
        do id=-nup_theta,nup_theta
           jc = icd_c(id, ic, itor)
           rhs_stream = rhs_stream &
                -rval*dtheta(id,ic,itor)*cap_h_c(jc,iv_loc,itor)  &
                -rval2*dtheta_up(id,ic,itor)*g_x(jc,iv_loc,itor)
        enddo

        rhs(ic,iv_loc,itor,ij) = rhs(ic,iv_loc,itor,ij) + rhs_stream
     enddo
   enddo
  enddo

#if (!defined(OMPGPU)) && defined(_OPENACC)
  call cgyro_rhs_comm_test
!$acc wait(1)
#endif

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc end data    

! g_x and h_x
!$acc end data 
#endif

  call timer_lib_out('str')

  ! updates rhs
  call cgyro_rhs_r_comm_async(ij)

  if (explicit_trap_flag == 1) then
     ! we should never get in here... should have failed during init
     ! hard abort if we somehow end up here
     call abort
  endif

end subroutine cgyro_rhs
