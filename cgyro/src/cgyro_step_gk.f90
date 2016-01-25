subroutine cgyro_step_gk

  use cgyro_globals

  implicit none

  ! RK4 time-advance for the distribution 
  !
  !           z e             vpar            z e  vperp^2
  !  h = H - ----- G0 ( phi - ----- Apar ) + ----- ---------- Gperp Bpar
  !            T               c               T   omega_a c
  !
  ! After time advance, we will have 
  !
  ! h    -> h_x
  ! H    -> cap_h_c
  ! phi  -> field(1)
  ! Apar -> field(2)
  ! Bpar -> field(3)

  h0_x = h_x
  
  ! Stage 1
  call cgyro_rhs(1)
  h_x = h0_x + 0.5 * delta_t * rhs(1,:,:)
  call cgyro_field_c

  ! Stage 2
  call cgyro_rhs(2)
  h_x = h0_x + 0.5 * delta_t * rhs(2,:,:)
  call cgyro_field_c

  ! Stage 3
  call cgyro_rhs(3)
  h_x = h0_x + delta_t * rhs(3,:,:)
  call cgyro_field_c

  ! Stage 4
  call cgyro_rhs(4)
  h_x = h0_x + delta_t/6.0 * &
       (rhs(1,:,:)+2.0*rhs(2,:,:)+2.0*rhs(3,:,:)+rhs(4,:,:))  
  call cgyro_field_c

  ! Filter special spectral components
  call cgyro_filter
  
end subroutine cgyro_step_gk
  
!==========================================================================

subroutine cgyro_rhs(ij)

  use timer_lib

  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: is, ir, it, ie, ix
  integer :: id, jt, jr, jc
  real :: rval
  complex :: rhs_stream
  complex :: rhs_ij(size(rhs,2),size(rhs,3))

  ! Prepare suitable distribution (g, not h) for conservative upwind method
  g_x(:,:) = h_x(:,:)
  if (n_field > 1) then
!$omp  parallel do collapse(2) &
!$omp& private(iv,ic,iv_loc,is,ir,it)
     do iv=nv1,nv2
     do ic=1,nc
        ! iv_loc = iv_locv(iv)
        iv_loc = iv-nv1+1
        is = is_v(iv)
        do ic=1,nc
           g_x(ic,iv_loc) = g_x(ic,iv_loc)+ & 
                z(is)/temp(is)*jvec_c(2,ic,iv_loc)*field(2,ic)
        enddo
     enddo
  endif

  call timer_lib_in('str_comm')
  call cgyro_hsym
  call timer_lib_out('str_comm')

  call timer_lib_in('str')

<<<<<<< HEAD
!$acc data  &
!$acc& pcopyout(rhs_ij) &
!$acc& pcopyin(h_x,field,cap_h_c) &
!$acc& present(is_v,ix_v,ie_v,ir_c,it_c) &
!$acc& present(omega_cap_h,omega_h,omega_s) &
!$acc& present(omega_stream,energy,xi) &
!$acc& present(thcyc,rcyc,ic_c,dtheta,dtheta_up)

!$acc kernels
!  rhs(ij,:,:) = (0.0,0.0)
   rhs_ij(:,:) = (0.0,0.0)
!$acc end kernels


#ifdef _OPENACC
!$acc parallel
!$acc loop collapse(2) gang vector &
!$acc& private(iv,ic,iv_loc,is,ix,ie,ir,it,rval,rhs_stream,id,jt,jr,jc)
#else
!$omp parallel private(ic,iv_loc,is,ix,ie,rval,rhs_stream,jc,id)
!$omp do 
#endif
  do iv=nv1,nv2
  do ic=1,nc

     ! iv_loc = iv_locv(iv)
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)


        ! Diagonal terms
        rhs_ij(ic,iv_loc) = rhs_ij(ic,iv_loc)+&
             omega_cap_h(ic,iv_loc)*cap_h_c(ic,iv_loc)+&
             omega_h(ic,iv_loc)*h_x(ic,iv_loc)+&
             sum(omega_s(:,ic,iv_loc)*field(:,ic))

        if (implicit_flag == 0) then
           ! Parallel streaming with upwind dissipation 
           rval = omega_stream(it_c(ic),is)*sqrt(energy(ie))
           rhs_stream = 0.0

           do id=-nup_theta,nup_theta
              jc = icd_c(ic,id)
              rhs_stream = rhs_stream &
                   -rval*xi(ix)*dtheta(ic,id)*cap_h_c(jc,iv_loc)  &
                   -abs(rval)*dtheta_up(ic,id)*g_x(jc,iv_loc) 
           enddo

           rhs_ij(ic,iv_loc) = rhs_ij(ic,iv_loc)+rhs_stream

        endif
     enddo
  enddo
#ifdef _OPENACC
!$acc end parallel

!$acc end data
#else
!$omp end do
!$omp end parallel
#endif
  rhs(ij,:,:) = rhs_ij(:,:)

  call timer_lib_out('str')

  ! Nonlinear evaluation [f,g]

  if (nonlinear_flag == 1) then     
     if (nonlinear_method == 1) then
        call cgyro_nl_direct(ij)
     else
        call cgyro_nl_fftw(ij)
     endif
  endif

end subroutine cgyro_rhs

!==========================================================================

subroutine cgyro_filter

  use cgyro_globals

  implicit none

  integer :: ir
  
  if (n == 0 .and. zf_test_flag == 0) then
     do ic=1,nc
        ir = ir_c(ic) 
        if (ir == 1 .or. px(ir) == 0) then
           h_x(ic,:)     = 0.0
           cap_h_c(ic,:) = 0.0
        endif
     enddo
  endif

end subroutine cgyro_filter
