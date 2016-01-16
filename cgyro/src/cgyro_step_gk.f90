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
  real :: rval,rval2
  complex :: rhs_stream
  complex, dimension(:,:), allocatable :: g_x

  call timer_lib_in('str_comm')
  if (upconserve_flag == 1) call cgyro_hsym
  call timer_lib_out('str_comm')

  call timer_lib_in('str')
  rhs(ij,:,:) = (0.0,0.0)

  allocate(g_x(nc,nv_loc))
  g_x(:,:) = h_x(:,:) 

  ! Address cancellation problem
  if (n_field > 1) then
     do ic=1,nc
        ir = ir_c(ic) 
        it = it_c(ic)
        g_x(ic,:) = g_x(ic,:)+z(is)/temp(is)*j0_c(ic,:)*field(ir,it,2)*efac(:,2)
     enddo
  endif

!$omp parallel private(ic,iv_loc,is,ix,ie,ir,it,rval,rval2,rhs_stream,jt,jr,jc)
!$omp do
  do iv=nv1,nv2

     iv_loc = iv_locv(iv)
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        ! Diagonal terms
        rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc)+&
             omega_cap_h(ic,iv_loc)*cap_h_c(ic,iv_loc)+&
             omega_h(ic,iv_loc)*h_x(ic,iv_loc)+&
             sum(omega_s(:,ic,iv_loc)*field(ir,it,:))

        ! Parallel streaming with upwind dissipation 

        rval = omega_stream(it,is)*sqrt(energy(ie))*xi(ix) 
        if (upconserve_flag == 1) then
           rval2 = omega_stream(it,is)*sqrt(energy(ie)) 
        else
           rval2 = 0.0
        endif
        rhs_stream = 0.0

        if (implicit_flag == 0) then
           do id=-nup_theta,nup_theta
              jt = thcyc(it+id)
              jr = rcyc(ir,it,id)
              jc = ic_c(jr,jt)
              rhs_stream = rhs_stream &
                   -rval*dtheta(ir,it,id)*cap_h_c(jc,iv_loc)  &
                   -abs(rval)*dtheta_up(ir,it,id)*g_x(jc,iv_loc) &
                   +abs(rval2)*dtheta_up(ir,it,id)*h_xs(jc,iv_loc)
           enddo
        endif

        rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc)+rhs_stream

     enddo

  enddo
!$omp end do
!$omp end parallel

  deallocate(g_x)

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
