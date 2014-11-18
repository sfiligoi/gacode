subroutine cgyro_step_gk

  use cgyro_globals

  implicit none

  ! RK4 time-advance for the distribution 
  !
  !           z e             vpar
  !  h = H - ----- G ( phi - ----- Apar )
  !            T               c
  !
  ! After time advance, we will have 
  !
  ! h    -> h_x
  ! H    -> cap_h_c
  ! phi  -> field(1)
  ! Apar -> field(2)

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

end subroutine cgyro_step_gk
  
!==========================================================================

subroutine cgyro_rhs(ij)

  use timer_lib

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer, intent(in) :: ij
  integer :: is, ir, it, ie, ix, jx
  integer :: id, jt, jr, jc
  real    :: rval

  call timer_lib_in('rhs')

  rhs(ij,:,:) = (0.0,0.0)

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        ! Parallel streaming with upwind dissipation
        rval = omega_stream(it,is)*sqrt(energy(ie))*xi(ix) 
        do id=-2,2
           jt = thcyc(it+id)
           jr = rcyc(ir,it+id)
           jc = ic_c(jr,jt)
           rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
                -rval*dtheta(it,id)*cap_h_c(jc,iv_loc)  &
                -abs(rval)*dtheta_up(it,id)*h_x(jc,iv_loc)
        enddo

        ! Diagonal terms
        rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc)+&
             omega_cap_h(ic,iv_loc)*cap_h_c(ic,iv_loc)+& 
             omega_h(ic,iv_loc)*h_x(ic,iv_loc)+&
             sum(omega_s(:,ic,iv_loc)*field(ir,it,:)) 

     enddo
  enddo

  ! TRAPPING TERM

  if (collision_model == 0) call cgyro_rhs_trap(ij)

  call timer_lib_out('rhs')

  ! Nonlinear evaluation [f,g]

  call timer_lib_in('rhs_nl')

  !call parallel_slib_f(h_x,h_nl)
  ! if (nonlinear_flag == 1) call cgyro_rhs_nl(ij)
  !call parallel_slib_r(h_nl,h_x)

  call timer_lib_out('rhs_nl')

end subroutine cgyro_rhs

!==========================================================================

subroutine cgyro_rhs_trap(ij)

  use parallel_lib

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer, intent(in) :: ij
  integer :: is,ir,it,ie,ix,jx
  complex :: val

  call parallel_lib_r(transpose(cap_h_c),cap_h_v)
  cap_h_v_prime(:,:) = (0.0,0.0)
  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1
     it = it_c(ic)
     ir = ir_c(ic)
     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        do jx=1,n_xi
           cap_h_v_prime(ic_loc,iv) = cap_h_v_prime(ic_loc,iv) &
                +xi_deriv_mat(ix,jx)*cap_h_v(ic_loc,iv_v(ie,jx,is))
        enddo
     enddo
  enddo

  ! Now have cap_h_v(ic_loc,iv)   

  call parallel_lib_f(cap_h_v_prime,cap_h_ct)
  cap_h_c = transpose(cap_h_ct)

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        val = omega_trap(it,is)*sqrt(energy(ie))*(1.0-xi(ix)**2) 

        rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
             -val*cap_h_c(ic,iv_loc)
     enddo
  enddo

end subroutine cgyro_rhs_trap

!==========================================================================

subroutine cgyro_rhs_nl(ij)

  use parallel_lib

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer, intent(in) :: ij


end subroutine cgyro_rhs_nl
