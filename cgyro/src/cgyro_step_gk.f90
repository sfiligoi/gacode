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
  integer :: is, ir, it, ie, ix
  integer :: id, jt, jr, jc
  real :: rval
  complex :: rhs_stream,fhp,fhm
  real, dimension(3) :: gamma=(/0.1,0.6,0.3/),beta,w,wt
  complex, dimension(3) :: fp,fm
  complex, dimension(-3:3) :: f

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
        rhs_stream = 0.0

        if (weno_flag == 0) then

           ! Upwind3

           do id=-2,2
              jt = thcyc(it+id)
              jr = rcyc(ir,it,id)
              jc = ic_c(jr,jt)
              rhs_stream = rhs_stream &
                   -rval*dtheta(ir,it,id)*cap_h_c(jc,iv_loc)  &
                   -abs(rval)*dtheta_up(ir,it,id)*h_x(jc,iv_loc)
           enddo

        else

           ! WENO5

           do id=-3,3
              jt = thcyc(it+id)
              jr = rcyc(ir,it,id)
              jc = ic_c(jr,jt)
              ! Multiply by appropriate phase factor
              f(id) = h_x(jc,iv_loc)*dtheta(ir,it,id)
              ! 5-point centered theta-derivative of fields
              rhs_stream = rhs_stream-rval*cderiv(id)*dtheta(ir,it,id)* &
                   (cap_h_c(jc,iv_loc)-h_x(jc,iv_loc))
           enddo

          if (rval > 0.0) then

              fp(1) = ( 2*f(-2)-7*f(-1)+11*f(0))/6.0           
              fp(2) = (-1*f(-1)+5*f( 0)+ 2*f(1))/6.0           
              fp(3) = ( 2*f( 0)+5*f( 1)- 1*f(2))/6.0           

              beta(1) = 13.0/12.0*abs(f(-2)-2*f(-1)+f(0))**2 &
                   +0.25*abs(f(-2)-4*f(-1)+3*f(0))**2
              beta(2) = 13.0/12.0*abs(f(-1)-2*f(0)+f(1))**2 &
                   +0.25*abs(f(-1)-f(1))**2
              beta(3) = 13.0/12.0*abs(f(0)-2*f(1)+f(2))**2 &
                   +0.25*abs(3*f(0)-4*f(1)+f(2))**2

              wt(:) = gamma(:)/(1e-6+beta(:))**2
              w(:)  = wt(:)/sum(wt)         

              fhp = sum(w(:)*fp(:))

              fm(1) = ( 2*f(-3)-7*f(-2)+11*f(-1))/6.0           
              fm(2) = (-1*f(-2)+5*f(-1)+ 2*f( 0))/6.0           
              fm(3) = ( 2*f(-1)+5*f( 0)- 1*f( 1))/6.0           

              beta(1) = 13.0/12.0*abs(f(-3)-2*f(-2)+f(-1))**2 &
                   +0.25*abs(f(-3)-4*f(-2)+3*f(-1))**2
              beta(2) = 13.0/12.0*abs(f(-2)-2*f(-1)+f(0))**2 &
                   +0.25*abs(f(-2)-f(0))**2
              beta(3) = 13.0/12.0*abs(f(-1)-2*f(0)+f(1))**2 &
                   +0.25*abs(3*f(-1)-4*f(0)+f(1))**2

              wt(:) = gamma(:)/(1e-6+beta(:))**2
              w(:)  = wt(:)/sum(wt)         

              fhm = sum(w(:)*fm(:))

           else

              ! fm(i) -> fp(-i) , fp(i) -> fm(-i)

              fm(1) = ( 2*f(+2)-7*f(+1)+11*f( 0))/6.0           
              fm(2) = (-1*f(+1)+5*f( 0)+ 2*f(-1))/6.0           
              fm(3) = ( 2*f( 0)+5*f(-1)- 1*f(-2))/6.0           

              beta(1) = 13.0/12.0*abs(f( 2)-2*f( 1)+f(0))**2 &
                   +0.25*abs(f( 2)-4*f( 1)+3*f(0))**2
              beta(2) = 13.0/12.0*abs(f( 1)-2*f(0)+f(-1))**2 &
                   +0.25*abs(f(1)-f(-1))**2
              beta(3) = 13.0/12.0*abs(f(0)-2*f(-1)+f(-2))**2 &
                   +0.25*abs(3*f(0)-4*f(-1)+f(-2))**2

              wt(:) = gamma(:)/(1e-6+beta(:))**2
              w(:)  = wt(:)/sum(wt)         

              fhm = sum(w(:)*fm(:))

              fp(1) = ( 2*f( 3)-7*f( 2)+11*f( 1))/6.0           
              fp(2) = (-1*f( 2)+5*f( 1)+ 2*f( 0))/6.0           
              fp(3) = ( 2*f( 1)+5*f( 0)- 1*f(-1))/6.0           

              beta(1) = 13.0/12.0*abs(f( 3)-2*f( 2)+f( 1))**2 &
                   +0.25*abs(f( 3)-4*f( 2)+3*f( 1))**2
              beta(2) = 13.0/12.0*abs(f( 2)-2*f( 1)+f(0))**2 &
                   +0.25*abs(f(2)-f(0))**2
              beta(3) = 13.0/12.0*abs(f(1)-2*f(0)+f(-1))**2 &
                   +0.25*abs(3*f(1)-4*f(0)+f(-1))**2

              wt(:) = gamma(:)/(1e-6+beta(:))**2
              w(:)  = wt(:)/sum(wt)         

              fhp = sum(w(:)*fp(:))

           endif

           rhs_stream = rhs_stream-rval*(fhp-fhm)/d_theta

        endif

        ! Diagonal terms
        rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc)+&
             rhs_stream+&
             omega_cap_h(ic,iv_loc)*cap_h_c(ic,iv_loc)+& 
             omega_h(ic,iv_loc)*h_x(ic,iv_loc)+&
             sum(omega_s(:,ic,iv_loc)*field(ir,it,:)) 

     enddo
  enddo

  ! TRAPPING TERM
  if (collision_model == 0) call cgyro_rhs_trap(ij)

  call timer_lib_out('rhs')

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

subroutine filter(f)

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: ie,ir,is,it,ix

  complex, dimension(nc,nv_loc), intent(in) :: f
  complex, dimension(n_radial):: fs

  fs = (0.0,0.0)

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        fs(ir) = fs(ir)+f(ic,iv_loc)

     enddo
  enddo

  if (n==0) print *
  do ir=1,n_radial
     if (n == 0) print *,px(ir),fs(ir)
  enddo

end subroutine filter
