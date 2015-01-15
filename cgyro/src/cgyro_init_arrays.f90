subroutine cgyro_init_arrays

  use timer_lib
  use mpi

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  real, external :: BESJ0
  real :: arg,ang
  integer :: ir,it,is,ie,ix
  integer :: jr,jt,id
  complex :: thfac
  real, dimension(n_radial,n_theta) :: sum_loc
  real, dimension(nv_loc) :: vfac

  call timer_lib_in('init_arrays')

  !-------------------------------------------------------------------------
  ! Distributed Bessel-function Gyroaverages

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        arg = k_perp(it,ir)*rho*vth(is)*mass(is)/(z(is)*Bmag(it)) &
             *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)

        j0_c(ic,iv_loc) = BESJ0(abs(arg))

     enddo
  enddo

  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1

     it = it_c(ic)
     ir = ir_c(ic)

     do iv=1,nv

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        arg = k_perp(it,ir)*rho*vth(is)*mass(is)/(z(is)*Bmag(it)) &
             *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)

        j0_v(ic_loc,iv) = BESJ0(abs(arg))

     enddo
  enddo
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Field equation prefactors, sums.
  !
  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     vfac(iv_loc) = 0.5*w_xi(ix)*w_e(ie)*z(is)**2/temp(is)*dens(is)

  enddo

  sum_den_h = 0.0
  do is=1,n_species
     do ie=1,n_energy
        do ix=1,n_xi
           sum_den_h = sum_den_h+0.5*w_xi(ix)*w_e(ie)*z(is)**2/temp(is)*dens(is)
        enddo
     enddo
  enddo

  if (ae_flag == 1) then
     sum_den_h = sum_den_h+dens_ele/temp_ele
  endif

  allocate(sum_den_x(n_radial,n_theta))
  sum_loc(:,:) = 0.0

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        sum_loc(ir,it) = sum_loc(ir,it)+vfac(iv_loc)*(1.0-j0_c(ic,iv_loc)**2) 

     enddo
  enddo

  call MPI_ALLREDUCE(sum_loc,&
       sum_den_x,&
       size(sum_den_x),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (ae_flag == 1) then
     sum_den_x(:,:) = sum_den_x(:,:) + dens_ele / temp_ele
  endif

  ! Pre-factors for Ampere eqn

  if (n_field > 1) then

     allocate(sum_cur_x(n_radial,n_theta))
     sum_loc(:,:)  = 0.0

     iv_loc = 0
     do iv=nv1,nv2

        iv_loc = iv_loc+1

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        do ic=1,nc

           ir = ir_c(ic) 
           it = it_c(ic)

           sum_loc(ir,it) = sum_loc(ir,it)+vfac(iv_loc) &
                *xi(ix)**2*2.0*energy(ie)*vth(is)**2*j0_c(ic,iv_loc)**2 
        enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          sum_cur_x,&
          size(sum_cur_x),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

  endif
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Zonal flow with adiabatic electrons:
  !
  if (ae_flag == 1) then

     allocate(hzf(n_radial,n_theta,n_theta))
     hzf(:,:,:) = 0.0      
     do ir=1,n_radial
        do it=1,n_theta
           hzf(ir,it,it) = k_perp(it,ir)**2 * lambda_debye**2 &
                * dens_ele / temp_ele + sum_den_h
           do jt=1,n_theta
              hzf(ir,it,jt) = hzf(ir,it,jt) &
                   - dens_ele / temp_ele * w_theta(jt)
           enddo
        enddo
     enddo

     allocate(work(n_theta))
     allocate(i_piv(n_theta))
     do ir=1,n_radial
        call DGETRF(n_theta,n_theta,hzf(ir,:,:),n_theta,i_piv,info)
        call DGETRI(n_theta,hzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
     enddo
     deallocate(i_piv)
     deallocate(work)

     allocate(xzf(n_radial,n_theta,n_theta))
     xzf(:,:,:) = 0.0     
     do ir=1,n_radial
        do it=1,n_theta
           xzf(ir,it,it) = k_perp(it,ir)**2 * lambda_debye**2 &
                * dens_ele / temp_ele + sum_den_x(ir,it)
           do jt=1,n_theta
              xzf(ir,it,jt) = xzf(ir,it,jt) &
                   - dens_ele / temp_ele * w_theta(jt)
           enddo
        enddo
     enddo

     allocate(work(n_theta))
     allocate(i_piv(n_theta))
     do ir=1,n_radial
        call DGETRF(n_theta,n_theta,xzf(ir,:,:),n_theta,i_piv,info)
        call DGETRI(n_theta,xzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
     enddo
     deallocate(i_piv)
     deallocate(work)

     allocate(pvec_in(n_theta))
     allocate(pvec_outr(n_theta))
     allocate(pvec_outi(n_theta))

  endif
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Streaming arrays
  !
  ! cyclic index (for theta-periodicity)
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo
  ! coefficients for 4th order centered derivative
  cderiv(-3) =  0.0
  cderiv(-2) =  1.0 / (12.0 * d_theta)
  cderiv(-1) = -8.0 / (12.0 * d_theta)
  cderiv(0)  =  0.0 / (12.0 * d_theta)
  cderiv(1)  =  8.0 / (12.0 * d_theta)
  cderiv(2)  = -1.0 / (12.0 * d_theta)
  cderiv(3)  = 0.0
  ! coefficients for 4th order filter for 3rd order upwinded derivative
  uderiv(-3) = 0.0
  uderiv(-2) =  1.0 / (12.0 * d_theta)
  uderiv(-1) = -4.0 / (12.0 * d_theta)
  uderiv(0)  =  6.0 / (12.0 * d_theta)
  uderiv(1)  = -4.0 / (12.0 * d_theta)
  uderiv(2)  =  1.0 / (12.0 * d_theta)
  uderiv(3)  = 0.0

  ! Indices for parallel streaming with upwinding
  if (zf_test_flag == 1) then

     ! 4th order upwind for zonal flow test

     do ir=1,n_radial
        do it=1,n_theta
           do id=-3,3
              dtheta(ir,it,id)    = cderiv(id)
              dtheta_up(ir,it,id) = uderiv(id)*up_theta
              rcyc(ir,it,id)      = ir
           enddo
        enddo
     enddo

  else if (weno_flag == 0) then

     ! 4th order upwind 

     do ir=1,n_radial
        do it=1,n_theta
           do id=-3,3
              jt = thcyc(it+id)
              if (it+id < 1) then
                 thfac = exp(2*pi*i_c*k_theta*rmin)
                 jr = ir-n*box_size
                 if (jr < 1) then
                    jr = jr+n_radial
                 endif
              else if (it+id > n_theta) then
                 thfac = exp(-2*pi*i_c*k_theta*rmin)
                 jr = ir+n*box_size
                 if (jr > n_radial) then
                    jr = jr-n_radial
                 endif
              else
                 thfac = (1.0,0.0)
                 jr = ir
              endif
              dtheta(ir,it,id)    = cderiv(id)*thfac
              dtheta_up(ir,it,id) = uderiv(id)*thfac*up_theta
              rcyc(ir,it,id)      = jr
           enddo
        enddo
     enddo

  else

     ! 5th order WENO

     do ir=1,n_radial
        do it=1,n_theta
           do id=-3,3
              jt = thcyc(it+id)
              if (it+id < 1) then
                 thfac = exp(2*pi*i_c*k_theta*rmin)
                 jr = ir-n*box_size
                 if (jr < 1) then
                    jr = jr+n_radial
                 endif
              else if (it+id > n_theta) then
                 thfac = exp(-2*pi*i_c*k_theta*rmin)
                 jr = ir+n*box_size
                 if (jr > n_radial) then
                    jr = jr-n_radial
                 endif
              else
                 thfac = (1.0,0.0)
                 jr = ir
              endif
              dtheta(ir,it,id) = thfac
              rcyc(ir,it,id)   = jr
           enddo
        enddo
     enddo

  endif


  ! Streaming coefficients (for speed optimization)

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        ! omega_rdrift
        omega_cap_h(ic,iv_loc) = -omega_rdrift(it,is)*energy(ie)*&
             (1.0 + xi(ix)**2)*(2.0*pi*i_c*px(ir)/length) 

        ! omega_dalpha
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) &
             -omega_adrift(it,is)*energy(ie)*(1.0 + xi(ix)**2)*i_c*k_theta

        ! omega_dalpha - pressure component
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) &
             -omega_aprdrift(it,is)*energy(ie)*xi(ix)**2*i_c*k_theta

        ! radial upwind
        omega_h(ic,iv_loc) = &
             -abs(omega_rdrift(it,is))*energy(ie)*(1.0 + xi(ix)**2)*up_radial & 
             *(2.0*px(ir)/(1.0*n_radial))**(up_radial_n-1.0) &
             *(2.0*pi*px(ir)/length)

        ! omega_star
        omega_s(1,ic,iv_loc) = &
             -i_c*k_theta*rho*(dlnndr(is)+dlntdr(is)*(energy(ie)-1.5)) &
             *j0_c(ic,iv_loc)

        if (n_field > 1) then
           omega_s(2,ic,iv_loc) = -omega_s(1,ic,iv_loc)* &
                xi(ix)*sqrt(2.0*energy(ie))*vth(is)
        endif

     enddo
  enddo

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Initial conditions
  !
  h_x(:,:) = (0.0,0.0)
  !
  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        if (n == 0) then

           ! Zonal-flow initial condition

           if (zf_test_flag == 1) then
              if (is == 1 .and. abs(px(ir)) == 1) then
                 h_x(ic,iv_loc) = 1e-6
              endif
           else
              ! CAUTION: Need f(p) = conjg[ f(-p) ] for n=0
              arg = abs(px(ir))/real(n_radial)
              h_x(ic,iv_loc) = arg*rho*exp(-4.0*arg)
              if (ir == 1) h_x(ic,iv_loc) = (0.0,0.0)
           endif

        else 

           ! Exponential in ballooning angle.

           if (n_toroidal == 1) then
              if (is == 1) then
                 ang = theta(it)+2*pi*px(ir)
                 h_x(ic,iv_loc) = rho*exp(-(ang/2)**2) 
              endif
           else
              h_x(ic,iv_loc) = amp*rho*exp(-px(ir)*4.0/n_radial) 
           endif

        endif

     enddo
  enddo

  call cgyro_field_c

  field_old = field
  !-------------------------------------------------------------------------

  call timer_lib_out('init_arrays')

end subroutine cgyro_init_arrays
