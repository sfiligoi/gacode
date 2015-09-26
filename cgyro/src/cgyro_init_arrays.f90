
subroutine cgyro_init_arrays

  use timer_lib
  use mpi
  use cgyro_globals

  implicit none

  real, external :: BESJ0
  real :: arg
  integer :: ir,it,is,ie,ix
  integer :: jr,jt,id
  complex :: thfac
  real, dimension(n_radial,n_theta) :: sum_loc
  real, dimension(nv_loc) :: vfac
  real, dimension(n_radial) :: u

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
  if (n == 0 .and. ae_flag == 1) then

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

  allocate(cderiv(-nup_theta:nup_theta))
  allocate(uderiv(-nup_theta:nup_theta))

  select case (nup_theta)

  case (1)

     ! 1st-order UPWIND

     ! 2nd-order centered derivative
     cderiv(-1) = -1.0 / (2.0 * d_theta)
     cderiv(0)  =  0.0 / (2.0 * d_theta)
     cderiv(1)  =  1.0 / (2.0 * d_theta)

     ! 2nd-derivative filter
     uderiv(-1) = -1.0 / (2.0 * d_theta)
     uderiv(0)  =  2.0 / (2.0 * d_theta)
     uderiv(1)  = -1.0 / (2.0 * d_theta)

  case (2)

     ! 3rd-order UPWIND

     ! 4th-order centered derivative
     cderiv(-2) =  1.0 / (12.0 * d_theta)
     cderiv(-1) = -8.0 / (12.0 * d_theta)
     cderiv(0)  =  0.0 / (12.0 * d_theta)
     cderiv(1)  =  8.0 / (12.0 * d_theta)
     cderiv(2)  = -1.0 / (12.0 * d_theta)

     ! 4th-derivative filter 
     uderiv(-2) =  1.0 / (12.0 * d_theta)
     uderiv(-1) = -4.0 / (12.0 * d_theta)
     uderiv(0)  =  6.0 / (12.0 * d_theta)
     uderiv(1)  = -4.0 / (12.0 * d_theta)
     uderiv(2)  =  1.0 / (12.0 * d_theta)

  case (3)

     ! 5th-order UPWIND

     ! 6th-order centered derivative
     cderiv(-3) =  -1.0 / (60.0 * d_theta)
     cderiv(-2) =   9.0 / (60.0 * d_theta)
     cderiv(-1) = -45.0 / (60.0 * d_theta)
     cderiv(0)  =   0.0 / (60.0 * d_theta)
     cderiv(1)  =  45.0 / (60.0 * d_theta)
     cderiv(2)  =  -9.0 / (60.0 * d_theta)
     cderiv(3)  =   1.0 / (60.0 * d_theta)

     ! 6th-derivative filter 
     uderiv(-3) =  -1.0 / (60.0 * d_theta)
     uderiv(-2) =   6.0 / (60.0 * d_theta)
     uderiv(-1) = -15.0 / (60.0 * d_theta)
     uderiv(0)  =  20.0 / (60.0 * d_theta)
     uderiv(1)  = -15.0 / (60.0 * d_theta)
     uderiv(2)  =   6.0 / (60.0 * d_theta)
     uderiv(3)  =  -1.0 / (60.0 * d_theta)

  end select

  allocate(spec_uderiv(n_radial))
  u(:) = (2.0*pi/n_radial)*px(:)

  select case(nup_radial)

  case (1)

     ! 2nd order spectral dissipation
     spec_uderiv(:) = 1-cos(u)

  case (2)

     ! 4th order spectral dissipation
     spec_uderiv(:) = (3-4*cos(u)+cos(2*u))/6

  case (3)

     ! 6th order spectral dissipation
     spec_uderiv(:) = (20-30*cos(u)+12*cos(2*u)-2*cos(3*u))/60

  case (4)

     ! 8th order spectral dissipation
     spec_uderiv(:) = (70-112*cos(u)+56*cos(2*u)-16*cos(3*u)+2*cos(4*u))/280

  end select

  ! Indices for parallel streaming with upwinding
  do ir=1,n_radial
     do it=1,n_theta
        do id=-nup_theta,nup_theta
           jt = thcyc(it+id)
           if (it+id < 1) then
              thfac = exp(2*pi*i_c*k_theta*rmin)
              jr = ir-n*box_size*ipccw*btccw
              if (jr < 1) then
                 jr = jr+n_radial
              endif
              if (jr > n_radial) then
                 jr = jr-n_radial
              endif
           else if (it+id > n_theta) then
              thfac = exp(-2*pi*i_c*k_theta*rmin)
              jr = ir+n*box_size*ipccw*btccw
              if (jr > n_radial) then
                 jr = jr-n_radial
              endif
              if(jr < 1) then
                 jr = jr+n_radial
              end if
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

        ! omega_dalpha
        omega_cap_h(ic,iv_loc) = &
             -omega_adrift(it,is)*energy(ie)*(1.0 + xi(ix)**2)*i_c*k_theta

        ! omega_dalpha - pressure component
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) &
             -omega_aprdrift(it,is)*energy(ie)*xi(ix)**2*i_c*k_theta

        ! omega_cdrift - mach component
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) &
             -omega_cdrift(it,is)*sqrt(energy(ie))*xi(ix)*i_c*k_theta

        ! omega_rdrift
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) & 
             -omega_rdrift(it,is)*energy(ie)*&
             (1.0 + xi(ix)**2)*(2.0*pi*i_c*px(ir)/length) 

        ! radial upwind
        !up_radial_n=n_radial
        !omega_h(ic,iv_loc) = &
        !     -abs(omega_rdrift(it,is))*energy(ie)*(1.0 + xi(ix)**2)*up_radial &
        !     *(2.0*px(ir)/(1.0*n_radial))**(up_radial_n-1.0) &
        !     *(2.0*pi*px(ir)/length)

        omega_h(ic,iv_loc) = &
             -abs(omega_rdrift(it,is))*energy(ie)*(1.0 + xi(ix)**2) &
             *up_radial *(n_radial/length) * spec_uderiv(ir)

        ! omega_star and rotation shearing
        omega_s(1,ic,iv_loc) = &
             -i_c*k_theta*rho*(dlnndr(is)+dlntdr(is)*(energy(ie)-1.5)) &
             *j0_c(ic,iv_loc)

        omega_s(1,ic,iv_loc) = omega_s(1,ic,iv_loc) &
             -i_c*k_theta*rho*(sqrt(2.0*energy(ie))*xi(ix)/vth(is) &
             *omega_gammap(it)) * j0_c(ic,iv_loc)

        if (n_field > 1) then
           omega_s(2,ic,iv_loc) = -omega_s(1,ic,iv_loc)* &
                xi(ix)*sqrt(2.0*energy(ie))*vth(is)
        endif

     enddo
  enddo

  !-------------------------------------------------------------------------

end subroutine cgyro_init_arrays
