subroutine cgyro_init_arrays

  use timer_lib
  use mpi
  use cgyro_globals
  use parallel_lib

  implicit none

  real :: arg
  real :: efac
  integer :: ir,it,is,ie,ix
  integer :: jr,jt,id, ccw_fac
  complex :: thfac, carg
  real, dimension(nc) :: sum_loc
  real, dimension(n_species,nc) :: res_loc
  real, dimension(nv_loc) :: vfac
  real, dimension(n_radial) :: u
  real, dimension(:,:), allocatable :: jloc_c
  real, dimension(:), allocatable :: pb11,pb12,pb21,pb22

  !-------------------------------------------------------------------------
  ! Distributed Bessel-function Gyroaverages

  allocate(jloc_c(2,nc))

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        it = it_c(ic)

        arg = k_perp(ic)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
             *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)

        ! Need this for (Phi, A_parallel) terms in GK and field equations

        jloc_c(1,ic) = bessel_j0(abs(arg))

        ! Needed for B_parallel in GK and field equations

        jloc_c(2,ic) = 0.5*(jloc_c(1,ic) + bessel_jn(2,abs(arg)))/bmag(it)

     enddo

     efac = 1.0
     jvec_c(1,:,iv_loc) = efac*jloc_c(1,:)
     if (n_field > 1) then
        efac = -xi(ix)*sqrt(2.0*energy(ie))*vth(is)
        jvec_c(2,:,iv_loc) = efac*jloc_c(1,:)
        if (n_field > 2) then
           efac = 2.0*energy(ie)*(1-xi(ix)**2)*temp(is)/z(is)
           jvec_c(3,:,iv_loc) = efac*jloc_c(2,:)
        endif
     endif

  enddo
  deallocate(jloc_c)
  call parallel_lib_rtrans_real(jvec_c(1,:,:),jvec_v(1,:,:))
  if (n_field > 1) call parallel_lib_rtrans_real(jvec_c(2,:,:),jvec_v(2,:,:))
  if (n_field > 2) call parallel_lib_rtrans_real(jvec_c(3,:,:),jvec_v(3,:,:))
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Field equation prefactors, sums.
  !
  do iv=nv1,nv2

     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     vfac(iv_loc) = w_xi(ix)*w_e(ie)*z(is)**2/temp(is)*dens(is)

  enddo

  sum_den_h = 0.0
  do is=1,n_species
     do ie=1,n_energy
        do ix=1,n_xi
           sum_den_h = sum_den_h+w_xi(ix)*w_e(ie)*z(is)**2/temp(is)*dens(is)
        enddo
     enddo
  enddo

  if (ae_flag == 1) then
     sum_den_h = sum_den_h+dens_ele/temp_ele
  endif

  allocate(sum_den_x(nc))
  sum_loc(:) = 0.0

  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        sum_loc(ic) = sum_loc(ic)+vfac(iv_loc)*(1.0-jvec_c(1,ic,iv_loc)**2) 
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
     sum_den_x(:) = sum_den_x(:)+dens_ele/temp_ele
  endif
  !------------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Conservative upwind factor
  !
  allocate(res_norm(n_species,nc))

  res_loc(:,:) = 0.0

  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        res_loc(is,ic) = res_loc(is,ic)+w_xi(ix)*w_e(ie)*jvec_c(1,ic,iv_loc)**2 
     enddo
  enddo

  call MPI_ALLREDUCE(res_loc,&
       res_norm,&
       size(res_norm),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
  !------------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Field-solve coefficients (i.e., final numerical factors).
  !
  do ic=1,nc
     ir = ir_c(ic) 
     it = it_c(ic)
     if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
        fcoef(:,ic) = 0.0
     else
        fcoef(1,ic) = 1.0/(k_perp(ic)**2*lambda_debye**2*dens_ele/temp_ele+sum_den_h)
        if (n_field > 1) fcoef(2,ic) = 1.0/(-2.0*k_perp(ic)**2* &
             rho**2/betae_unit*dens_ele*temp_ele)
        if (n_field > 2) fcoef(3,ic) = -betae_unit/(2.0*dens_ele*temp_ele)
     endif
  enddo

  if (n_field > 1) then

     allocate(sum_cur_x(nc))
     sum_loc(:) = 0.0

     iv_loc = 0
     do iv=nv1,nv2
        iv_loc = iv_loc+1
        do ic=1,nc
           sum_loc(ic) = sum_loc(ic)+vfac(iv_loc)*jvec_c(2,ic,iv_loc)**2
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

  if (n_field == 1 .or. n_field == 2) then
     do ic=1,nc
        if (k_perp(ic) > 0.0) then
           gcoef(1,ic) = 1.0/(k_perp(ic)**2*lambda_debye**2*&
                dens_ele/temp_ele+sum_den_x(ic))
        endif
     enddo
  endif

  if (n_field > 1) then
     do ic=1,nc
        if (k_perp(ic) > 0.0) then
           gcoef(2,ic) = 1.0/(-2.0*k_perp(ic)**2*&
                rho**2/betae_unit*dens_ele*temp_ele-sum_cur_x(ic))
        endif
     enddo
  endif

  if (n_field > 2) then
     allocate(pb11(nc))
     allocate(pb12(nc))
     allocate(pb21(nc))
     allocate(pb22(nc))

     do ic=1,nc
        pb11(ic) = k_perp(ic)**2*lambda_debye**2* &
             dens_ele/temp_ele+sum_den_x(ic)
     enddo

     sum_loc(:)  = 0.0
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        do ic=1,nc
           ir = ir_c(ic) 
           it = it_c(ic)
           sum_loc(ic) = sum_loc(ic)-w_xi(ix)*w_e(ie)*dens(is) &
                *z(is)*jvec_c(1,ic,iv_loc)*jvec_c(3,ic,iv_loc) &
                *z(is)/temp(is)
        enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          pb12,&
          size(pb12),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     pb21(:) = pb12(:)*betae_unit/(-2*dens_ele*temp_ele)

     sum_loc(:)  = 0.0
     iv_loc = 0
     do iv=nv1,nv2
        iv_loc = iv_loc+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        do ic=1,nc
           ir = ir_c(ic) 
           it = it_c(ic)
           sum_loc(ic) = sum_loc(ic)+w_xi(ix)*w_e(ie)*dens(is) &
                *temp(is)*jvec_c(3,ic,iv_loc)**2 &
                *(z(is)/temp(is))**2
        enddo
     enddo
     call MPI_ALLREDUCE(sum_loc,&
          pb22,&
          size(pb22),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     pb22(:) = 1.0-pb22(:)*betae_unit/(-2*dens_ele*temp_ele) 

     ! Determinant
     do ic=1,nc
        if (k_perp(ic) > 0.0) then
           sum_loc(ic) = pb11(ic)*pb22(ic)-pb12(ic)*pb21(ic)
        else
           sum_loc(ic) = 1.0
        endif
     enddo

     pb11 = pb11/sum_loc
     pb12 = pb12/sum_loc
     pb21 = pb21/sum_loc
     pb22 = pb22/sum_loc

     gcoef(3,:) = pb11
     gcoef(1,:) = pb22
     gcoef(4,:) = -pb12
     gcoef(5,:) = -pb21

     deallocate(pb11)
     deallocate(pb12)
     deallocate(pb21)
     deallocate(pb22)
  endif

  ! Set selected zeros
  do ic=1,nc
     ir = ir_c(ic) 
     if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
        gcoef(:,ic) = 0.0
     endif
  enddo
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Zonal flow with adiabatic electrons:
  !
  if (n == 0 .and. ae_flag == 1) then

     allocate(hzf(n_radial,n_theta,n_theta))
     hzf(:,:,:) = 0.0      
     do ir=1,n_radial
        do it=1,n_theta
           hzf(ir,it,it) = k_perp(ic_c(ir,it))**2 * lambda_debye**2 &
                * dens_ele / temp_ele + sum_den_h
           do jt=1,n_theta
              hzf(ir,it,jt) = hzf(ir,it,jt) &
                   - dens_ele/temp_ele*w_theta(jt)
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
           xzf(ir,it,it) = k_perp(ic_c(ir,it))**2 * lambda_debye**2 &
                * dens_ele / temp_ele + sum_den_x(ic_c(ir,it))
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
  ! cyclic functions (for radial and theta periodicity)
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo
  do ir=1,n_radial
     rcyc(ir-n_radial) = ir
     rcyc(ir) = ir
     rcyc(ir+n_radial) = ir
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

  if (ipccw*btccw < 0) then
     ccw_fac = -1
  else
     ccw_fac = 1
  endif


  ! Indices for parallel streaming with upwinding
  do ir=1,n_radial
     do it=1,n_theta
        do id=-nup_theta,nup_theta
           jt = thcyc(it+id)
           if (it+id < 1) then
              thfac = exp(2*pi*i_c*k_theta*rmin)
              jr = rcyc(ir-n*box_size*ccw_fac)
           else if (it+id > n_theta) then
              thfac = exp(-2*pi*i_c*k_theta*rmin)
              jr = rcyc(ir+n*box_size*ccw_fac)
           else
              thfac = (1.0,0.0)
              jr = ir
           endif
           dtheta(ic_c(ir,it),id)    = cderiv(id)*thfac
           dtheta_up(ic_c(ir,it),id) = uderiv(id)*thfac*up_theta
           icd_c(ic_c(ir,it),id)     = ic_c(jr,thcyc(it+id))
        enddo
     enddo
  enddo
!$acc enter data copyin(dtheta,dtheta_up,icd_c)

  ! Streaming coefficients (for speed optimization)

  iv_loc = 0
!$omp parallel do collapse(2) &
!$omp& private(iv,ic,iv_loc,is,ix,ie,ir,it,carg)
  do iv=nv1,nv2
     do ic=1,nc

        ! iv_loc = iv_loc+1
        iv_loc = iv-nv1+1

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)


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
             -omega_cdrift(it,is)*vel(ie)*xi(ix)*i_c*k_theta

        ! omega_rdrift
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) & 
             -omega_rdrift(it,is)*energy(ie)*&
             (1.0 + xi(ix)**2)*(2.0*pi*i_c*px(ir)/length) 

        ! radial upwind
        omega_h(ic,iv_loc) = &
             -abs(omega_rdrift(it,is))*energy(ie)*(1.0 + xi(ix)**2) &
             *up_radial *(n_radial/length) * spec_uderiv(ir)

        ! omega_star and rotation shearing
        carg = -i_c*k_theta*rho*(dlnndr(is)+dlntdr(is)*(energy(ie)-1.5)) &
             -i_c*k_theta*rho*(sqrt(2.0*energy(ie))*xi(ix)/vth(is) &
             *omega_gammap(it))

        omega_s(:,ic,iv_loc) = carg*jvec_c(:,ic,iv_loc)

     enddo

  enddo
  !-------------------------------------------------------------------------
!$acc enter data copyin(omega_cap_h,omega_h,omega_s)
end subroutine cgyro_init_arrays
