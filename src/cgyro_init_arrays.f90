subroutine cgyro_init_arrays

  use mpi

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  real, external :: BESJ0
  real :: arg
  integer :: is,ir,it,ie,ix
  integer :: jt
  real, dimension(n_radial,n_theta) :: sum_loc


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

  sum_den_h = 0.0
  do is=1,n_species
     do ie=1,n_energy
        do ix=1,n_xi
           sum_den_h = sum_den_h &
                + 0.5 * w_xi(ix) &
                * z(is)**2/temp(is) *dens(is) * w_e(ie)
        enddo
     enddo
  enddo

  if (ae_flag == 1) then
     sum_den_h = sum_den_h + dens_ele / temp_ele
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

        sum_loc(ir,it) = sum_loc(ir,it) &
             +0.5*w_xi(ix)*w_e(ie)*z(is)**2/temp(is)*dens(is) &
             *(1.0-j0_c(ic,iv_loc)**2) 

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

  if (zf_test_flag == 1 .and. ae_flag == 1) then

     ! Zonal flow with adiabatic electrons:

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

           sum_loc(ir,it) = sum_loc(ir,it) &
                +0.5*w_xi(ix)*w_e(ie)*xi(ix)**2*2.0*energy(ie) &
                *vth(is)**2*z(is)**2/temp(is)*dens(is) & 
                *j0_c(ic,iv_loc)**2 
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

end subroutine cgyro_init_arrays
