!============================================================================================
! Velocity-space field solve

subroutine cgyro_field_v

  use mpi
  use timer_lib

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: is, ie, ix, ir, it
  complex :: fac

  call timer_lib_in('field_v')

  field_loc(:,:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of H

  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1

     it = it_c(ic)
     ir = ir_c(ic)

     do iv=1,nv

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        fac = w_e(ie)*0.5*w_xi(ix)*z(is)*dens(is)*&
             j0_v(ic_loc,iv)*cap_h_v(ic_loc,iv)

        field_loc(ir,it,1) = field_loc(ir,it,1)+fac 

        if (n_field > 1) then
           field_loc(ir,it,2) = field_loc(ir,it,2) + fac &
                *xi(ix)*sqrt(2.0*energy(ie))*vth(is)
        endif

     enddo
  enddo

  call MPI_ALLREDUCE(field_loc(:,:,:),&
       field(:,:,:),&
       size(field(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! Poisson LHS factors

  if (n == 0 .and. ae_flag == 1) then

     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(ir,:,1) = 0.0
        else
           pvec_in(:) = real(field(ir,:,1))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
           pvec_in(:) = imag(field(ir,:,1))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
           field(ir,:,1) = pvec_outr(:) + i_c * pvec_outi(:)
        endif
     enddo

  else

     do ir=1,n_radial
        if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(ir,:,1) = 0.0
        else
           do it=1,n_theta
              field(ir,it,1) = field(ir,it,1) &
                   /(k_perp(it,ir)**2*lambda_debye**2* &
                   dens_ele/temp_ele+sum_den_h)
           enddo
        endif
     enddo

  endif

  ! Ampere LHS factors

  if (n_field > 1) then
     do ir=1,n_radial
        if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(ir,:,2) = 0.0
        else
           do it=1,n_theta
              field(ir,it,2) = field(ir,it,2) &
                   /(2.0*k_perp(it,ir)**2*rho**2 &
                   /betae_unit*dens_ele*temp_ele)
           enddo
        endif
     enddo
  endif

  call timer_lib_out('field_v')

end subroutine cgyro_field_v

!============================================================================================
! Configuration-space field solve

subroutine cgyro_field_c

  use mpi
  use timer_lib

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: is, ie, ix, ir, it
  complex :: fac
  complex :: efac(n_field)

  call timer_lib_in('field_c')

  field_loc(:,:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic)
        it = it_c(ic)

        fac = w_e(ie)*0.5*w_xi(ix)*z(is)*dens(is)* &
             j0_c(ic,iv_loc)*h_x(ic,iv_loc)

        field_loc(ir,it,1) = field_loc(ir,it,1)+fac 

        if (n_field > 1) then
           field_loc(ir,it,2) = field_loc(ir,it,2) + &
                fac*xi(ix)*sqrt(2.0*energy(ie))*vth(is)
        endif

     enddo
  enddo

  call MPI_ALLREDUCE(field_loc(:,:,:),&
       field(:,:,:),&
       size(field(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! Poisson LHS factors

  if (n == 0 .and. ae_flag == 1) then

     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(ir,:,1) = 0.0
        else
           pvec_in(:) = real(field(ir,:,1))
           call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
           pvec_in(:) = imag(field(ir,:,1))
           call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
           field(ir,:,1) = pvec_outr(:) + i_c * pvec_outi(:)
        endif
     enddo

  else

     do ir=1,n_radial
        if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(ir,:,1) = 0.0
        else
           do it=1,n_theta
              field(ir,it,1) = field(ir,it,1) &
                   /(k_perp(it,ir)**2*lambda_debye**2* &
                   dens_ele/temp_ele+sum_den_x(ir,it))
           enddo
        endif
     enddo

  endif

  ! Ampere LHS factors

  if (n_field > 1) then
     do ir=1,n_radial
        if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(ir,:,2) = 0.0
        else
           do it=1,n_theta
              field(ir,it,2) = field(ir,it,2) &
                   /(2.0*k_perp(it,ir)**2*rho**2/betae_unit & 
                   *dens_ele*temp_ele+sum_cur_x(ir,it))
           enddo
        endif
     enddo
  endif

  ! Compute H given h and [phi(h), apar(h)]

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     efac(1) = 1.0
     if (n_field > 1) then
        efac(2) = -xi(ix)*sqrt(2.0*energy(ie))*vth(is)
     endif

     do ic=1,nc

        ir = ir_c(ic)
        it = it_c(ic)

        psi(ic,iv_loc) = j0_c(ic,iv_loc)*sum(efac(:)*field(ir,it,:))

        cap_h_c(ic,iv_loc) = h_x(ic,iv_loc)+z(is)*psi(ic,iv_loc)/temp(is)              

     enddo
  enddo

  call timer_lib_out('field_c')

end subroutine cgyro_field_c
