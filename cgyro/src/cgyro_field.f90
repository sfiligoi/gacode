!============================================================================================
! Velocity (configuration-distributed) field solve

subroutine cgyro_field_v

  use mpi
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is, ie, ix, ir, it
  complex :: fac

  call timer_lib_in('field_H')

  field_loc(:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of H

  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1

     do iv=1,nv

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        fac = w_e(ie)*w_xi(ix)*z(is)*dens(is)*cap_h_v(ic_loc,iv)
        field_loc(:,ic) = field_loc(:,ic)+fac*jvec_v(:,ic_loc,iv) 

     enddo
  enddo

  call MPI_ALLREDUCE(field_loc(:,:),&
       field(:,:),&
       size(field(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! Poisson LHS factors

  if (n == 0 .and. ae_flag == 1) then

     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
           pvec_in(:) = real(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
           pvec_in(:) = aimag(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
           field(1,ic_c(ir,:)) = pvec_outr(:) + i_c * pvec_outi(:)
        endif
     enddo

  else

     do ir=1,n_radial
        if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
           do it=1,n_theta
              field(1,ic_c(ir,it)) = field(1,ic_c(ir,it)) &
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
           field(2,ic_c(ir,:)) = 0.0
        else
           do it=1,n_theta
              field(2,ic_c(ir,it)) = field(2,ic_c(ir,it)) &
                   /(-2.0*k_perp(it,ir)**2*rho**2 &
                   /betae_unit*dens_ele*temp_ele)
           enddo
        endif
     enddo

     ! Ampere Bpar LHS factors
     
     if (n_field > 2) then
        do ir=1,n_radial
           if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
              field(3,ic_c(ir,:)) = 0.0
           else
              field(3,ic_c(ir,:)) = field(3,ic_c(ir,:))*(-0.5*betae_unit)/(dens_ele*temp_ele)
           endif
        enddo
     endif

  endif

  call timer_lib_out('field_H')

end subroutine cgyro_field_v

!============================================================================================
! Configuration (velocity-distributed) field solve

subroutine cgyro_field_c

  use mpi
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is, ie, ix, ir, it
  complex :: fac

  call timer_lib_in('field_h')

  field_loc(:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        fac = w_e(ie)*w_xi(ix)*z(is)*dens(is)*h_x(ic,iv_loc)
        field_loc(:,ic) = field_loc(:,ic)+jvec_c(:,ic,iv_loc)*fac
     enddo
  enddo

  call MPI_ALLREDUCE(field_loc(:,:),&
       field(:,:),&
       size(field(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (n_field > 1) then
     field(2,:) = -field(2,:)
  endif
  if (n_field > 2) then
     field(3,:) = field(3,:)*(-0.5*betae_unit)/(dens_ele*temp_ele)
  endif

  ! Poisson LHS factors

  if (n == 0 .and. ae_flag == 1) then

     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
           pvec_in(:) = real(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
           pvec_in(:) = aimag(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
           field(1,ic_c(ir,:)) = pvec_outr(:) + i_c * pvec_outi(:)
        endif
     enddo

  else

     do ir=1,n_radial

        if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
           if (n_field > 2) then
              field(3,ic_c(ir,:)) = 0.0
           endif
        else

           if (n_field == 3) then
              do it=1,n_theta
                 fac = field(1,ic_c(ir,it))

                 field(1,ic_c(ir,it)) =  poisson_pb22(ir,it)*field(1,ic_c(ir,it)) &
                      - poisson_pb12(ir,it)*field(3,ic_c(ir,it))

                 field(3,ic_c(ir,it)) =  -poisson_pb21(ir,it)*fac &
                      + poisson_pb11(ir,it)*field(3,ic_c(ir,it))
              enddo

           else
              do it=1,n_theta
                 field(1,ic_c(ir,it)) = field(1,ic_c(ir,it)) &
                      /(k_perp(it,ir)**2*lambda_debye**2* &
                      dens_ele/temp_ele+sum_den_x(ir,it))
              enddo
           endif

        endif
     enddo
  endif

  ! Ampere LHS factors

  if (n_field > 1) then
     do ir=1,n_radial
        if (n == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(2,ic_c(ir,:)) = 0.0
        else
           do it=1,n_theta
              field(2,ic_c(ir,it)) = field(2,ic_c(ir,it)) &
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

     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        cap_h_c(ic,iv_loc) = h_x(ic,iv_loc)+psi(ic,iv_loc)*z(is)/temp(is)
     enddo
  enddo

  call timer_lib_out('field_h')

end subroutine cgyro_field_c
