!-----------------------------------------------------------------
! cgyro_field.f90
!
! PURPOSE:
!  Perform field solves for data distributed in both the velocity 
!  and configuration indices:
!  
!  - cgyro_field_v() (configuration distributed)
!  - cgyro_field_c() (velocity distributed)
!-----------------------------------------------------------------

! Velocity (configuration-distributed) field solve

subroutine cgyro_field_v

  use mpi
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is, ie, ix, ir
  complex :: fac

  logical, parameter :: use_dgemv = .false.
  integer :: i,j
  real, dimension(n_theta) :: pvec_inr, pvec_ini

  call timer_lib_in('field_H')

!$omp workshare
  field_loc(:,:) = (0.0,0.0)
!$omp end workshare

  ! Poisson and Ampere RHS integrals of H

!$omp parallel do    &
!$omp& private(ic,ic_loc,iv,is,ix,ie,fac)
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        fac = (w_e(ie)*w_xi(ix)*z(is)*dens(is))*cap_h_v(ic_loc,iv)
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
!$omp parallel do &
!$omp& private(ir,i,j) &
!$omp& private(pvec_in,pvec_outr,pvec_outi) &
!$omp& private(pvec_ini,pvec_inr)
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
          if (use_dgemv) then
           pvec_in(:) = real(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
           pvec_in(:) = aimag(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
           field(1,ic_c(ir,:)) = pvec_outr(:) + i_c * pvec_outi(:)
          else
           do i=1,n_theta
            pvec_inr(i) = real(field(1,ic_c(ir,i)),kind=kind(hzf))
            pvec_ini(i) = aimag(field(1,ic_c(ir,i)))
           enddo

           do i=1,n_theta
            pvec_outr(i) = 0
            pvec_outi(i) = 0
           enddo

            do j=1,n_theta
            do i=1,n_theta
              pvec_outr(i) = pvec_outr(i) + hzf(ir,i,j)*pvec_inr(j)
              pvec_outi(i) = pvec_outi(i) + hzf(ir,i,j)*pvec_ini(j)
            enddo
            enddo

            do i=1,n_theta
              field(1,ic_c(ir,i)) = cmplx(pvec_outr(i),pvec_outi(i))
            enddo

          endif
        endif
     enddo

  else

!$omp workshare
     field(:,:) = fcoef(:,:)*field(:,:)
!$omp end workshare

  endif

  call timer_lib_out('field_H')

end subroutine cgyro_field_v


! Configuration (velocity-distributed) field solve

subroutine cgyro_field_c

  use mpi
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is, ie, ix, ir
  complex :: fac
  complex, dimension(nc) :: tmp
  integer :: i,j
  logical, parameter :: use_dgemv = .false.
  real, dimension(n_theta) :: pvec_inr,pvec_ini

  call timer_lib_in('field_h')

!$omp workshare
  field_loc(:,:) = (0.0,0.0)
!$omp end workshare

  ! Poisson and Ampere RHS integrals of h

!$omp parallel private(ic,iv_loc,is,ix,ie,fac)
!$omp do reduction(+:field_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        fac = w_e(ie)*w_xi(ix)*z(is)*dens(is)*h_x(ic,iv_loc)
        field_loc(:,ic) = field_loc(:,ic)+jvec_c(:,ic,iv_loc)*fac
     enddo
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLREDUCE(field_loc(:,:),&
       field(:,:),&
       size(field(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (n_field > 2) then
     field(3,:) = field(3,:)*fcoef(3,:)
  endif

  ! Poisson LHS factors

  if (n == 0 .and. ae_flag == 1) then
!$omp parallel do &
!$omp& private(ir,i,j) &
!$omp& private(pvec_in,pvec_inr,pvec_ini,pvec_outr,pvec_outi)
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_flag == 0) then
           field(1,ic_c(ir,:)) = 0.0
        else
          if (use_dgemv) then
           pvec_in(:) = real(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
           pvec_in(:) = aimag(field(1,ic_c(ir,:)))
           call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
           field(1,ic_c(ir,:)) = pvec_outr(:) + i_c * pvec_outi(:)
          else
            do i=1,n_theta
              pvec_outr(i) = 0
              pvec_outi(i) = 0
            enddo
            do i=1,n_theta
              pvec_inr(i) = real(field(1,ic_c(ir,i)),kind=kind(xzf))
              pvec_ini(i) = aimag(field(1,ic_c(ir,i)))
            enddo

            do j=1,n_theta
            do i=1,n_theta
              pvec_outr(i) = pvec_outr(i) + xzf(ir,i,j)*pvec_inr(j)
              pvec_outi(i) = pvec_outi(i) + xzf(ir,i,j)*pvec_ini(j)
            enddo
            enddo

            do i=1,n_theta
               field(1,ic_c(ir,i)) = cmplx( pvec_outr(i),pvec_outi(i))
            enddo

          endif
        endif
     enddo

  else

     if (n_field > 2) then
        tmp(:) = field(1,:)
        field(1,:) = gcoef(1,:)*field(1,:)+gcoef(4,:)*field(3,:)
        field(2,:) = gcoef(2,:)*field(2,:)
        field(3,:) = gcoef(3,:)*field(3,:)+gcoef(5,:)*tmp(:)
     else
!$omp workshare
        field(:,:) = gcoef(:,:)*field(:,:)
!$omp end workshare
     endif

  endif

!$omp parallel do  &
!$omp& private(iv,iv_loc,is,ix,ie,ic)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        cap_h_c(ic,iv_loc) = h_x(ic,iv_loc)+psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  call timer_lib_out('field_h')

end subroutine cgyro_field_c
