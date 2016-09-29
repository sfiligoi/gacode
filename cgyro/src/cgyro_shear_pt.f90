!---------------------------------------------------------
! cgyro_shear_pt.f90
!
! PURPOSE:
!  Spectral shear algorithm with linear wavenumber variation.

! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!---------------------------------------------------------

subroutine cgyro_shear_pt

  use cgyro_shear_interface

  use cgyro_globals
  use timer_lib
  use GEO_interface
  use parallel_lib
  use mpi

  implicit none

  integer :: ir,it,ix,ie,is
  real :: a,fac
  complex, dimension(n_theta,nv_loc) :: a1

  if (i_time == 1) then
     allocate(fcoef0(n_field,nc))
     allocate(gcoef0(n_field,nc))
     allocate(omega_s0(n_field,nc,nv_loc))
     allocate(omega_cap_h0(nc,nv_loc))
     fcoef0     = fcoef
     gcoef0     = gcoef
     omega_s0   = omega_s
     omega_cap_h0 = omega_cap_h
  endif

  a     = omega_eb*delta_t
  gtime = gtime+a

  ! Forward shearing
  if (gtime > 1.0) then

     gtime = gtime-1.0

     a1 = h_x(ic_c(1,:),:)

     do ir=2,n_radial
        h_x(ic_c(ir-1,:),:) = h_x(ic_c(ir,:),:)
     enddo

     h_x(ic_c(n_radial,:),:) = a1*gamma_e_decay

     call cgyro_field_c

  endif

  do ic=1,nc
     ir = ir_c(ic) 
     if (ir < n_radial) then 
        it = it_c(ic)
!        fcoef(:,ic) = fcoef0(:,ic_c(ir,it))*(1-gtime) + fcoef0(:,ic_c(ir+1,it))*gtime
        gcoef(:,ic) = gcoef0(:,ic_c(ir,it))*(1-gtime) + gcoef0(:,ic_c(ir+1,it))*gtime
        omega_cap_h(ic,:) = omega_cap_h0(ic_c(ir,it),:)*(1-gtime) &
             +omega_cap_h0(ic_c(ir+1,it),:)*gtime
        omega_s(:,ic,:) = omega_s0(:,ic_c(ir,it),:)*(1-gtime) &
             +omega_s0(:,ic_c(ir+1,it),:)*gtime
        jvec_c(:,ic,:) = jvec_c0(:,ic_c(ir,it),:)*(1-gtime) &
             +jvec_c0(:,ic_c(ir+1,it),:)*gtime
     endif
  enddo

  call parallel_lib_rtrans_real(jvec_c(1,:,:),jvec_v(1,:,:))

  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        fac = w_xi(ix)*w_e(ie)*z(is)**2/temp(is)*dens(is)
        field_loc(1,ic) = field_loc(1,ic)+(1-jvec_c(1,ic,iv_loc)**2)*fac
     enddo
  enddo

  call MPI_ALLREDUCE(field_loc(1,:),&
       gcoef(1,:),&
       size(field(1,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  gcoef = 1.0/gcoef

  call cgyro_field_c

end subroutine cgyro_shear_pt
