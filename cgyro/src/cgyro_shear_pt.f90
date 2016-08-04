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

  implicit none

  integer :: ir,it
  real :: a
  complex, dimension(n_theta,nv_loc) :: a1

  if (i_time == 1) then
     allocate(fcoef0(n_field,nc))
     allocate(gcoef0(n_field,nc))
     allocate(omega_s0(n_field,nc,nv_loc))
     allocate(omega_cap_h0(nc,nv_loc))
     allocate(jvec_c0(n_field,nc,nv_loc))
     fcoef0     = fcoef
     gcoef0     = gcoef
     omega_s0   = omega_s
     jvec_c0    = jvec_c
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
        fcoef(:,ic) = fcoef0(:,ic_c(ir,it))*(1-gtime) + fcoef0(:,ic_c(ir+1,it))*gtime
        gcoef(:,ic) = gcoef0(:,ic_c(ir,it))*(1-gtime) + gcoef0(:,ic_c(ir+1,it))*gtime
        omega_cap_h(ic,:) = omega_cap_h0(ic_c(ir,it),:)*(1-gtime) &
             +omega_cap_h0(ic_c(ir+1,it),:)*gtime
        omega_s(:,ic,:) = omega_s0(:,ic_c(ir,it),:)*(1-gtime) &
             +omega_s0(:,ic_c(ir+1,it),:)*gtime
        jvec_c(:,ic,:) = jvec_c0(:,ic_c(ir,it),:)*(1-gtime) &
             +jvec_c0(:,ic_c(ir+1,it),:)*gtime
     endif
  enddo

  call cgyro_field_c

end subroutine cgyro_shear_pt
