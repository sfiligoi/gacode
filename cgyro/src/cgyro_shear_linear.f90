!-----------------------------------------------------------------
! cgyro_shear_linear.f90
!
! PURPOSE:
!  Advanced spectral shear algorithm with linear wavenumber 
!  variation plus wavenumber shift.  This should be more 
!  accurate than the Hammett method.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!-----------------------------------------------------------------

subroutine cgyro_shear_linear

  use cgyro_shear_interface

  use cgyro_globals
  use parallel_lib

  implicit none

  integer :: ir,i_field
  integer :: ica,icb
  complex, dimension(n_theta,nv_loc) :: a1


  if (i_time == 1) then
     allocate(omega_s0(n_field,nc,nv_loc))
     allocate(omega_cap_h0(nc,nv_loc))
     allocate(jvec_c0(n_field,nc,nv_loc))
     allocate(k_perp0(nc))
     omega_s0     = omega_s
     omega_cap_h0 = omega_cap_h
     jvec_c0 = jvec_c
     k_perp0 = k_perp
  endif

  gtime = gtime+omega_eb*delta_t

  ! Forward shearing
  if (gtime > 1.0) then

     gtime = gtime-1.0

     a1 = h_x(ic_c(1,:),:)

     do ir=2,n_radial
        h_x(ic_c(ir-1,:),:) = h_x(ic_c(ir,:),:)
     enddo
     h_x(ic_c(n_radial,:),:) = a1*gamma_e_decay

  endif
 
!$omp parallel do private(ica,icb)
  do ic=1,nc
     ica = ica_c(ic)
     icb = icb_c(ic)
     ! Linear interpolation of kx-dependent arrays
     omega_cap_h(ic,:) = omega_cap_h0(ica,:)*(1-gtime)+omega_cap_h0(icb,:)*gtime
     omega_s(:,ic,:)   =   omega_s0(:,ica,:)*(1-gtime)  +omega_s0(:,icb,:)*gtime
     jvec_c(:,ic,:)    =    jvec_c0(:,ica,:)*(1-gtime)   +jvec_c0(:,icb,:)*gtime
     k_perp(ic)        =        k_perp0(ica)*(1-gtime)       +k_perp0(icb)*gtime
  enddo

  do i_field=1,n_field
     call parallel_lib_rtrans_real(jvec_c(i_field,:,:),jvec_v(i_field,:,:))
  enddo

  call cgyro_field_coefficients
  call cgyro_field_c

end subroutine cgyro_shear_linear
