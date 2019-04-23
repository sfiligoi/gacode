!------------------------------------------------------------
! prgen_map_inputprofiles.f90
!
! PURPOSE:
!  Map input profiles plus gfile data.
!------------------------------------------------------------

subroutine prgen_map_inputprofiles

  use prgen_globals
  use expro

  implicit none

  integer :: i
  real :: xnew(nx)
  
  if (efit_method > 1) then

     ! Compute rho, bref and arho based on GATO dpsi and q.
     call prgen_get_chi(nx,q,kappa,rmin,dpsi,rho,expro_torfluxa)

     ! Align with new (GATO) rho grid, not old rho grid, expro_rho:
     call cub_spline(expro_rho,expro_te,nx,rho,xnew,nx)
     expro_te = xnew
     call cub_spline(expro_rho,expro_ne,nx,rho,xnew,nx)
     expro_ne = xnew
     call cub_spline(expro_rho,expro_z_eff,nx,rho,xnew,nx)
     expro_z_eff = xnew
     call cub_spline(expro_rho,expro_w0,nx,rho,xnew,nx)
     expro_w0 = xnew
     call cub_spline(expro_rho,expro_flow_mom,nx,rho,xnew,nx)
     expro_flow_mom = xnew
     call cub_spline(expro_rho,expro_pow_e,nx,rho,xnew,nx)
     expro_pow_e = xnew
     call cub_spline(expro_rho,expro_pow_i,nx,rho,xnew,nx)
     expro_pow_i = xnew
     call cub_spline(expro_rho,expro_pow_ei,nx,rho,xnew,nx)
     expro_pow_ei = xnew
     call cub_spline(expro_rho,expro_flow_beam,nx,rho,xnew,nx)
     expro_flow_beam = xnew
     call cub_spline(expro_rho,expro_flow_wall,nx,rho,xnew,nx)
     expro_flow_wall = xnew
     call cub_spline(expro_rho,expro_ptot,nx,rho,xnew,nx)
     expro_ptot = xnew
     do i=1,10
        call cub_spline(expro_rho,expro_ni(i,:),nx,rho,xnew,nx)
        expro_ni(i,:) = xnew(:)
        call cub_spline(expro_rho,expro_ti(i,:),nx,rho,xnew,nx)
        expro_ti(i,:) = xnew(:)
        call cub_spline(expro_rho,expro_vtor(i,:),nx,rho,xnew,nx)
        expro_vtor(i,:) = xnew(:)
        call cub_spline(expro_rho,expro_vpol(i,:),nx,rho,xnew,nx)
        expro_vpol(i,:) = xnew(:)
     enddo

     !---------------------------------------------------------
     ! Map profile data into expro interface.
     ! NOTE: expro_alloc already called in prgen_read_inputprofiles
     !
     expro_rho(:)  = rho(:)
     expro_rmin(:) = rmin
     expro_rmaj(:) = rmaj(:)
     ! COORDINATES: set sign of q
     expro_q(:)     = abs(q(:))*ipccw*btccw
     expro_kappa(:) = kappa(:)
     expro_delta(:) = delta(:)
     expro_zeta(:)  = zeta(:)
     expro_zmag(:)  = zmag(:)
     ! COORDINATES: set sign of poloidal flux
     expro_polflux = abs(dpsi(:))*(-ipccw)

  endif

  !---------------------------------------------------------
  ! Notification about CER file
  !
  if (file_cer /= "null") then
    print '(a)','INFO: (prgen_map_inputprofiles) IGNORING cer file'
  endif
  !---------------------------------------------------------

end subroutine prgen_map_inputprofiles

