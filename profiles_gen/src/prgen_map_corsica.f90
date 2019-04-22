!------------------------------------------------------------
! prgen_map_corsica.f90
!
! PURPOSE:
!  Map native corsica data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_corsica

  use prgen_globals
  use expro
  
  implicit none

  ! Compute rho, bref and arho:
  ! Ignore corsica input "rho", "rmin", and "q"; use gato and corsica poloidal
  ! flux
  call prgen_get_chi(nx,q,kappa,rmin,dpsi,rho,corsica_torfluxa)

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  expro_n_exp = corsica_nvals
  expro_n_ion = 2
  call expro_init(1)
  !
  expro_rho(:)       = rho(:)
  expro_rmin(:)      = rmin(:)
  expro_rmaj(:)      = rmaj(:)
  expro_q(:)         = q(:)
  expro_kappa(:)     = kappa(:)
  expro_delta(:)     = delta(:)
  expro_te(:)        = corsica_te(:)
  expro_ne(:)        = corsica_ne(:)*10.
  expro_z_eff(:)     = corsica_zeff(:)
  expro_w0(:)        = 0.0      ! omega
  expro_flow_mom(:)  = 0.0      ! flow_mom
  expro_pow_e(:)     = 0.0      ! pow_e
  expro_pow_i(:)     = 0.0      ! pow_i 
  expro_pow_ei(:)    = 0.0      ! pow_ei_exp
  expro_zeta(:)      = zeta(:)
  expro_flow_beam(:) = 0.0      ! flow_beam
  expro_flow_wall(:) = 0.0      ! flow_wall_exp
  expro_zmag(:)      = zmag(:)  
  expro_ptot(:)      = p_tot(:)
  expro_polflux(:)   = dpsi(:)

  ! Construct ion densities and temperatures assuming corsica ion species
  ! (n D+T) is two species, each with 1/2 n_corsica and same temperature

  expro_mass(:) = 2.0
  expro_z(:) = 1.0
  expro_type(:) = type_therm
  expro_name(:) = 'D'
  
  ! ni
  expro_ni(1,:) = 0.5 * corsica_ndt(:)*10.0
  expro_ni(2,:) = 0.5 * corsica_ndt(:)*10.0
  
  ! ti
  expro_ti(1,:) = corsica_ti(:)
  expro_ti(2,:) = corsica_ti(:)

  ! vphi
  expro_vtor(:,:) = 0.0

  ! vpol
  expro_vpol(:,:) = 0.0
  
end subroutine prgen_map_corsica
