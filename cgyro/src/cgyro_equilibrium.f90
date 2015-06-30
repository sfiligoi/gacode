subroutine cgyro_equilibrium
  
  use cgyro_globals
  use GEO_interface
  
  implicit none
  
  integer :: it, ir, is
  
  ! parameters needed for Miller equilibrium
  ! geo_numeq_flag, geo_ny, and geo_yin already set    
  GEO_rmin_in      = rmin
  GEO_rmaj_in      = rmaj
  GEO_drmaj_in     = shift
  GEO_zmag_in      = zmag
  GEO_dzmag_in     = s_zmag
  GEO_q_in         = q
  GEO_s_in         = s
  GEO_kappa_in     = kappa
  GEO_s_kappa_in   = s_kappa
  GEO_delta_in     = delta
  GEO_s_delta_in   = s_delta
  GEO_zeta_in      = zeta
  GEO_s_zeta_in    = s_zeta
  GEO_beta_star_in = beta_star
  GEO_fourier_in(:,:) = geo_yin(:,:)
  
  call GEO_do()  
  
  do it=1,n_theta
     
     call GEO_interp(theta(it))
     
     do is=1,n_species
        
        omega_stream(it,is) = sqrt(2.0) * vth(is) / (q * rmaj * GEO_g_theta)
        
        omega_trap(it,is) = -0.5*sqrt(2.0) * vth(is) &
             * (GEO_dbdt / GEO_b) / (q * rmaj * GEO_g_theta) 
        
        omega_rdrift(it,is) = -rho * vth(is)**2 * mass(is)/(Z(is)*GEO_b) &
             * GEO_grad_r / rmaj * GEO_gsin 
        
        omega_adrift(it,is) = -rho * vth(is)**2 * mass(is)/(Z(is)*GEO_b) &
             * GEO_gq / rmaj &
             * (GEO_gcos1 + GEO_gcos2 + GEO_captheta * GEO_gsin) 
        
        omega_aprdrift(it,is) = 2.0 * rho * vth(is)**2 &
             * mass(is)/(Z(is)*GEO_b) * GEO_gq / rmaj * GEO_gcos2
        
     enddo
     
     Bmag(it) = GEO_b
     
     ! flux-surface average weights
     w_theta(it) = GEO_g_theta / GEO_b
     
     do ir=1,n_radial
        k_perp(it,ir) = sqrt((2.0*pi*px(ir)*GEO_grad_r/length &
             + k_theta*GEO_gq*GEO_captheta)**2 &
             + (k_theta*GEO_gq)**2) 
     enddo
     
  enddo
  
  w_theta(:) = w_theta(:)/sum(w_theta) 
  
end subroutine cgyro_equilibrium


