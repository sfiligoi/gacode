subroutine cgyro_equilibrium

  use cgyro_globals
  use GEO_interface

  implicit none

  integer :: it,ir,is
  real :: gtheta_ave,gtheta0,err
  real, dimension(n_theta+1) :: x,y

  real, parameter :: tol=1e-14

  ! Parameters needed for equilibrium
  ! geo_numeq_flag, geo_ny, and geo_yin already set 

  GEO_signb_in     = -btccw
  GEO_rmin_in      = rmin
  GEO_rmaj_in      = rmaj
  GEO_drmaj_in     = shift
  GEO_zmag_in      = zmag
  GEO_dzmag_in     = dzmag
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

  !-----------------------------------------------------------------
  ! Generate theta-grid (equally-spaced or constant-wind-speed)
  !
  d_theta = (2*pi/n_theta)
  do it=1,n_theta+1
     y(it) = -pi+(it-1)*d_theta
     x(it) = y(it)
  enddo

  call GEO_interp(0.0)
  bigR_th0   = GEO_bigr
  bigR_r_th0 = GEO_bigr_r
  
  if (constant_stream_flag == 1) then

     gtheta_ave = GEO_g_theta  ! at theta=0
     err = 1e4

     do while (err > tol)
        do it=1,n_theta
           call GEO_interp(0.5*(y(it)+y(it+1)))
           x(it+1) = x(it)+d_theta/GEO_g_theta
        enddo
        gtheta0 = gtheta_ave
        gtheta_ave  = (2*pi)/(x(n_theta+1)-x(1))
        y   = (x-x(1))*gtheta_ave-pi
        err = abs((gtheta_ave-gtheta0)/gtheta_ave)
     enddo

  endif

  theta(:) = y(1:n_theta)

  ! Theta location of field output:
  if (zf_test_flag == 0) then
  ! Location of theta=0
     it0 = n_theta/2+1
  else
     it0 = n_theta/3+1
  endif
  !-----------------------------------------------------------------

  do ir=1,n_radial/box_size
     do it=1,n_theta
        if (ipccw*btccw > 0) then
           thetab(ir,it) = theta(it)+2*pi*(ir-1-n_radial/2/box_size)
        else
           thetab(ir,it) = -theta(it)+2*pi*(-ir+n_radial/2/box_size)
        endif
     enddo
  enddo
  
  do it=1,n_theta

     call GEO_interp(theta(it))     

     bigR(it)   = GEO_bigr
     
     do is=1,n_species

        if (constant_stream_flag == 0) then
           g_theta(it) = GEO_g_theta
        else
           g_theta(it) = gtheta_ave
        endif
        omega_stream(it,is) = sqrt(2.0)*vth(is)/(q*rmaj*g_theta(it))
 
        omega_trap(it,is) = -0.5*sqrt(2.0)*vth(is) &
             *(GEO_dbdt/GEO_b)/(q*rmaj*GEO_g_theta) 

        omega_rdrift(it,is) = -rho*vth(is)**2*mass(is)/(Z(is)*GEO_b) &
             *GEO_grad_r/rmaj*GEO_gsin 

        omega_adrift(it,is) = -rho*vth(is)**2*mass(is)/(Z(is)*GEO_b) &
             *GEO_gq/rmaj*(GEO_gcos1+GEO_gcos2+GEO_captheta*GEO_gsin) 

        omega_aprdrift(it,is) = 2.0*rho*vth(is)**2 &
             *mass(is)/(Z(is)*GEO_b)*GEO_gq/rmaj*GEO_gcos2

        ! Finite-Mach drift (Coriolis)
        omega_cdrift(it,is) = -2.0*sqrt(2.0)*rho*vth(is) &
             * mass(is)/(Z(is)*GEO_b)*GEO_gq/rmaj &
             *(GEO_ucos+GEO_captheta*GEO_usin)*mach

        omega_cdrift_r(it,is) = -2.0*sqrt(2.0)*rho*vth(is) &
             * mass(is)/(Z(is)*GEO_b)*GEO_grad_r/rmaj*GEO_usin*mach

        ! Partial Finite-Mach centrifugal terms
        ! These are re-set to 0 in cgyro_init_rot if cf_model=0

        ! bhat dot grad lambda
        ! Add phi_rot term in cgyro_init_rotation
        dlambda_rot(it,is) = - (mach / rmaj /vth(is))**2 * GEO_bigr &
             * GEO_bigr_t / (q*rmaj*GEO_g_theta)

        omega_rot_drift(it,is) = -(mach/rmaj)**2 * GEO_bigr * rho &
             * mass(is)/(z(is)*GEO_b) * GEO_gq &
             * (GEO_captheta*GEO_usin*GEO_bt/GEO_b + GEO_ucos*GEO_b/GEO_bt)
        
        omega_rot_drift_r(it,is) = -(mach/rmaj)**2 * GEO_bigr * rho &
             * mass(is)/(z(is)*GEO_b) * GEO_usin * GEO_grad_r &
             * GEO_bt / GEO_b

        ! Add phi_rot term in cgyro_init_rotation
        omega_rot_star(it,is) = -dlntdr(is) * (0.5*(mach/vth(is))**2 &
             * (GEO_bigr**2 - bigR_th0**2)/rmaj**2)  &
             + mach*gamma_p/vth(is)**2 * (GEO_bigr**2 - bigR_th0**2)/rmaj**2 &
             + (mach/vth(is))**2 * bigR_th0/rmaj**2 * bigR_r_th0 
        
        ! Multiply pressure theta derivative in cgyro_init_rotation
        omega_rot_prdrift(it,is) = -betae_unit * rho*vth(is)**2 &
             *mass(is)/(Z(is)*GEO_b) * GEO_bt/GEO_bp/GEO_b * GEO_captheta &
             / GEO_grad_r 

        ! Multiply pressure theta derivative in cgyro_init_rotation
        omega_rot_prdrift_r(it,is) = -betae_unit * rho*vth(is)**2 &
             *mass(is)/(Z(is)*GEO_b) * GEO_bt/GEO_b**2 *q*rmaj/rmin 

        ! Multiply phi_rot theta derivative in cgyro_init_rotation
        omega_rot_edrift(it,is) = -rho * GEO_bt/GEO_bp &
             * GEO_captheta / GEO_grad_r 

        ! Multiply phi_rot theta derivative in cgyro_init_rotation
        omega_rot_edrift_r(it,is) = -rho * GEO_bt/GEO_b * q * GEO_bigR/rmin 
        
     enddo

     ! Used in cgyro_init_rotation
     omega_rot_edrift_0(it) = (mach/rmaj)**2 * (GEO_bigr * GEO_bigr_r &
          - bigR_th0 * bigR_r_th0) + (mach/rmaj) * (-gamma_p/rmaj) &
          * (GEO_bigr**2 - bigR_th0**2)
     
     ! Rotation shear (GAMMA_P)
     omega_gammap(it) = GEO_bt/GEO_b*GEO_bigr/rmaj*gamma_p

     bmag(it) = GEO_b

     ! flux-surface average weights
     w_theta(it) = g_theta(it)/GEO_b

     do ir=1,n_radial
        k_perp(ic_c(ir,it)) = sqrt((2.0*pi*px(ir)*GEO_grad_r/length &
             + k_theta*GEO_gq*GEO_captheta)**2 &
             + (k_theta*GEO_gq)**2) 
     enddo

  enddo
!$acc enter data copyin(energy,xi,vel,omega_stream)

  w_theta(:) = w_theta(:)/sum(w_theta) 
  
  call cgyro_init_rotation
  
end subroutine cgyro_equilibrium


