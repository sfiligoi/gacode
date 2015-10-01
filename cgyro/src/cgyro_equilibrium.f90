subroutine cgyro_equilibrium

  use cgyro_globals
  use GEO_interface

  implicit none

  integer :: it,ir,is
  real :: ga,ga0,err
  real, dimension(n_theta+1) :: x,y

  real, parameter :: tol=1e-14

  ! Parameters needed for equilibrium
  ! geo_numeq_flag, geo_ny, and geo_yin already set 

  GEO_signb_in     = -btccw
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

  !-----------------------------------------------------------------
  ! Generate theta-grid (equally-spaced or constant-wind-speed)
  !
  d_theta = (2*pi/n_theta)
  do it=1,n_theta+1
     y(it) = -pi+(it-1)*d_theta
     x(it) = y(it)
  enddo

  if (constant_wind_flag == 1) then

     call GEO_interp(0.0)
     ga = GEO_g_theta
     err = 1e4

     do while (err > tol)
        do it=1,n_theta
           call GEO_interp(0.5*(y(it)+y(it+1)))
           x(it+1) = x(it)+d_theta/GEO_g_theta
        enddo
        ga0 = ga
        ga  = (2*pi)/(x(n_theta+1)-x(1))
        y   = (x-x(1))*ga-pi
        err = abs((ga-ga0)/ga)
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
        if(ipccw*btccw > 0) then
           thetab(ir,it) = theta(it)+2*pi*(ir-1-n_radial/2/box_size)
        else
           thetab(ir,it) = -theta(it)+2*pi*(-ir+n_radial/2/box_size)
        endif
     enddo
  enddo

  do it=1,n_theta

     call GEO_interp(theta(it))     

     do is=1,n_species

        if (constant_wind_flag == 0) then
           omega_stream(it,is) = sqrt(2.0) * vth(is) / (q * rmaj * GEO_g_theta)
        else
           omega_stream(it,is) = sqrt(2.0) * vth(is) / (q * rmaj * ga)
        endif

        omega_trap(it,is) = -0.5*sqrt(2.0) * vth(is) &
             * (GEO_dbdt / GEO_b) / (q * rmaj * GEO_g_theta) 

        omega_rdrift(it,is) = -rho * vth(is)**2 * mass(is)/(Z(is)*GEO_b) &
             * GEO_grad_r / rmaj * GEO_gsin 

        omega_adrift(it,is) = -rho * vth(is)**2 * mass(is)/(Z(is)*GEO_b) &
             * GEO_gq / rmaj &
             * (GEO_gcos1 + GEO_gcos2 + GEO_captheta * GEO_gsin) 

        omega_aprdrift(it,is) = 2.0 * rho * vth(is)**2 &
             * mass(is)/(Z(is)*GEO_b) * GEO_gq / rmaj * GEO_gcos2

        omega_cdrift(it,is) = -2.0 * sqrt(2.0) * rho * vth(is) &
             * mass(is)/(Z(is)*GEO_b) * GEO_gq / rmaj &
             * (GEO_ucos + GEO_captheta * GEO_usin) * mach

     enddo

     omega_gammap(it) = GEO_bt/GEO_b * GEO_bigr/rmaj * gamma_p

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


