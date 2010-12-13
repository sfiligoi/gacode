module neo_theory

  ! Contains routines to compute analytical neoclassical theory values.
  ! NOTE: The theories assume the diamagnetic ordering.  Thus, rotation
  ! effects are not included here.

  implicit none

  public :: THEORY_alloc, THEORY_do, THEORY_ZF_do
  real :: pflux_HH, efluxi_HH, efluxe_HH, jpar_HH, kpar_HH, uparB_HH, &
       vpol_ion_HH, &
       efluxi_CH1, efluxi_CH2, efluxi_TG1, efluxi_TG2, &
       jpar_S, kpar_S, uparB_S, vpol_ion_S, phi_HR
  real, dimension(:), allocatable :: pflux_multi_HS1, eflux_multi_HS1, &
       pflux_multi_HS2, eflux_multi_HS2
  
  real, private :: eps ! r/rmaj
  real, private :: EparB_avg     ! <E_|| B>
  real, private :: sigma_spitzer ! Spitzer resistivity
  real, private :: nui_HH, nue_HH, nui_star_HH, nue_star_HH
  integer, private :: is_ele, is_ion
  integer, private :: ir_global, is_global, ietype
  integer, private :: zf_type, zf_ion1, zf_ion2
  integer, parameter, private :: io=40
  logical, private :: initialized = .false.

contains

  subroutine THEORY_alloc(flag)
    use neo_globals
    use neo_nclass_dr
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: is

    if(flag == 1) then
       if(initialized) return
       if(adiabatic_ele_model == 0) then
          do is=1, n_species
             if(Z(is) == -1) then
                is_ele = is
                exit
             endif
          enddo
       endif
       ! assume primary ion species is is=1
       is_ion = 1

       allocate(pflux_multi_HS1(n_species))
       allocate(eflux_multi_HS1(n_species))
       allocate(pflux_multi_HS2(n_species))
       allocate(eflux_multi_HS2(n_species))

       if(write_out_mode > 0) then
          open(io,file='theory.out',status='replace')
       end if

       call NCLASS_DR_alloc(1)
       
       initialized = .true.

    else
       if(.NOT. initialized) return
       deallocate(pflux_multi_HS1)
       deallocate(eflux_multi_HS1)
       deallocate(pflux_multi_HS2)
       deallocate(eflux_multi_HS2)
       if(write_out_mode > 0) then
          close(io)
       end if
       call NCLASS_DR_alloc(0)
       initialized = .false.
    end if

  end subroutine THEORY_alloc

  ! Driver routine for computing theory results
  ! Theories assume no rotation -- so use input dens(psi) as dens
  ! (i.e. not with poloidal part)
  subroutine THEORY_do(ir)
    use neo_globals
    use neo_equilibrium
    use neo_nclass_dr
    implicit none
    integer, intent (in) :: ir
    integer :: is, it

    ! inverse aspect ratio
    eps = r(ir) / rmaj(ir)
    
    ! theory coll freqs and Spitzer resistivity
    nui_HH = nu(is_ion,ir) * (4.0/3.0) / sqrt(2.0*pi)
    nui_star_HH = nui_HH * rmaj(ir) * abs(q(ir)) &
         / (sqrt(eps)*sqrt(eps)*sqrt(eps) * vth(is_ion,ir))
    if(adiabatic_ele_model == 1) then
       nue_HH        = 0.0
       nue_star_HH   = 0.0
       sigma_spitzer = 0.0
    else
       nue_HH = nu(is_ion,ir) * (4.0/3.0) / sqrt(pi) &
            * sqrt(mass(is_ion)/mass(is_ele)) &
            * (temp(is_ion,ir) / temp(is_ele,ir))**1.5 &
            * 1.0 / (1.0* Z(is_ion))**2
       nue_star_HH = nue_HH * rmaj(ir) * abs(q(ir)) &
            / (sqrt(eps)*sqrt(eps)*sqrt(eps) * vth(is_ele,ir))
       sigma_spitzer = dens(is_ele,ir) / ( mass(is_ele) * nue_HH &
            * (0.29 + 0.46/(1.08+z(is_ion))) )
    endif
    
    ! inductive electric field
    ! < E_|| B>
    EparB_avg = 0.0
    do it=1, n_theta
       EparB_avg = EparB_avg + w_theta(it) * epar0(ir) * Bmag(it) 
    enddo

    ! compute the theory
    call compute_HH(ir,pflux_HH, efluxi_HH, efluxe_HH, &
         jpar_HH, kpar_HH, uparB_HH)
    call compute_vpol_ion(ir,kpar_HH, vpol_ion_HH)
    call compute_CH(ir,efluxi_CH1, 1)
    call compute_CH(ir,efluxi_CH2, 2)
    call compute_TG(ir,efluxi_TG1, 1)
    call compute_TG(ir,efluxi_TG2, 2)
    call compute_Sauter(ir,jpar_S, kpar_S, uparB_S)
    call compute_vpol_ion(ir,kpar_S, vpol_ion_S)
    call compute_HR(ir,phi_HR)
    call compute_HS(ir,pflux_multi_HS1, eflux_multi_HS1, 1)
    call compute_HS(ir,pflux_multi_HS2, eflux_multi_HS2, 2)
    
    if(write_out_mode > 0) then
       write(io,'(e16.8,$)') r(ir)
       write(io,'(e16.8,$)') pflux_HH
       write(io,'(e16.8,$)') efluxi_HH
       write(io,'(e16.8,$)') efluxe_HH
       write(io,'(e16.8,$)') jpar_HH
       write(io,'(e16.8,$)') kpar_HH
       write(io,'(e16.8,$)') uparB_HH
       write(io,'(e16.8,$)') vpol_ion_HH
       write(io,'(e16.8,$)') efluxi_CH1
       write(io,'(e16.8,$)') efluxi_CH2
       write(io,'(e16.8,$)') efluxi_TG1
       write(io,'(e16.8,$)') efluxi_TG2
       write(io,'(e16.8,$)') jpar_S
       write(io,'(e16.8,$)') kpar_S
       write(io,'(e16.8,$)') uparB_S
       write(io,'(e16.8,$)') vpol_ion_S
       write(io,'(e16.8,$)') phi_HR
       do is=1, n_species
          write(io,'(e16.8,$)') pflux_multi_HS1(is)
          write(io,'(e16.8,$)') eflux_multi_HS1(is)
       enddo
       do is=1, n_species
          write(io,'(e16.8,$)') pflux_multi_HS2(is)
          write(io,'(e16.8,$)') eflux_multi_HS2(is)
       enddo
       write (io,*)

    end if

    call NCLASS_DR_do(ir)

  end subroutine THEORY_do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Hinton-Rosenbluth potential
  ! Phys. Fluids, vol. 16, 836 (1973)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_HR(ir,phi)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: phi
    integer :: it
    real :: thavg_fac=0.5        ! theta-avg factor <sin^2> = 0.5
    real ::  nu_crit, phi_1, phi_2, fac, temp_ele

    if(adiabatic_ele_model == 1) then
       temp_ele = te_ade(ir)
    else
       temp_ele = temp(is_ele,ir)
    end if

    fac = (q(ir)/eps) * eps * dlntdr(is_ion,ir) * rho(ir) &
         * (sqrt(temp(is_ion,ir) * mass(is_ion)) / (1.0 * Z(is_ion))) &
         / (1.0*Z(is_ion) / temp(is_ion,ir) + 1.0/temp_ele)
    
    ! banana regime -- phi_1 = phi_1 * nui_HH
    phi_1 = 1.81 * nui_star_HH * fac / nui_HH
    ! plateau regime
    phi_2 = sqrt(pi/2.0) * fac
  
    nu_crit = phi_2 / phi_1

    if(nui_HH > nu_crit) then
       phi = phi_2 * phi_2 * thavg_fac
    else
       phi = (phi_1 * nui_HH) * (phi_1 * nui_HH) * thavg_fac
    endif
    
  end subroutine compute_HR
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Hinton-Hazeltine fluxes and flows
  ! Rev. Mod. Phys., vol. 48, 239 (1976)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_HH(ir,pflux, efluxi, efluxe, jpar, kpar, uparB)
    use neo_globals
    use neo_equilibrium, only : I_div_psip, Bmag2_avg, Bmag2inv_avg
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: pflux, efluxi, efluxe, jpar, kpar, uparB
    real :: k0, a0, b0, c0, fac, beta1, A1e, K11, K12, K22, K2, K13, K23, K33,&
         gi, qi

    fac = (1.17-0.35*sqrt(nui_star_HH)) / (1+0.7*sqrt(nui_star_HH))
    beta1 = (fac - 2.1*nui_star_HH*nui_star_HH*eps*eps*eps) &
         / (1.0 + nui_star_HH*nui_star_HH*eps*eps*eps)
    
    kpar = beta1
    
    ! Note: only jpar and upar B use I_div_psip rather than q/eps
    uparB  = I_div_psip * rho(ir) * temp(is_ion,ir) / (z(is_ion)*1.0) &
         * (dlnndr(is_ion,ir) &
         - (z(is_ion)*1.0)/temp(is_ion,ir) * dphi0dr(ir) &
         + (1.0 - kpar) * dlntdr(is_ion,ir))

    k0=0.66; a0=1.03; b0=0.31; c0=0.74
    K2 = k0 * ( 1.0/(1 + a0*sqrt(nui_star_HH) + b0*nui_star_HH) &
         + eps*eps*eps * (c0*c0/b0) * nui_star_HH &
         / (1 + c0 * sqrt(eps*eps*eps) * nui_star_HH) )
    
    ! Note: ion Q get additional part below if not adiabatic ele
    qi = (q(ir)/eps)**2 &
         * dens(is_ion,ir) * temp(is_ion,ir) &
         * (2.0 * mass(is_ion) * temp(is_ion,ir) / (1.0 *Z(is_ion)**2)) &
         * rho(ir)**2 * sqrt(eps) * nui_HH * K2 * dlntdr(is_ion,ir)
    
    if(adiabatic_ele_model == 1) then
       pflux  = 0.0
       efluxe = 0.0
       jpar   = 0.0
       efluxi = qi
       
    else
       
       A1e = -(dlnndr(is_ele,ir) + dlntdr(is_ele,ir)) &
            + 2.5 * dlntdr(is_ele,ir) &
            + temp(is_ion,ir) / (temp(is_ele,ir) * Z(is_ion)) &
              * ( -(dlnndr(is_ion,ir) + dlntdr(is_ion,ir)) &
                  + beta1 / (1.0 + nue_star_HH*nue_star_HH*eps*eps) &
                    * dlntdr(is_ion,ir) )

       k0=1.04; a0=2.01; b0=1.53; c0=0.89
       K11 = k0 * ( 1.0/(1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) &
            + eps*eps*eps * (c0*c0/b0) * nue_star_HH  &
            / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )
       
       k0=1.20; a0=0.76; b0=0.67; c0=0.56
       K12 =  k0 * ( 1.0/(1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) &
            + eps*eps*eps * (c0*c0/b0) * nue_star_HH &
            / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )

       k0=2.55; a0=0.45; b0=0.43; c0=0.43
       K22 = k0 * ( 1.0/(1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) &
            + eps*eps*eps * (c0*c0/b0) * nue_star_HH &
            / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )

       k0=2.30; a0=1.02; b0=1.07; c0=1.07
       K13 = k0 * ( 1.0 / (1 + a0*sqrt(nue_star_HH) + b0* nue_star_HH) ) &
            * ( 1.0 / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )

       k0=4.19; a0=0.57; b0=0.61; c0=0.61
       K23 = k0 * ( 1.0 / (1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) ) &
            * ( 1.0 / (1 + c0  * sqrt(eps*eps*eps) * nue_star_HH) )
       
       k0=1.83; a0=0.68; b0=0.32; c0=0.66
       K33 = k0 * ( 1.0 / (1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) ) &
            * ( 1.0 / (1 + c0  * sqrt(eps*eps*eps) * nue_star_HH) )

       gi  = -(q(ir)/eps)**2 * dens(is_ele,ir) &
            * (2.0*mass(is_ele)*temp(is_ele,ir)) *rho(ir)**2 &
            * sqrt(eps) * nue_HH * (K11*A1e - K12*dlntdr(is_ele,ir)) &
            - EparB_avg  * dens(is_ele,ir) * sqrt(eps) * rho(ir) * I_div_psip &
            * K13

       pflux  = gi

       efluxi = qi &
            - beta1 * temp(is_ion,ir) / (1.0 *Z(is_ion)) * gi &
            / (1.0 + nue_star_HH*nue_star_HH*eps*eps*eps) &
            + 2.5 * temp(is_ion,ir) * gi

       efluxe = -(q(ir)/eps)**2 &
            * dens(is_ele,ir) * temp(is_ele,ir) &
            * (2.0*mass(is_ele)*temp(is_ele,ir)) *rho(ir)**2 &
            * sqrt(eps) * nue_HH * (K12*A1e - K22*dlntdr(is_ele,ir)) &
            - EparB_avg  * dens(is_ele,ir) * sqrt(eps) * rho(ir) * I_div_psip &
            * temp(is_ele,ir) * K23
       

       ! Note: only jpar and upar B and inductive E-field 
       ! use I_div_psip rather than q/eps

       jpar   = -(I_div_psip) * dens(is_ele,ir) * temp(is_ele,ir) &
            * sqrt(eps) * rho(ir) * (K13 * A1e - K23*dlntdr(is_ele,ir)) &
            - EparB_avg * sigma_spitzer * sqrt(eps) * K33 &
            + EparB_avg * sigma_spitzer * Bmag2_avg * Bmag2inv_avg

    end if

  end subroutine compute_HH


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Chang-Hinton ion heat flux (q_i) (Q = q_i + 2.5 * particle flux)
  ! Phys. Plasmas, vol. 25, 1493 (1982)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_CH(ir,efluxi, geo_type)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, Bmag2_avg, Bmag2inv_avg
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: efluxi
    integer, intent (in) :: geo_type
    real :: k0, a0, b0, c0, K2, K2_star, F2, &
         CH_Bmag2inv_avg, CH_Bmag2avg_inv, CH_I_div_psip

    k0=0.66; a0=1.03; b0=0.31; c0=0.74

    if(geo_type == 1) then
       ! original Chang-Hinton
       CH_Bmag2inv_avg = (1.0 + 1.5*(eps*eps + eps*shift(ir)) &
            + 0.375 * eps*eps*eps * shift(ir)) &
            / (1.0 + 0.5*eps*shift(ir))
       CH_Bmag2avg_inv = (sqrt(1.0-eps*eps) * (1.0 + 0.5*eps*shift(ir))) &
            / (1.0 + (shift(ir)/eps) * (sqrt(1.0 - eps*eps) - 1.0))
       CH_I_div_psip   = q(ir)/eps

    else if(geo_type == 2) then
       ! modified Chang-Hinton (for shaped plasmas)
       CH_Bmag2inv_avg = Bmag2inv_avg
       CH_Bmag2avg_inv = 1.0 / Bmag2_avg
       CH_I_div_psip   = I_div_psip
    endif

    F2 = (0.5/sqrt(eps)) * (CH_Bmag2inv_avg - CH_Bmag2avg_inv)

    K2_star = (k0 + 1.88*sqrt(eps) - 1.54*eps) * CH_Bmag2inv_avg

    K2 = k0 *( (K2_star/k0) / (1 + a0*sqrt(nui_star_HH) + b0*nui_star_HH) &
         + sqrt(eps*eps*eps) * (c0*c0/b0) * nui_star_HH * F2 &
         / (1.0 + c0 * sqrt(eps*eps*eps) * nui_star_HH) )

    efluxi = CH_I_div_psip**2 &
         * dens(is_ion,ir) * temp(is_ion,ir) &
         * (2.0*mass(is_ion)*temp(is_ion,ir) / (1.0 *Z(is_ion)**2)) &
         * rho(ir)**2 * sqrt(eps) * nui_HH * K2 * dlntdr(is_ion,ir)
    
  end subroutine compute_CH

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Taguchi ion heat flux (q_i) (Q = q_i + 2.5 * particle flux)
  ! modified with Chang-Hinton collisional interpolation factor
  ! PPCF, vol. 30, 1897 (1988)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_TG(ir,efluxi, geo_type)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, Bmag2_avg, Bmag2inv_avg, ftrap
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: efluxi
    integer, intent (in) :: geo_type
    integer :: it
    real :: k0, a0, b0, c0, K2, K2_star, F2, &
         CH_Bmag2inv_avg, CH_Bmag2avg_inv, CH_I_div_psip, fc

    ! Geo factors     
    if(geo_type == 1) then
       ! original Chang-Hinton
       CH_Bmag2inv_avg = (1.0 + 1.5*(eps*eps + eps*shift(ir)) &
            + 0.375 * eps*eps*eps * shift(ir)) &
            / (1.0 + 0.5*eps*shift(ir))
       CH_Bmag2avg_inv = (sqrt(1.0-eps*eps) * (1.0 + 0.5*eps*shift(ir))) &
            / (1.0 + (shift(ir)/eps) * (sqrt(1.0 - eps*eps) - 1.0))
       CH_I_div_psip   = q(ir)/eps

    else if(geo_type == 2) then
       ! modified Chang-Hinton (for shaped plasmas)
       CH_Bmag2inv_avg = Bmag2inv_avg
       CH_Bmag2avg_inv = 1.0 / Bmag2_avg
       CH_I_div_psip   = I_div_psip
    endif

    ! Taguchi thermal conductivity
    fc = 1.0 - ftrap

    K2_star = (CH_Bmag2inv_avg - CH_Bmag2avg_inv) &
         + CH_Bmag2avg_inv * (1.0-fc)/(1.0+1.17*fc)

    ! Chang-Hinton collisional interpolation factor
    k0=0.66; a0=1.03; b0=0.31; c0=0.74

    F2 = (0.5/sqrt(eps)) * (CH_Bmag2inv_avg - CH_Bmag2avg_inv)
    
    K2 = k0 *( (K2_star/k0) / (1 + a0*sqrt(nui_star_HH) + b0*nui_star_HH) &
         + sqrt(eps*eps*eps) * (c0*c0/b0) * nui_star_HH * F2 &
         / (1.0 + c0 * sqrt(eps*eps*eps) * nui_star_HH) )

    ! sqrt(eps) factor is already in K2_star
    efluxi = CH_I_div_psip**2 &
         * dens(is_ion,ir) * temp(is_ion,ir) &
         *   (2.0*mass(is_ion)*temp(is_ion,ir) / (1.0 *Z(is_ion)**2)) &
         * rho(ir)**2 * nui_HH * K2 * dlntdr(is_ion,ir)
    
  end subroutine compute_TG

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Sauter bootstrap current model
  ! Phys. Plasmas, vol. 6, 2834 (1999)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_Sauter(ir,jpar, kpar, uparB)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, ftrap
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: jpar, kpar, uparB
    real :: X31, L31_S, X32e, F32_ee, X32i, F32_ei, L32_S, alpha_0, alpha_S, &
         X34, L34_S, Rpe, X33, sigma_S

    alpha_0 = -1.17*(1.0-ftrap) / (1.0 - 0.22*ftrap - 0.19*ftrap*ftrap)
    
    alpha_S = ( (alpha_0 + 0.25*(1.0-ftrap*ftrap) * sqrt(nui_star_HH)) &
         / (1.0 + 0.5 * sqrt(nui_star_HH)) &
         + (0.315) * nui_star_HH * nui_star_HH  &
         * ftrap*ftrap*ftrap*ftrap*ftrap*ftrap ) &
         / (1.0 + 0.15 *nui_star_HH*nui_star_HH &
         * ftrap*ftrap*ftrap*ftrap*ftrap*ftrap)
    
    kpar = -1.0 * alpha_S

    uparB  = I_div_psip * rho(ir) * temp(is_ion,ir) / (z(is_ion)*1.0) &
         * (dlnndr(is_ion,ir) &
         - (z(is_ion)*1.0)/temp(is_ion,ir) * dphi0dr(ir) &
         + (1.0 - kpar) * dlntdr(is_ion,ir))

    if(adiabatic_ele_model == 1) then
       jpar = 0.0

    else

       X31    = ftrap / (1.0 + (1.0-0.1*ftrap) * sqrt(nue_star_HH)  &
            + 0.5*(1.0-ftrap) * nue_star_HH / (1.0*Z(is_ion)))
       
       L31_S  = (1.0 + 1.4/(Z(is_ion)+1)) * X31 &
            - (1.9/(Z(is_ion)+1)) * X31*X31 &
            + (0.3/(Z(is_ion)+1)) * X31*X31*X31 &
            + (0.2/(Z(is_ion)+1)) * X31*X31*X31*X31
       
       X32e   = ftrap / (1.0 + 0.26*(1-ftrap) * sqrt(nue_star_HH) &
            + 0.18*(1.0-0.37*ftrap) * nue_star_HH / sqrt(1.0*Z(is_ion)))
       
       F32_ee = (0.05 + 0.62*Z(is_ion)) &
            / (Z(is_ion)*(1+0.44*Z(is_ion))) * (X32e - X32e*X32e*X32e*X32e) &
            + 1.0/(1+0.22*Z(is_ion)) * (X32e*X32e - X32e*X32e*X32e*X32e &
            - 1.2*(X32e*X32e*X32e - X32e*X32e*X32e*X32e)) &
            + 1.2/(1+0.5*Z(is_ion)) * X32e*X32e*X32e*X32e
       
       X32i   = ftrap / (1.0 + (1+0.6*ftrap) * sqrt(nue_star_HH)  &
            + 0.85*(1-0.37*ftrap) * nue_star_HH*(1.0+Z(is_ion)))
       
       F32_ei = -(0.56 + 1.93*Z(is_ion)) & 
            /(Z(is_ion)*(1+0.44*Z(is_ion))) * (X32i - X32i*X32i*X32i*X32i)  &
            + 4.95/(1+2.48*Z(is_ion)) * (X32i*X32i - X32i*X32i*X32i*X32i &
            - 0.55*(X32i*X32i*X32i - X32i*X32i*X32i*X32i)) &
            + (-1.2)/(1+0.5*Z(is_ion)) * X32i*X32i*X32i*X32i;
       
       L32_S  = F32_ee + F32_ei
       
       X34   = ftrap / (1.0 + (1.0-0.1*ftrap) * sqrt(nue_star_HH) &
            + 0.5*(1.0-0.5*ftrap) * nue_star_HH/(1.0*Z(is_ion)))
       
       L34_S = (1 + 1.4/(Z(is_ion)+1.0)) * X34 &
            - (1.9/(Z(is_ion)+1.0)) * X34*X34 &
            + (0.3/(Z(is_ion)+1.0)) * X34*X34*X34 &
            + (0.2/(Z(is_ion)+1.0)) * X34*X34*X34*X34
       
       X33 = ftrap / (1.0 + (0.55-0.1*ftrap) * sqrt(nue_star_HH) &
            + 0.45*(1.0-ftrap) * nue_star_HH/(1.0*Z(is_ion)**1.5))

       sigma_S = sigma_spitzer / (0.58+0.74/(0.76+Z(is_ion))) &
            * (1.0 - (1.0+0.36/(1.0*Z(is_ion))) * X33 &
            + 0.59/(1.0*Z(is_ion)) * X33 * X33 &
            - 0.23/(1.0*Z(is_ion)) * X33 * X33 * X33)

       Rpe = dens(is_ele,ir) * temp(is_ele,ir) &
            / (dens(is_ele,ir) * temp(is_ele,ir) & 
            + dens(is_ion,ir) * temp(is_ion,ir))
       
       jpar = I_div_psip * rho(ir) &
            * (dens(is_ele,ir) * temp(is_ele,ir) &
            + dens(is_ion,ir)*temp(is_ion,ir)) &
            * (L31_S * dlnndr(is_ele,ir) &
            + Rpe * (L31_S + L32_S) * dlntdr(is_ele,ir) &
            + (1.0 - Rpe) * (L31_S + alpha_S * L34_S) * dlntdr(is_ion,ir)) &
            + sigma_S * EparB_avg

    end if
    
  end subroutine compute_Sauter

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Hirshman-Sigmar fluxes
  ! Phys. Fluids, vol. 20, 418 (1977)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_HS(ir,pflux_multi, eflux_multi, geo_type)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, Bmag2_avg, ftrap
    implicit none
    integer, intent (in) :: ir
    real, intent (inout) :: pflux_multi(n_species), eflux_multi(n_species)
    integer, intent (in) :: geo_type
    real, dimension(:), allocatable :: nux0, nux2, nux4
    real :: A1, A2, omega_fac, HS_I_div_psip, &
         L11, L12, L21, L22, L_a, L_b, sum_nm, eii_val
    integer :: Nx=100
    integer, parameter :: integ_order = 64
    integer :: js

    if(geo_type == 1) then
       omega_fac = (sqrt(1.0-eps*eps) * (1.0 + 0.5*eps*shift(ir))) &
            / (1.0 + (shift(ir)/eps) * (sqrt(1.0 - eps*eps) - 1.0))
       HS_I_div_psip = q(ir)/eps
    else if(geo_type == 2) then
       omega_fac = 1.0 / Bmag2_avg
       HS_I_div_psip = I_div_psip
    endif

    allocate(nux0(n_species))
    allocate(nux2(n_species))
    allocate(nux4(n_species))

    ir_global = ir
    do is_global=1, n_species
       do ietype=1, 3
          call gauss_integ(-1.0,1.0,myHSenefunc,integ_order,Nx,eii_val)
          if(ietype == 1) then
             nux0(is_global) = eii_val * 4.0 / (3.0 * sqrt(pi))
          else if(ietype == 2) then
             nux2(is_global) = eii_val * 4.0 / (3.0 * sqrt(pi))
          else if(ietype == 3) then
             nux4(is_global) = eii_val * 4.0 / (3.0 * sqrt(pi))
          endif
       enddo
    enddo

    sum_nm = 0.0
    do is_global=1, n_species
       sum_nm = sum_nm + mass(is_global)*dens(is_global,ir)*nux0(is_global)
    enddo

    do is_global=1, n_species
             
       A1 = -dlnndr(is_global,ir) + 1.5*dlntdr(is_global,ir)
       A2 = -dlntdr(is_global,ir)
       pflux_multi(is_global) = 0.0
       eflux_multi(is_global) = 0.0

       L_a = nux0(is_global) * omega_fac * HS_I_div_psip**2 &
            * rho(ir)**2 * ftrap * dens(is_global,ir) &
            * temp(is_global,ir) * mass(is_global) &
            / (Z(is_global) * Z(is_global) * 1.0)

       do js=1, n_species

          L_b = nux0(js) * omega_fac * HS_I_div_psip**2 &
               * rho(ir)**2 * ftrap * dens(js,ir) &
               * temp(js,ir) * mass(js) / (Z(js) * Z(js) * 1.0)

          if(is_global == js) then
             L11 = -L_a * (sum_nm &
                  - mass(is_global) * dens(is_global,ir) * nux0(is_global)) &
                  / sum_nm
             L12 = L11 * (nux2(is_global) / nux0(is_global))
             L22 = (nux2(is_global) / nux0(is_global)) &
                  * (L12 &
                  + L_a * (nux2(is_global) / nux0(is_global) &
                  - nux4(is_global) / nux2(is_global)))
             L21 = L12
             
          else       
             L11 = L_a * (Z(is_global) * temp(js,ir)) &
                  / (Z(js) * temp(is_global,ir)) &
                  * mass(js) * dens(js,ir) * nux0(js) &
                  / sum_nm
             L12 = (nux2(js) / nux0(js)) * L11
             L21 = (nux2(is_global) / nux0(is_global)) &
                  * Z(js)/(1.0*Z(is_global)) &
                  * (mass(is_global) * dens(is_global,ir) * nux0(is_global) &
                  / sum_nm) * L_b
             L22 = (nux2(is_global) / nux0(is_global)) * L12
             
          end if


          pflux_multi(is_global) = pflux_multi(is_global) &
               + L11 * A1 + L12 * A2
          eflux_multi(is_global) = eflux_multi(is_global) &
               + L21 * A1 + L22 * A2
          
       enddo
    enddo

    deallocate(nux0)
    deallocate(nux2)
    deallocate(nux4)

  end subroutine compute_HS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Ion Poloidal Flow at theta=0 (from kpar)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_vpol_ion(ir, kpar, vpol_ion)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, Bmag2_avg, Bpol
    implicit none
    real, intent (out) :: vpol_ion
    real, intent (in) :: kpar
    integer, intent (in) :: ir
    integer :: it, jt, m_theta
    real :: bigk, sum
    real, dimension(:), allocatable :: vpol_ion_th

    bigk = -kpar *  (dlntdr(is_ion,ir) * I_div_psip * rho(ir) &
         * temp(is_ion,ir) / (z(is_ion)*1.0)) &
         / (Bmag2_avg / dens(is_ion,ir))
    
    allocate(vpol_ion_th(n_theta))
    do it=1,n_theta
       vpol_ion_th(it) = bigk * Bpol(it) / dens(is_ion,ir)
    end do

    vpol_ion = 0.0
    m_theta = (n_theta-1)/2-1
    do jt=0, m_theta
       ! compute the fourier cosine series coefficient c_n:
       ! vpol = sum_n c_n * cos(n*theta)
       sum = 0.0
       do it=1, n_theta
          sum = sum +  vpol_ion_th(it) * cos(jt * theta(it))
       enddo
       if(jt == 0) then
          sum = sum / (1.0*n_theta)
       else
          sum = sum / (0.5*n_theta)
       end if
       ! vpol(theta=0) = sum_n c_n
       vpol_ion = vpol_ion + sum
    enddo

    deallocate(vpol_ion_th)

  end subroutine compute_vpol_ion


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Braun and Helander theory for zf test polarization coefficient
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine THEORY_ZF_do(ir, zf_Pcoeff_BH)
    use neo_globals
    use neo_equilibrium, only : ftrap
    implicit none
    real, intent (inout) :: zf_Pcoeff_BH
    integer, intent (in) :: ir
    real :: fc, eii_val, nufac1, nufac2, nufac3, nufac4
    integer :: is
    integer :: Nx=100
    integer, parameter :: integ_order = 64

    eps = r(ir) / rmaj(ir)
    fc = 1.0 - ftrap

    zf_type = 0
    do is=1, n_species
       if(Z(is) > 0) then
          if(zf_type == 0) then
             zf_ion1 = is
             zf_type = 1
          else if(zf_type == 1) then
             zf_ion2 = is
             zf_type = 2
             exit
          endif
       endif
    enddo

    ir_global = ir
    do ietype=1, 4
       call gauss_integ(-1.0,1.0,myZFenefunc,integ_order,Nx,eii_val)
       if(ietype == 1) then
          nufac1 = eii_val * 4.0 / (3.0 * sqrt(pi))
       else if(ietype == 2) then
          nufac2 = eii_val * 4.0 / (3.0 * sqrt(pi))
       else if(ietype == 3) then
          nufac3 = eii_val * 4.0 / (3.0 * sqrt(pi))
       else if(ietype == 4) then
          nufac4 = eii_val * 4.0 / (3.0 * sqrt(pi))
       endif
    enddo

    zf_Pcoeff_BH = fc * (zf_time) * (nufac1 &
         + fc * nufac2**2 / (nufac3 - fc * nufac4) )

    !print *, nufac1*nu(1,ir), nufac2, &
    !     nufac3/nu(1,ir), nufac4/nu(1,ir)

    !print *, nufac1, fc * nufac2**2 / (nufac3 - fc * nufac4)
    !print *, fc, nufac2, nufac3, nufac4

  end subroutine THEORY_ZF_do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function Defininitions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! uses ir_global, is_global, ietype, eps
  real*8 function myHSenefunc(x)
    use neo_globals
    use neo_energy_grid, only : get_coll_freqs
    use neo_equilibrium, only : ftrap
    implicit none
    integer :: js
    real :: x, xa, xb
    real :: ene, de, val, energy_min
    real :: ft_fac, ft_star, nu_d_tot
    real :: nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p

    ! x grid: -1 -> 1 (equally spaced) for int 0..emax
    ! x = 2 (v/vth)/sqrt(2*emax)-1
    energy_min=0.0
    xa = 2.0 / (1.0 - sqrt(energy_min / energy_max))
    xb = - (1.0 + sqrt(energy_min / energy_max)) &
         / (1.0 - sqrt(energy_min / energy_max))
    ene = energy_max * ( (x-xb) / xa )**2
    de = 2.0 * sqrt(energy_max) / xa
    val = de * exp(-ene) ! weight w/o ene factor
    
    ! nu_d * ene 
    nu_d_tot = 0.0
    do js=1, n_species
       call get_coll_freqs(0, ir_global,is_global,js,ene, &
            nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p)
       nu_d_tot = nu_d_tot + nu_d * nu(is_global,ir_global)
    enddo
    
    ft_star = (3.0*pi/16.0) * eps**2 * vth(is_global,ir_global) &
         * sqrt(2.0) * ene**1.5 &
         / (rmaj(ir_global) * q(ir_global) * nu_d_tot)
    ft_fac  = 1.0 / (1.0 + ftrap / ft_star)
    
    ! interpolated
    if(ietype == 1) then
       myHSenefunc = val * nu_d_tot * ene * ft_fac
    else if(ietype == 2) then
       myHSenefunc = val * nu_d_tot * ene * ene * ft_fac
    else if(ietype == 3) then
       myHSenefunc = val * nu_d_tot * ene * ene * ene * ft_fac

    endif
    
  end function myHSenefunc

  ! uses ietype, ir_global, zf_ion1, zf_ion2, zf_type
  real*8 function myZFenefunc(x)
    use neo_globals
    use neo_energy_grid, only : get_coll_freqs
    implicit none
    real :: x, xa, xb
    real :: ene, de, val, energy_min
    real, external :: derf
    real :: nu_dii, nu_diz
    real :: nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p
    integer :: num, is, js
    
    ! x grid: -1 -> 1 (equally spaced) for int 0..emax
    ! x = 2 (v/vth)/sqrt(2*emax)-1
    energy_min=0.0
    xa = 2.0 / (1.0 - sqrt(energy_min / energy_max))
    xb = - (1.0 + sqrt(energy_min / energy_max)) &
         / (1.0 - sqrt(energy_min / energy_max))
    ene = energy_max * ( (x-xb) / xa )**2
    de = 2.0 * sqrt(energy_max) / xa
    val = de * exp(-ene) ! weight w/o ene factor
    
    ! "nu_d" is nu_d * ene
    is=zf_ion1
    do num=1, zf_type
       if(num == 1) then
          js=zf_ion1
       else
          js=zf_ion2
       endif
       call get_coll_freqs(0, ir_global,is,js,ene, &
            nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p)
       nu_d = nu_d * nu(is,ir_global)
       if(num == 1) then
          nu_dii = nu_d
          if(zf_type == 1) then
             nu_diz = 0.0
          endif
       else
          nu_diz = nu_d
       endif
    enddo
    
    ! func * ene
    if(ietype == 1) then
       nu_d = ene**2 / (nu_dii + nu_diz)
    else if(ietype == 2) then
       nu_d = ene * nu_dii / (nu_dii + nu_diz)
    else if(ietype == 3) then
       nu_d = nu_dii
    else if(ietype == 4) then
       nu_d = nu_dii**2 / (nu_dii + nu_diz)
    end if
       
    ! int is (weight * ene) * func * ene
    myZFenefunc = val * ene * nu_d
    
  end function myZFenefunc

end module neo_theory
