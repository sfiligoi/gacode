      SUBROUTINE tglf_startup
!
      USE tglf_pkg  
      USE tglf_tg
!  
      IMPLICIT NONE
!
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
!
      INTEGER :: k
      real*8 qm,rmajm,rminm,rhom
      real*8 dr_loc,b_unit
!
! read the namelist 
!
      OPEN (unit=3,file='tglfin',status='old')
      READ(3,nml=tglfin)
      CLOSE(3)
!
      if(i_proc.eq.0)write(*,*)"sign_bt=",sign_Bt_exp
!
      CALL put_rare_switches(theta_trapped_tg,park_tg,ghat_tg,
     > gchat_tg, wd_zero_tg,Linsker_factor_tg,
     > gradB_factor_tg,filter_tg,damp_psi_tg,damp_sig_tg)
!
      CALL put_switches(iflux_tg,use_bper_tg,use_bpar_tg,
     > use_mhd_rule_tg,use_bisection_tg,ibranch_tg,
     > nmodes_tg,nbasis_max_tg,nbasis_min_tg,nxgrid_tg,nky_tg)
!
      CALL put_species(ns_tg,zs_tg,mass_tg)
!
      CALL put_model_parameters(adiabatic_elec_tg,alpha_e_tg,alpha_p_tg,
     > alpha_n_tg,alpha_t_tg,alpha_kx_e_tg,alpha_kx_p_tg,    
     > alpha_kx_n_tg,alpha_kx_t_tg,alpha_quench_tg,xnu_factor_tg,
     > debye_factor_tg,etg_factor_tg,sat_rule_tg,kygrid_model_tg,
     > xnu_model_tg,vpar_model_tg,vpar_shear_model_tg)
!
      CALL put_kys(ky_tg)
!
      CALL put_gaussian_width(width_max_tg,width_min_tg,nwidth_tg    
     > ,find_width_tg)
!
      CALL put_gradients(rlns_tg,rlts_tg,vpar_shear_tg,vexb_shear_tg)
!
      CALL put_averages(taus_tg,as_tg,vpar_tg,betae_tg,xnue_tg,
     > zeff_tg,debye_tg)
!
      sign_Bt_tg = sign_bt_exp
      sign_It_tg=1.0
      CALL put_signs(sign_Bt_tg,sign_It_tg)
!
      a_unit_exp = rmin_exp(mxgrid)
      if(igeo_tg.eq.0)a_unit_exp=arho_exp
!
      do k=0,mxgrid-1
      jm = k
      qm =(q_exp(jm+1)+q_exp(jm))/2.D0
      rmajm=(rmaj_exp(jm+1)+rmaj_exp(jm))/2.D0
      rhom=arho_exp*(rho(jm+1)+rho(jm))/2.D0
c recompute drhodr just to make sure
      drhodr(jm) = arho_exp*(rho(jm+1)-rho(jm))/
     >             (rmin_exp(jm+1)-rmin_exp(jm))
      rminm=(rmin_exp(jm+1)+rmin_exp(jm))/2.D0
c   local magnetic field unit 
      if(igeo_tg.eq.0)then
        b_unit = bt_exp
c        if(bt_flag.gt.0)b_unit = bteff_exp(jm)
      else
        b_unit = bt_exp*(rhom/rminm)*drhodr(jm)
      endif
       if(igeo_tg.eq.0)then
c s-alpha geometry
        rmin_tg = rminm/a_unit_exp
        rmaj_tg = rmajm/a_unit_exp
        q_tg = qm
        shat_tg = shat_exp(jm)
        alpha_tg=xalpha*alpha_exp(jm)
        if(ialphastab.gt.0) then
          alpha_tg=xalpha*alpha_m(jm)
        endif
        CALL put_s_alpha_geometry(rmin_tg,rmaj_tg,q_tg,shat_tg,alpha_tg, 
     >    xwell_tg,theta0_tg,b_model_tg,ft_model_tg)
c
       elseif(igeo_tg.eq.1)then  ! miller geometry
c
         dr_loc = (rmin_exp(jm+1)-rmin_exp(jm))/a_unit_exp
         rmin_tg = rminm/a_unit_exp
         rmaj_tg = rmajm/a_unit_exp
         zmaj_tg=0.0
         q_tg = qm
         q_prime_tg = (q_tg/rmin_tg)*(q_exp(jm+1)-q_exp(jm))/dr_loc
         p_prime_tg = (q_tg/rmin_tg)*(1.6022D-4/b_unit**2)*
     >                (ptot_exp(jm+1)-ptot_exp(jm))/dr_loc
         if(ialphastab.eq.1)then
           p_prime_tg = (q_tg/rmin_tg)*(1.6022D-4/b_unit**2)*
     > (nem*gradtem+tem*gradnem+(nim+nzm)*gradtim+tim*(gradnim+gradnzm))
         endif
         drmindx_tg=1.0
         drmajdx_tg = (rmaj_exp(jm+1)-rmaj_exp(jm))/(dr_loc*a_unit_exp)
         dzmajdx_tg=0.0
         kappa_tg = 0.5*(elong_exp(jm+1)+elong_exp(jm))
         s_kappa_tg = (rmin_tg/kappa_tg)*
     >                (elong_exp(jm+1)-elong_exp(jm))/dr_loc
         delta_tg = 0.5*(delta_exp(jm+1)+delta_exp(jm))
         s_delta_tg = rmin_tg*(delta_exp(jm+1)-delta_exp(jm))/dr_loc
         zeta_tg = 0.0
         s_zeta_tg=0.0
        CALL put_Miller_geometry(rmin_tg,rmaj_tg,zmaj_tg,drmindx_tg,
     >  drmajdx_tg,dzmajdx_tg,kappa_tg,s_kappa_tg,delta_tg,s_delta_tg,
     >  zeta_tg,s_zeta_tg,q_tg,q_prime_tg,p_prime_tg)
c
c
       else
        write(*,*)"igeo_tg invalid",igeo_tg
        stop
       endif
!     this is first call to tglf so it will precompute hermite basis functions
!     nb_max_tg should be the maximum number of basis functions that will 
!     be used for all future calls
!
        CALL tglf_setup_geometry
!
        xr2_exp(k) = get_R2_ave()*a_unit_exp**2
        xb2_exp(k) = get_B2_ave()*B_unit**2
        f_exp(k) = rmajor_exp*bt_exp/(get_RBt_ave()*B_unit*a_unit_exp)
        c_par(k) = xb2_exp(k)/Bt_exp**2
        c_tor(k) = xr2_exp(k)/rmajor_exp**2
        c_per(k) = 1.0/f_exp(k)
        a_pol(k) = get_a_pol()*B_unit/bt_exp
        a_tor(k) = get_a_tor()*rmaj_tg*a_unit_exp/rmajor_exp
        Bp0(k) = get_Bp0()*B_unit
        interchange_DR_m(k) = get_DR()
      enddo
      k=mxgrid
      xr2_exp(k) = xr2_exp(k-1)
      xb2_exp(k) = xb2_exp(k-1)
      f_exp(k) = f_exp(k-1)
      c_par(k) = c_par(k-1)
      c_tor(k) = c_tor(k-1)
      c_per(k) = c_per(k-1)
      a_pol(k) = a_pol(k-1)
      a_tor(k) = a_tor(k-1)
      Bp0(k) = Bp0(k-1)
      interchange_DR_m(k) = 0.0
      k=0
      xr2_exp(k) = xr2_exp(k+1)
      xb2_exp(k) = xb2_exp(k+1)
      f_exp(k) = f_exp(k+1)
      c_par(k) = c_par(k+1)
      c_tor(k) = c_tor(k+1)
      c_per(k) = c_per(k+1)
      a_pol(k) = a_pol(k+1)
      a_tor(k) = a_tor(k+1)
      Bp0(k) = 0.5*Bp0(k+1)
!
      return
      end
!

