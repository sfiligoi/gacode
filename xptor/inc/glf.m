c@glf.m J. Kinsey, General Atomics
c 30-may-06 added cnormi_tg, cnorme_tg, cnormd_tg, chii_ib_tg, chie_ib_tg
c           and chii_eb_tg, chie_eb_tg, diff_ib_tg, diff_eb_tg
c 23-may-06 added particle_flux_gf(2,2),energy_flux_gf(2,2)
c 09-may-06 added cnorme_gf, cnormi_gf norms for sep. branches
c           also added diff_ib_gf, diff_eb_gf, chii_ib_gf, chii_eb_gf
c           chie_ib_gf, chie_eb_gf
c 03-may-06 added iroot_ion, iroot_elec
c 17-dec-04 added igauss_gf  for glf2d_663j.F
c 14-sep-04 added ibranch_gf for glf2d_max.f
c 04-may-04 added theta0_gf to be used w/ glf2d_S663.F
c 21-oct-03 added n_tilda_gf
c 18-jun-03 added cnorm_p_gf for part. diffusivity
c 13-jun-02 added gamma_r_gf for radial mode damping rate
c 16-apr-01 added iglf
c 05-mar-01 changed 20 to nmode
c 23-aug-00 aligned common block, added xky_gf
c 14-june-00 added ngrow_k_gf
c 13-june-00 added ipert_gf
c 20-dec-99 added igfac
c 03-aug-99 added jeigen
c---------------------------------------------------------------------
      integer nmode
      parameter (nmode=20)
      integer iflagin_gf(30), ngrow_k_gf(0:nmode)
      integer nroot_gf, nbasis_gf, nx_gf, jeigen, igfac
     & , lprint_gf, ikymax_gf, eigen_gf
     & , i_err, i_proc, n_proc, first_order_gf, ipert_gf
     & , iglf, ipade_gf, ibranch_gf, igauss_gf, iflux_gf, igeo_gf
     & , iroot_ion(2), iroot_elec(2)

      real*8 yparam_k_gf(nmode,nmode)
     & , gamma_k_gf(1:4,nmode),freq_k_gf(1:4,nmode)
     & , phi_norm_k_gf(1:4,nmode)
     & , xparam_gf(30)
     & , xkyf_k_gf(nmode),diff_k_gf(nmode)
     & , diff_im_k_gf(nmode),chii_k_gf(nmode)
     & , chie_k_gf(nmode),exch_k_gf(nmode)
     & , eta_par_k_gf(nmode),eta_per_k_gf(nmode),eta_phi_k_gf(nmode)
     & , chie_e_k_gf(nmode), n_tilda2_k_gf(nmode), yparam_gf(nmode)
     & , gamma_gf(1:4),freq_gf(1:4),phi_norm_gf(1:4),xky_gf(1:4)
     & , particle_flux_gf(2,2),energy_flux_gf(2,2)
     & , xky0_gf,rms_theta_gf,theta0_gf,rlti_gf
      real*8 rlte_gf,rlne_gf,rlni_gf,rlnimp_gf,dil_gf,apwt_gf
     & , aiwt_gf,taui_gf,rmin_gf,rmaj_gf,q_gf,xnu_gf,betae_gf
     & , shat_gf,alpha_gf,elong_gf,xwell_gf,park_gf,phii_gf,ghat_gf
     & , gchat_gf,adamp_gf,alpha_star_gf,gamma_star_gf
     & , alpha_e_gf,gamma_e_gf,alpha_mode_gf,gamma_mode_gf
     & , alpha_p_gf,gamma_p_gf,gamma_r_gf,xkdamp_gf,xkyf_gf
     & , diff_gf,diff_im_gf,chii_gf,chie_gf
     & , diff_ib_gf, diff_eb_gf, chii_ib_gf, chii_eb_gf
     & , chie_ib_gf, chie_eb_gf, exch_gf,eta_par_gf
     & , eta_per_gf,eta_phi_gf,chie_e_gf,n_tilda_gf
     & , cnorm_p_gf,cnorm_gf,cnorme_gf,cnormi_gf
     & , xkymin_gf,xkymax_gf,amassgas_gf,amassimp_gf,zimp_gf
     & , cnormi_tg, cnorme_tg, cnormd_tg
     & , chii_ib_tg, chie_ib_tg
     & , chii_eb_tg, chie_eb_tg, diff_ib_tg, diff_eb_tg

      complex*16 zevec_k_gf(nmode,12,12), zomega_k_gf(nmode,12)

      common /glf/ zevec_k_gf, zomega_k_gf
     & , yparam_k_gf,gamma_k_gf,freq_k_gf,phi_norm_k_gf
     & , xparam_gf,xkyf_k_gf,diff_k_gf
     & , diff_im_k_gf,chii_k_gf,chie_k_gf,exch_k_gf
     & , eta_par_k_gf,eta_per_k_gf,eta_phi_k_gf
     & , chie_e_k_gf,n_tilda2_k_gf,yparam_gf
     & , gamma_gf,freq_gf,phi_norm_gf,xky_gf
     & , particle_flux_gf,energy_flux_gf
     & , xky0_gf,rms_theta_gf,theta0_gf
     & , rlti_gf,rlte_gf,rlne_gf,rlni_gf,rlnimp_gf,dil_gf
     & , apwt_gf,aiwt_gf,taui_gf,rmin_gf,rmaj_gf,q_gf,xnu_gf
     & , betae_gf,shat_gf,alpha_gf,elong_gf,xwell_gf,park_gf
     & , phii_gf,ghat_gf,gchat_gf,adamp_gf,alpha_star_gf
     & , gamma_star_gf,alpha_e_gf,gamma_e_gf,alpha_mode_gf
     & , gamma_mode_gf,alpha_p_gf,gamma_p_gf,gamma_r_gf,xkdamp_gf
     & , xkyf_gf,diff_gf,diff_im_gf,chii_gf,chie_gf,exch_gf
     & , eta_par_gf,eta_per_gf,eta_phi_gf,chie_e_gf
     & , diff_ib_gf, diff_eb_gf, chii_ib_gf, chii_eb_gf
     & , chie_ib_gf, chie_eb_gf
     & , n_tilda_gf,cnorm_gf,cnorm_p_gf,cnorme_gf,cnormi_gf
     & , cnormi_tg, cnorme_tg, cnormd_tg, chii_ib_tg, chie_ib_tg
     & , chii_eb_tg, chie_eb_tg, diff_ib_tg, diff_eb_tg
     & , xkymin_gf,xkymax_gf
     & , amassgas_gf,amassimp_gf,zimp_gf
     & , iflagin_gf,ngrow_k_gf,iroot_ion,iroot_elec
     & , nroot_gf,nbasis_gf,nx_gf,lprint_gf,ikymax_gf
     & , i_err,i_proc,n_proc,eigen_gf,igfac
     & , first_order_gf,ipert_gf, iglf, ipade_gf 
     & , ibranch_gf, igauss_gf, iflux_gf, igeo_gf

