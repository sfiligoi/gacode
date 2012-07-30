c@ptor.m 05-dec-02 J. Kinsey, General Atomics
c 05-dec-02 added ngin for innermost grid pt
c 03-apr-02 added time1_pow-time3_pow for time-dep power
c 27-nov-01 added taupin for default ptcle confinement time, time_neut
c 21-sep-01 added Prads, radref for synchrotron radiation
c 25-apr-01 added ineutp for neutrals, charge-exchange, recombination
c 15-mar-01 added ifusmodel for choice of fusion power calculation
c           and idt to use n_D and n_T from ufiles, and fuscale
c           to scale up fusion power to match experiment
c 06-mar-01 added lprint_time_cdf to store and write exp time data
c           and time3_pt, time_ipulse_pt, yipulse_pt for ion pulse
c 18-jan-01 added sigma, vexb_bc, Psour_wall
c 01-nov-00 added dznpro
c 28-sep-00 added istep_glf to log step no. for diag. in model.f
c           and lprint_glfstep, lprint_glfgrid for selecting step, grid pt.
c 11-may-00 added flux and flux_s, changed mxgrd from 50 to 301
c 04-apr-00 added time_upd and upsrc(8) to update sources in ptor
c 22-mar-00 added and ichiep, ktime_chiep to prescribe Gaussian chie
c           at a given time-step set by ktime_chiep
c 11-feb-00 added mxflds, endtime_pt, restart_pt, sramp_pt, smult
c           s0_pt, s1_pt
c---------------------------------------------------------------------
      integer mxgrd, mxflds
      parameter (mxgrd=301, mxflds=5)
c..1505 elements
      real*8 flux(mxflds,mxgrd),flux_s(mxflds,mxgrd)
      real*8 r(mxgrd,2), vprime(mxgrd,2), dr(mxgrd,2)
c..301 elements
      real*8 qpro(mxgrd), denpro(mxgrd), dinpro(mxgrd)
      real*8 dznpro(mxgrd)
      real*8 Tipro(mxgrd), Tepro(mxgrd), Zeffpro(mxgrd)
      real*8 vparpro(mxgrd), vperpro(mxgrd), vphipro(mxgrd)
      real*8 Peaux(mxgrd), Piaux(mxgrd), Psour(mxgrd)
      real*8 Psour_wall(mxgrd)
      real*8 Mphi(mxgrd), Mpar(mxgrd), Mper(mxgrd), psume(mxgrd)
      real*8 psumi(mxgrd), Pe_alpha(mxgrd), Pi_alpha(mxgrd)
      real*8 Pradb(mxgrd), Prads(mxgrd), Pohpro(mxgrd), chii(mxgrd)
      real*8 chie(mxgrd), chiei(mxgrd), chiie(mxgrd)
      real*8 convTi(mxgrd), convTe(mxgrd), Dpart(mxgrd)
      real*8 Vpart(mxgrd), convV(mxgrd)
      real*8 eta_par(mxgrd), eta_per(mxgrd), eta_phi(mxgrd)
      real*8 conv_phi(mxgrd), nu_ie(mxgrd), nu_ei(mxgrd)
      real*8 chie_ptor(mxgrd), chii_ptor(mxgrd)
      real*8 v2_bar(mxgrd)
c..20 elements
      real*8 xparam_pt(1:30)
c..5 elements
      real*8 sramp_pt(mxflds), smult(mxflds)
      real*8 s0_pt(mxflds), s1_pt(mxflds)
c..1 element
      real*8 Rmaj,amin,kappa,Btor,qa,q0,Paux_e,Paux_i,pfusion_max
      real*8 chiip_mult,dt_v,delt_v,d_art
      real*8 time_pulse_pt,time_ipulse_pt,plse_pt,plsi_pt
      real*8 time_del_pt,duty_pt,ancnt_pt,xon_pt
      real*8 tau_sd_pt,p_pulse_pt,p_pulse2_pt,t_pulse_pt
      real*8 epulse_pt,ipulse_pt,time1_pt,ypulse_pt,yipulse_pt
      real*8 time1_pow, time2_pow, time3_pow
      real*8 pbescale2, pbiscale2, pbescale3, pbiscale3
      real*8 prfscale2, prfscale3
      real*8 temp_te,temp_ti,time2_pt,time3_pt,pi,kJpereV
      real*8 time,time_upd,Pfusion,Pradb_tot,Prads_tot,Pohmic_tot
      real*8 diag_off,diag_ae,diag_ai,diag_an,diag_na,diag_ne
      real*8 diag_ni,diag_ea,diag_ei,diag_en,diag_ia,diag_ie,diag_in
      real*8 chiep_a, chiep_b, chiep_c
      real*8 endtime_pt, restart_pt, sigma, vexb_bc
      real*8 fuscale, radref, taupin, time_neut, norm_s, norm_dt
c
      integer iparam_pt(1:20),upsrc(8),itport_pt(5)
      integer mxgrid,mxfields,ngin,ngrid,nsteps_v,dvflag,time_series
      integer iohm,ialpha,iexch,irad,ineutp,ifusmodel,idt
      integer k1_pt,k2_pt,k3_pt,k4_pt,ifix,izero_m,idvloop_max
      integer ichiep, ktime_chiep
      integer lprint_glf2, istep_glf, lprint_glfstep, lprint_glfgrid,
     &        lprint_time_cdf
      integer doppler_shear_model, j_write_state
c
      common /ptorcmn/ flux, flux_s
     & , r, vprime, dr
     & , qpro, denpro, dinpro, dznpro, Tipro
     & , Tepro, Zeffpro, vparpro
     & , vperpro, vphipro, Peaux
     & , Piaux, Psour, Psour_wall, Mpar
     & , Mphi, Mper, psume, psumi
     & , Pe_alpha, Pi_alpha, Pradb, Prads
     & , Pohpro, chii, chie
     & , chiei, chiie, convTi
     & , convTe, Dpart, Vpart
     & , convV, eta_par, eta_per
     & , eta_phi, conv_phi, nu_ie
     & , nu_ei, chie_ptor, chii_ptor, xparam_pt
     & , sramp_pt, smult, s0_pt, s1_pt
     & , Rmaj, amin, kappa, Btor, qa, q0, Paux_e, Paux_i
     & , pfusion_max, chiip_mult, dt_v, delt_v, d_art
     & , time_pulse_pt, time_ipulse_pt
     & , plse_pt, plsi_pt, time_del_pt
     & , duty_pt, ancnt_pt, xon_pt, tau_sd_pt, p_pulse_pt
     & , p_pulse2_pt, t_pulse_pt, epulse_pt, ipulse_pt
     & , time1_pt, ypulse_pt, yipulse_pt, temp_te, temp_ti
     & , time2_pt, time3_pt, time1_pow, time2_pow
     & , time3_pow, pbescale2, pbiscale2 
     & , pbescale3, pbiscale3, prfscale2, prfscale3
     & , pi, kJpereV, time, time_upd
     & , Pfusion, Pradb_tot, Prads_tot, Pohmic_tot
     & , diag_off, diag_ae, diag_ai, diag_an, diag_na
     & , diag_ne, diag_ni, diag_ea, diag_ei, diag_en
     & , diag_ia, diag_ie, diag_in, chiep_a, chiep_b, chiep_c
     & , endtime_pt, restart_pt, sigma, vexb_bc
     & , fuscale, radref, taupin, time_neut, norm_s, norm_dt
     & , v2_bar, iparam_pt, upsrc, itport_pt
     & , mxgrid, mxfields, ngin, ngrid, nsteps_v, dvflag, time_series
     & , iohm, ialpha, iexch, irad, ineutp, ifusmodel, idt
     & , k1_pt, k2_pt, k3_pt, k4_pt
     & , ifix, izero_m, idvloop_max
     & , ichiep, ktime_chiep
     & , istep_glf, lprint_glfstep, lprint_glfgrid, lprint_glf2
     & , lprint_time_cdf, doppler_shear_model, j_write_state


