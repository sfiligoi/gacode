c@input.m XPTOR, General Atomics
c 14-feb-11 removed igks
c 17-nov-05 (jek) removed _sc variables
c 01-sep-04 added idatzero for new 0D ufiles
c 12-sep-02 added cped for multiplier on pedestal scaling
c 19-aug-02 added ibound, t_ped(8) to use pedestal scaling for BC
c 10-jun-02 added ipfst switch for pfast in alpha_m
c 12-may-02 added iptotr to read in ptotr from iterdb files
c 16-nov-01 added iglf2d to call glf2d from driver
c 16-oct-01 added wallneutp to scale neutral power flow
c 19-sep-01 added igks to call GKS code
c 29-jun-01 added ifilter
c 22-feb-01 added zgas_exp
c 18-jan-01 added iscan_exb, iscan_alpha, iscan_cmod, cmodel_scan
c 12-jan-01 added itorq
c 01-nov-00 added ave_field array over fields, added xsdot
c 08-sep-00 added iptot,ipramp to print Ip vs time for current ramps
c 29-aug-00 added ave_ve for spatial averaging of ve
c 10-mar-00 added lprint_pflow
c 22-feb-00 added didledge_vphi
c 11-feb-00 added ilog used in ptor_dv.f
c 05-nov-99 added bt_flag switch to allow user to change from
c           bt_exp to bteff_exp(j) (new variable in tport.m)
c 29-sept-99 added ncl_flag to allow reading in NCLASS variables
c---------------------------------------------------------------------
      real*8 t_ped(8), ave_field(5)
c
      real*8 d_m,dpiv_m,dxe_m,dxi_m,dxn_m,limin,limout,didledge_te
      real*8 didledge_ti,didledge_n,didledge_vphi,rsigma
      real*8 zeff_sc,zeff_e
      real*8 zpmte_hold,zpmti_hold,zpmne_hold,zpmni_hold,betae_hold
      real*8 xnu_hold,rhos_hold,q_hold,shat_hold,eps_hold,rmaj_hold
      real*8 tiote_hold,zeff_hold,elong_hold
      real*8 wallneut,wallneutp,corot,xwdot,xsdot,arho_exp,elonga_exp
      real*8 deltaa_exp,rmajor_exp,zimp_exp,amassimp_exp
      real*8 bt_exp,sign_bt_exp,bp0_exp,amassgas_exp,zgas_exp,te0_exp
      real*8 ti0_exp,alfte_exp,ftea_exp,alfti_exp,ftia_exp,ne0_exp
      real*8 alfne_exp,xpnq_exp,fnea_exp,frad_exp
      real*8 alfrad_exp,powe0_exp,alfpowe_exp,powi0_exp
      real*8 alfpowi_exp,flow0_exp,alfflow_exp
      real*8 arfe_exp,brfe_exp,wrfe_exp,arfi_exp
      real*8 brfi_exp,wrfi_exp,q0_exp,qa_exp
      real*8 alfj_exp,qmin_exp,rho_qm_exp,betarad,alfarad,epsrad
      real*8 temult, timult, xmult, adiff_dv, ave_dv, ave_ve
      real*8 adiffphi_dv, curtot, cped
c
      integer tauscale(4)
      integer istk,nstk,istkmax
      integer ishoot,imethod,jmaxm,jin_m,jout_m,jzmn_out
      integer jzmn_in,jd_m,nmaxwn,iscan_m,iscanm_m
      integer idata,nprimd,nneud,niond
      integer nplasbdryd,ntime_d,ntimemax,islice_d
      integer ncl_flag, bt_flag, ilog, lprint_pflow
      integer itorque, iptotr, ipfst, ipramp
      integer ifilter, ifixeta, iglf2d, ibound
      integer idatzero, igyro
c
      common /input/ t_ped, ave_field
     & , d_m, dpiv_m, dxe_m, dxi_m, dxn_m, limin, limout
     & , didledge_te, didledge_ti, didledge_n, didledge_vphi
     & , rsigma, zeff_sc,zeff_e, zpmte_hold, zpmti_hold
     & , zpmne_hold, zpmni_hold, betae_hold, xnu_hold
     & , rhos_hold, q_hold, shat_hold, eps_hold, rmaj_hold
     & , tiote_hold, zeff_hold, elong_hold
     & , wallneut, wallneutp, corot, sign_bt_exp
     & , xwdot, xsdot, arho_exp, elonga_exp, deltaa_exp
     & , rmajor_exp, zimp_exp,amassimp_exp, bt_exp, bp0_exp
     & , amassgas_exp, zgas_exp, te0_exp, ti0_exp, alfte_exp
     & , ftea_exp, alfti_exp, ftia_exp, ne0_exp, alfne_exp
     & , xpnq_exp, fnea_exp, frad_exp, alfrad_exp, powe0_exp
     & , alfpowe_exp, powi0_exp, alfpowi_exp, flow0_exp
     & , alfflow_exp, arfe_exp, brfe_exp, wrfe_exp, arfi_exp
     & , brfi_exp, wrfi_exp, q0_exp, qa_exp, alfj_exp
     & , qmin_exp, rho_qm_exp, betarad, alfarad, epsrad
     & , temult, timult, xmult, adiff_dv, ave_dv, ave_ve
     & , adiffphi_dv, curtot, cped, tauscale
     & , istk, nstk, istkmax
     & , ishoot, imethod, jmaxm
     & , jin_m, jout_m, jzmn_out, jzmn_in, jd_m
     & , nmaxwn, iscan_m, iscanm_m
     & , idata, idatzero, nprimd, nneud, niond
     & , nplasbdryd, ntime_d, ntimemax
     & , islice_d, ncl_flag, bt_flag, ilog, lprint_pflow
     & , itorque, iptotr, ipfst, ipramp, ifilter
     & , ifixeta, iglf2d, ibound, igyro

