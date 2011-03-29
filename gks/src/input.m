c@input.m 03-Jun-02 G. Staebler, General Atomics
c 25-feb-11 added sign_bt_exp, bp0_exp, iexb, ialpha, lprint_pflow
c           wallneut,wallneutp,xwdot,xsdot,curtot,idatzero,ipfst,
c           itest_ntcc, q0_exp, qa_exp, qmin_exp, rho_qm_exp
c           removed endtime_pt, time_series, mxgrid
c---------------------------------------------------------------------
c
      real*8 wallneut,wallneutp,corot,xwdot,xsdot
      real*8 arho_exp,elonga_exp
      real*8 deltaa_exp,rmajor_exp,zimp_exp,amassimp_exp
      real*8 bt_exp,sign_bt_exp, bp0_exp, amassgas_exp,zgas_exp
      real*8 q0_exp, qa_exp, qmin_exp, rho_qm_exp
      real*8 temult, timult, xmult
      real*8 curtot
      real*8 xvarmin,xvarmax,zvarmin,zvarmax
      real*8 echconv
c
      integer istk,nstk,istkmax
      integer jmaxm,jin_m,jout_m
      integer idata,nprimd,nneud,niond
      integer nplasbdryd,ntime_d,ntimemax,islice_d
      integer bt_flag,ncl_flag,igks, lprint_cdf
      integer itorque, iptotr, ipfst, lprint_pflow
      integer idatzero,itest_ntcc
      integer iscan_exb, iscan_alpha
      integer ivar,jvar,irunmax,jrunmax,igraph
      integer ipptot
      integer irot1,irot2,irotstab
      integer ialphastab,igeo_m,i_dengrad
      integer ismooth_data,gks_defaults,tglf_defaults
      integer igks_model,nbasis_min,nbasis_max
      integer ibranch, iexb, ialpha
c
      logical save_tglf,overwrite_tglf
c
      common /input/ wallneut, wallneutp, corot, sign_bt_exp 
     & , xwdot, xsdot, xvarmin, xvarmax, zvarmin, zvarmax, echconv
     & , arho_exp, elonga_exp, deltaa_exp
     & , rmajor_exp, zimp_exp, amassimp_exp, bt_exp, bp0_exp
     & , amassgas_exp, zgas_exp, q0_exp, qa_exp
     & , qmin_exp, rho_qm_exp
     & , temult, timult, xmult, curtot
     & , istk, nstk, istkmax, jmaxm
     & , jin_m, jout_m, ibranch, iexb, ialpha
     & , idata, nprimd, nneud, niond
     & , nplasbdryd, ntime_d, ntimemax,islice_d
     & , bt_flag, ncl_flag, igks, lprint_cdf
     & , itorque, iptotr, ipfst, lprint_pflow, idatzero, itest_ntcc
     & , iscan_exb, iscan_alpha, ivar, jvar, igraph
     & , irunmax, jrunmax, ipptot
     & , irot1,irot2,irotstab,ialphastab,igeo_m
     & , i_dengrad,ismooth_data, igks_model
     & , nbasis_min,nbasis_max
     & , gks_defaults,tglf_defaults
     & , save_tglf,overwrite_tglf

