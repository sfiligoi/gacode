c@data_t.m 30-may-02 G. Staebler, General Atomics
c common block for _t variables for GKS stand alone code
c parameters are defined in data_exp.m
c..725,801 elements
      real*8 te_t(0:jmaxmt,0:kmaxmt), ti_t(0:jmaxmt,0:kmaxmt)
      real*8 te_exp_t(0:jmaxmt,0:kmaxmt), ti_exp_t(0:jmaxmt,0:kmaxmt)
      real*8 vphi_exp_t(0:jmaxmt,0:kmaxmt)
      real*8 ni_t(0:jmaxmt,0:kmaxmt), ne_t(0:jmaxmt,0:kmaxmt)
      real*8 nz_t(0:jmaxmt,0:kmaxmt)
      real*8 vper_t(0:jmaxmt,0:kmaxmt), vphi_t(0:jmaxmt,0:kmaxmt)
      real*8 xte_t(0:jmaxmt,0:kmaxmt), xti_t(0:jmaxmt,0:kmaxmt)
      real*8 vexb_t(0:jmaxmt,0:kmaxmt), vetor_t(0:jmaxmt,0:kmaxmt)
      real*8 vepol_t(0:jmaxmt,0:kmaxmt), vstar_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_t(0:jmaxmt,0:kmaxmt),anratem_t(0:jmaxmt,0:kmaxmt)
      real*8 etaphi_t(0:jmaxmt,0:kmaxmt), gammap_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_vphi_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_vpol_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_vstar_t(0:jmaxmt,0:kmaxmt)
      real*8 chie_t(0:jmaxmt,0:kmaxmt), chii_t(0:jmaxmt,0:kmaxmt)
      real*8 pech_t(0:jmaxmt,0:kmaxmt)
      real*8 powe_t(0:jmaxmt,0:kmaxmt), powi_t(0:jmaxmt,0:kmaxmt)
c..2501 elements
      real*8 time_t(0:kmaxmt), pechtot_t(0:kmaxmt)
      real*8 powe_exp_t(0:kmaxmt), powi_exp_t(0:kmaxmt)
      real*8 powewdot_exp_t(0:kmaxmt), powiwdot_exp_t(0:kmaxmt)
      integer ntime_t
      common /data_t/ te_t, ti_t, te_exp_t, ti_exp_t, vphi_exp_t
     & , ne_t, ni_t, nz_t, powe_t, powi_t, vphi_t, vper_t
     & , xte_t, xti_t, vexb_t, vetor_t, vepol_t, vstar_t, egamma_t 
     & , egamma_vphi_t, egamma_vpol_t, egamma_vstar_t, anratem_t
     & , chie_t, chii_t, etaphi_t, gammap_t, pech_t, time_t
     & , pechtot_t, powe_exp_t, powi_exp_t 
     & , powewdot_exp_t, powiwdot_exp_t
     & , ntime_t
