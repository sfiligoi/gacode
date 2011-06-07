      subroutine setup
c
c setup for glf23_v4_1_10
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      integer j
c
c setup for glf23_v4_1_10
c
      nroot_gf=8
      iflagin_gf(1)=0
      iflagin_gf(2)=1
      iflagin_gf(3)=1
      iflagin_gf(5)=3
      xparam_gf(3)=0.7D0
      xparam_gf(6)=0.D0
      xparam_gf(13)=.2D0
      xparam_gf(14)=1.D0
      xparam_gf(15)=-.1D0
      xparam_gf(16)=0.D0
      xparam_gf(17)=.1D0
      xparam_gf(18)=0.D0
      xky0_gf= .2D0
      rms_theta_gf=3.14159265/3.D0
      park_gf  =0.7D0
      ghat_gf  =1.D0
      gchat_gf =1.D0
      adamp_gf=.75D0
      alpha_star_gf=0.D0
c       alpha_e_gf=0.D0
c       xparam_gf(10)=1.D0
      gamma_e_gf=0.D0
      xkdamp_gf=0.D0
      alpha_p_gf=0.5D0
      xparam_gf(7)=1
      xparam_gf(9)=1.D0
      cbetae=1.D-6
      cnorm_gf=100.D0
      cnorm_p_gf=100.D0
      ikymax_gf=10
      xkymin_gf=.02D0
      xkymax_gf=.5D0
      cmodel=1.D0
      xalpha=1.D0
      ialphastab=0
      ineo=-2
      adiffphi_dv=0.D0
      xchi=0.1D0
      chiaddexp=1.D0
      iexch_m=1
      iohm=0
      irad=0
      ineutp=0
      idt=0
      jin_m=1
      i_dengrad=2
      iexp_imp=1
      igeo_m=3
      igfac=0
c       irot1=1
c       irot2=1
      adamp_gf=.5D0
      xparam_gf(23)=1
      lprint_gf=0
      iparam_pt(6)=1
      mscale=1.D0
      nfscale=1.D0
      btscale=1.D0
      prfscale=1.D0
      prfescale=1.D0
      prfiscale=1.D0
      pbescale=1.D0
      pbiscale=1.D0
      ipfst=0
      fuscale=1.D0
      itest_ntcc=0
      ialpha=0
      ifusmodel=0
      pfusion_max=2000.D0
      ibound=0
      cped=1.D0
      ipade_gf=0
      ibranch_gf=0
      theta0_gf=0.0
      nbasis_gf=4
      xoh_exp=0.D0
      xfus_exp=0.D0
c
      return
      end
