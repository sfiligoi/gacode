       module gks_var
********* start parameters
!variable descripter file replaces common.inc
       integer, parameter :: negmax=10
       integer, parameter :: ntm=80
!ntm = ntheta/2 + (nperiod-1)*ntheta
!ntm=48 for nperiod<=2, ntheta=32
!ntm=80 for nperiod<=3, ntheta=32
!ntm=192 for nperiod=2, ntheta=128
!ntm=384 for nperiod=2, ntheta=256
!ntm=224 for nperiod=4, ntheta=64
        integer, parameter :: ntml=ntm
!nlm = 2*ngauss + ntheta/2 + 1
!ngauss <=5
        integer, parameter :: nlm=27
! for ntheta=32 nlm=27
! for ntheta=64 nlm=43
! for ntheta=128 nlm=75
! for ntheta=256 nml=139
! note nspecmax must be at least 3
        integer, parameter :: nspecmax=5
        integer, parameter :: nlmt2=2*nlm
        integer, parameter :: ntmp1=ntm+1
        integer, parameter :: ntmp1t2=2*ntmp1
        integer, parameter :: ntmp2p1t2=2*(2*ntm+1)
        integer, parameter :: ntmp1t3=3*ntmp1
        integer, parameter :: ntmp2p1t3=3*(2*ntm+1)
        integer, parameter :: nkgrid=10
        integer, parameter :: ipgrid=1
        integer, parameter :: ntime=6000
********* end parameters
*
********* start complex variables:
       complex gkss(147,147)                                 !complex
       complex  gfss(147)                                         !complex
       complex  uss(147)                                          !complex
       complex phisav(-ntml:ntm)                            !complex
       complex   aparsav(-ntml:ntm)                       !complex
       complex   bparsav(-ntml:ntm)                       !complex
       complex   gsav(-ntml:ntm,nlm,negmax,2,nspecmax)           !complex
       complex phi(-ntml:ntm)                                   !complex
       complex phig(-ntml:ntm)                                 !complex
       complex apar(-ntml:ntm)                                 !complex
       complex bpar(-ntml:ntm)                                 !complex
       complex aparg(-ntml:ntm)                               !complex
       complex bparg(-ntml:ntm)                               !complex
       complex antot(-ntml:ntm)                                 !complex
       complex antota(-ntml:ntm)                               !complex
       complex antotb(-ntml:ntm)                               !complex
       complex bpar_mhd(-ntml:ntm)                             !complex
       complex term(negmax,nspecmax)                             !complex
       complex termi(negmax,nspecmax)                            !complex
       complex geint(-ntml:ntm,nspecmax)                          !complex
       complex g0eint(-ntml:ntm,nspecmax)                         !complex
       complex g1eint(-ntml:ntm,nspecmax)                         !complex
       complex ginti(-ntml:ntm,negmax,nspecmax)                   !complex
       complex gint(-ntml:ntm,negmax,nspecmax)                    !complex
       complex g1int(-ntml:ntm,negmax,nspecmax)                   !complex
       complex g2int(-ntml:ntm,negmax,nspecmax)                   !complex
       complex g3int(-ntml:ntm,negmax,nspecmax)                   !complex
       complex r(-ntml:ntm,nlm,negmax,nspecmax)                   !complex
       complex ainv(-ntml:ntm,nlm,negmax,nspecmax)                !complex
       complex a(-ntml:ntm,nlm,negmax,nspecmax)                   !complex
       complex b(-ntml:ntm,nlm,negmax,nspecmax)                   !complex
       complex g(-ntml:ntm,nlm,negmax,2,nspecmax)                 !complex
       complex gnew(-ntml:ntm,nlm,negmax,2,nspecmax)              !complex
       complex g0(-ntml:ntm,nlm,negmax,2,nspecmax)                !complex
       complex g1(-ntml:ntm,nlm,negmax,2,nspecmax)                !complex
       complex g2(-ntml:ntm,nlm,negmax,2,nspecmax)                !complex
       complex am(ntmp2p1t3,ntmp2p1t3)                        !complex
       complex transfer                                    !complex
       complex beta1                                       !complex
       complex add                                          !complex
       complex zi                                              !complex
       complex omega                                     !complex
       complex omegaprev                              !complex
       complex omega0                                    !complex
       complex omm1                                        !complex
       complex alphal                                       !complex
       complex alphapl                                     !complex
       complex ala                                             !complex
       complex al1a                                           !complex
       complex omegacom                              !complex
       complex c                                                 !complex
       complex d(nlmt2)                                    !complex
       complex x(nlmt2)                                    !complex
       complex delta(nlmt2)                             !complex
       complex alpha                                        !complex
       complex  alpha1                                     !complex
       complex  alp1                                          !complex
       complex  al1a1                                        !complex
       complex  phim(-ntml:ntm)                                  !complex
       complex  aparm(-ntml:ntm)                                 !complex
       complex  bparm(-ntml:ntm)                                 !complex
       complex  phinorm(-ntml:ntm,0:nkgrid)                      !complex
       complex  aparnorm(-ntml:ntm,0:nkgrid)                     !complex
       complex  bparnorm(-ntml:ntm,0:nkgrid)                     !complex
       complex  weightchkt(-ntml:ntm)                            !complex
       complex  weightchkp(-ntml:ntm)                            !complex
       complex  weightchk(-ntml:ntm)                             !complex
       complex  colratedet(-ntml:ntm,nlm,negmax,2)               !complex
       complex  colratedetm(-ntml:ntm,nlm,negmax,2)              !complex
       complex  colratedetcm(-ntml:ntm,nlm,negmax,2)             !complex
       complex  colratetrap(-ntml:ntm)                           !complex
       complex  colratepass(-ntml:ntm)                           !complex
       complex  colratetot(-ntml:ntm)                            !complex
       complex  colratetrapn(-ntml:ntm)                          !complex
       complex  colratepassn(-ntml:ntm)                          !complex
       complex  colratetotn(-ntml:ntm)                           !complex
       complex  colratetrapm(-ntml:ntm)                          !complex
       complex  colratepassm(-ntml:ntm)                          !complex
       complex  colratetrapnm(-ntml:ntm)                         !complex
       complex  colratepassnm(-ntml:ntm)                         !complex
       complex  colratetrapcm(-ntml:ntm)                         !complex
       complex  colratepasscm(-ntml:ntm)                         !complex
       complex  colratetrapncm(-ntml:ntm)                        !complex
       complex  colratepassncm(-ntml:ntm)                        !complex
       complex  respfunc(-ntml:ntm,nspecmax)                        !complex
       complex  dennorm(-ntml:ntm,nspecmax)                         !complex
       complex  dengyronorm(-ntml:ntm,nspecmax)                     !complex
       complex  phf(-ntm:ntm)                                           !complex
       complex  denf(-ntm:ntm,nspecmax)                     !complex
       complex  respff(-ntm:ntm,nspecmax)                   !complex
       complex  curtrap(0:ntime,-ntml:ntm)                     !complex
       complex  curpass(0:ntime,-ntml:ntm)                    !complex
       complex  wave(0:ntime)                                 !complex
       complex pmhd(0:ntm)                                    !complex
**********end complex variables
*
**********start real variables
       real x_bpar                                       !real
       real   y_bpar                                     !real
       real   xy_bpar                                   !real
!rew add signomega +1 ion freq > 0 flow to -x
!rew                 real 1          < 0         +x
       real signomega                                   !real
!rew icontinue=0(gt 0) do (don't) re-initialize to adiabatic
         real e(negmax)                                   !real
         real w(negmax)                                   !real
         real al(nlm)                                     !real
         real bmag(-ntml:ntm)                                 !real
         real theta(-ntml:ntm)                                !real
         real delthet(-ntml:ntm)                              !real
         real gamtot(-ntml:ntm)                               !real
         real gamtotpb(-ntml:ntm)                             !real
         real gamtota(-ntml:ntm)                              !real
         real gamtotb(-ntml:ntm)                              !real
         real gamtotbb(-ntml:ntm)                             !real
!rew added for quasilinear
         real poldrift(-ntml:ntm,nspecmax)                       !real     	
         real aintnorm(-ntml:ntm,nspecmax)                       !real
         real vnew(negmax,nspecmax)                          !real
         real wstar(negmax,nspecmax)                         !real
         real wl(-ntml:ntm,nlm)                               !real
         real sq(-ntml:ntm,nlm,2)                             !real
         real vpa(-ntml:ntm,nlm,negmax,2)                     !real
         real vpac(-ntml:ntm,nlm,negmax,2)                    !real
! rew 9/1/98 add vper**2
         real vper2(-ntml:ntm,nlm,negmax)                     !real
         real aj0(-ntml:ntm,nlm,negmax,nspecmax)                 !real
         real azj1(-ntml:ntm,nlm,negmax,nspecmax)                 !real
         real wdrift(-ntml:ntm,nlm,negmax,nspecmax)              !real
         real aa(-ntml:ntm,nlm,negmax,nspecmax)                  !real
         real bb(-ntml:ntm,nlm,negmax,nspecmax)                  !real
         real cc(-ntml:ntm,nlm,negmax,nspecmax)                  !real
         real a1(nlmt2)                                   !real
         real b1(nlmt2)                                   !real
         real c1(nlmt2)                                   !real
         real betaa(nlmt2)                                !real
         real ql(nlmt2)                                   !real
         real vpar(-ntml:ntm,nlm,negmax,2,nspecmax)              !real
         real bparophi_mhd(-ntml:ntm)                         !real
         real pi                                          !real
         real fprimm(nspecmax)                               !real
         real ecut                                        !real
! /mag/
         real bmin                                        !real
         real bmax                                        !real
         real eps                                         !real
         real shift                                       !real
!rew added
         real xwell                                       !real
! /plaspar/
         real vnew0(nspecmax)                             !real
         real shat                                        !real
         real aky                                         !real
         real pp                                          !real
         real uprim(nspecmax)                             !real
         real fprim(nspecmax)                             !real
         real tprim(nspecmax)                             !real
         real epsl                                        !real
         real beta                                        !real
         real gridfac                                     !real
!2/3/98 rew added debyelorhos
         real debyelorhos                                 !real
         real cdebye                                      !real
! /itpar/
         real delt                                        !real
         real fexp(nspecmax)                              !real
         real width0                                      !real
         real phydif                                      !real
         real power                                       !real
         real diff                                        !real
         real tcomp                                       !real
         real fv                                          !real
         real test1                                       !real
         real test2                                       !real
! /species/
         real z(nspecmax)                                 !real
         real amass(nspecmax)                             !real
         real temp(nspecmax)                      !real
         real an(nspecmax)                                !real
         real bakdif(nspecmax)                            !real
         real zeff                                        !real
         real teti                                        !real
!rew xteti fix
         real xteti                                       !real
         real dbeam                                       !real
         real zeff5                                       !real
! /omval/
        real absom                                       !real
        real pflux(nspecmax)                              !real
        real qheat(nspecmax)                              !real
        real adamp                                        !real
! /wstardat/
         real tguess1                                      !real
         real tguess3                                      !real
         real tp1                                          !real
         real tp3                                          !real
         real time                                         !real
         real gammause                                     !real
!rew  an1 left out?
         real an1                                          !real
         real an2                                          !real
         real an4                                          !real
         real an5                                          !real
         real akperp                                       !real
         real theta0                                       !real
         real bt                                           !real
         real ne                                           !real
         real ni                                           !real
         real nz                                           !real
         real nfast                                        !real
         real nfastp                                       !real
         real te                                           !real
         real ti                                           !real
         real q                                            !real
         real lne                                          !real
         real lni                                          !real
         real lnz                                          !real
         real lte                                          !real
         real lti                                          !real
         real ls                                           !real
         real rmin                                         !real
         real rsurf                                        !real
         real rmaj                                         !real
         real teprf0                                       !real
         real tiprf0                                       !real
         real neprf0                                       !real
         real fteprfa                                      !real
         real ftiprfa                                      !real
         real fneprfa                                      !real
         real alfprfte                                     !real
         real alfprfti                                     !real
         real alfprfne                                     !real
         real alfprfj                                      !real
         real qprf0                                        !real
         real qprfa                                        !real
         real qextra                                       !real
         real deltprf0                                     !real
         real deltprfgam                                   !real
         real deltprf(0:ipgrid)                                !real
         real teprf(0:ipgrid)                                  !real
         real tiprf(0:ipgrid)                                  !real
         real neprf(0:ipgrid)                                  !real
         real rlteprf(0:ipgrid)                                !real
         real rltiprf(0:ipgrid)                                !real
         real rlneprf(0:ipgrid)                                !real
         real qprf(0:ipgrid)                                   !real
         real rprf(0:ipgrid)                                   !real
         real shatprf(0:ipgrid)                                !real
         real zerop(0:ipgrid)                                  !real
         real betaprf(0:ipgrid)                                !real
         real vnewstareprf(0:ipgrid)                           !real
         real vnewstariprf(0:ipgrid)                           !real
         real tetiprf(0:ipgrid)                                !real
         real dgyrobohms(0:ipgrid)                             !real
         real dmxlsgyrobohm(0:ipgrid)                          !real
         real dmxlpsgyrobohm(0:ipgrid)                          !real
         real di_effcgyrobohm(0:ipgrid)                        !real
         real di_effsgyrobohm(0:ipgrid)                        !real
         real di_effwgyrobohm(0:ipgrid)                        !real
         real de_effcgyrobohm(0:ipgrid)                        !real
         real de_effsgyrobohm(0:ipgrid)                        !real
         real de_effwgyrobohm(0:ipgrid)                        !real
         real chii_effcgyrobohm(0:ipgrid)                      !real
         real chii_effsgyrobohm(0:ipgrid)                      !real
         real chii_effwgyrobohm(0:ipgrid)                      !real
         real chie_effcgyrobohm(0:ipgrid)                      !real
         real chie_effsgyrobohm(0:ipgrid)                      !real
         real chie_effwgyrobohm(0:ipgrid)                      !real
         real agammasprf(0:nkgrid,0:ipgrid)                    !real
         real dgammasprf(0:nkgrid,0:ipgrid)                    !real
         real afreqsprf(0:nkgrid,0:ipgrid)                     !real
         real dfreqsprf(0:nkgrid,0:ipgrid)                     !real
         real dmxlsprf(0:nkgrid,0:ipgrid)                      !real
         real dmxlpsprf(0:nkgrid,0:ipgrid)                      !real
         real di_effcsprf(0:nkgrid,0:ipgrid)                   !real
         real di_effssprf(0:nkgrid,0:ipgrid)                   !real
         real di_effwsprf(0:nkgrid,0:ipgrid)                   !real
         real de_effcsprf(0:nkgrid,0:ipgrid)                   !real
         real de_effssprf(0:nkgrid,0:ipgrid)                   !real
         real de_effwsprf(0:nkgrid,0:ipgrid)                   !real
         real chii_effcsprf(0:nkgrid,0:ipgrid)                   !real
         real chii_effssprf(0:nkgrid,0:ipgrid)                   !real
         real chii_effwsprf(0:nkgrid,0:ipgrid)                   !real
         real chie_effcsprf(0:nkgrid,0:ipgrid)                   !real
         real chie_effssprf(0:nkgrid,0:ipgrid)                   !real
         real chie_effwsprf(0:nkgrid,0:ipgrid)                   !real
         real gammadelprf(0:nkgrid,0:ipgrid)                   !real
         real freqdelprf(0:nkgrid,0:ipgrid)                    !real
         real tp1crtprf(0:nkgrid,0:ipgrid)                     !real
         real tp3crtprf(0:nkgrid,0:ipgrid)                     !real
         real maxagammasprf(0:ipgrid)                          !real
         real maxdgammasprf(0:ipgrid)                          !real
         real maxafreqsprf(0:ipgrid)                           !real
         real maxdmxlsprf(0:ipgrid)                            !real
         real maxdmxlpsprf(0:ipgrid)                            !real
         real maxdi_effcsprf(0:ipgrid)                         !real
         real maxdi_effssprf(0:ipgrid)                         !real
         real maxdi_effwsprf(0:ipgrid)                         !real
         real maxde_effcsprf(0:ipgrid)                         !real
         real maxde_effssprf(0:ipgrid)                         !real
         real maxde_effwsprf(0:ipgrid)                         !real
         real maxchii_effcsprf(0:ipgrid)                       !real
         real maxchii_effssprf(0:ipgrid)                       !real
         real maxchii_effwsprf(0:ipgrid)                       !real
         real maxchie_effcsprf(0:ipgrid)                       !real
         real maxchie_effssprf(0:ipgrid)                       !real
         real maxchie_effwsprf(0:ipgrid)                       !real
         real maxagammasprfprev(0:ipgrid)                      !real
         real maxdgammasprfprev(0:ipgrid)                      !real
********* Basiscom:
       real dgyrobohmnorm                               !real
       real  t(3,3)                                       !real
       real  ratio(3,3)                                   !real
       real  akyarr(16)                                   !real
       real  kysi                                         !real
       real  kys(0:nkgrid)                                !real
       real  freqsj(-ntml:ntm,0:nkgrid)                       !real
       real  gammasj(-ntml:ntm,0:nkgrid)                      !real
       real  afreqs(0:nkgrid)                             !real
       real  dfreqs(0:nkgrid)                             !real
       real  agammas(0:nkgrid)                            !real
       real  dgammas(0:nkgrid)                            !real
       real  dtgammas(0:nkgrid)                           !real
       real  afreq(0:nkgrid)                              !real
       real  agamma(0:nkgrid)                             !real
       real  agammasprev(0:nkgrid)                        !real
       real  dgammasprev(0:nkgrid)                        !real
       real  gammadel(0:nkgrid)                             !real
       real  freqdel(0:nkgrid)                            !real
       real  twidths(0:nkgrid)                            !real
       real  dmxls(0:nkgrid)                              !real
       real  dmxlps(0:nkgrid)                              !real
       real  di_effc(0:nkgrid)                            !real
       real  di_effs(0:nkgrid)                            !real
       real  di_effw(0:nkgrid)                            !real
       real  de_effc(0:nkgrid)                            !real
       real  de_effs(0:nkgrid)                            !real
       real  de_effw(0:nkgrid)                            !real
       real  chii_effc(0:nkgrid)                            !real
       real  chii_effs(0:nkgrid)                            !real
       real  chii_effw(0:nkgrid)                            !real
       real  chie_effc(0:nkgrid)                            !real
       real  chie_effs(0:nkgrid)                            !real
       real  chie_effw(0:nkgrid)                            !real
       real  di_effgbc(0:nkgrid)                            !real
       real  di_effgbs(0:nkgrid)                            !real
       real  di_effgbw(0:nkgrid)                            !real
       real  de_effgbc(0:nkgrid)                            !real
       real  de_effgbs(0:nkgrid)                            !real
       real  de_effgbw(0:nkgrid)                            !real
       real  chii_effgbc(0:nkgrid)                            !real
       real  chii_effgbs(0:nkgrid)                            !real
       real  chii_effgbw(0:nkgrid)                            !real
       real  chie_effgbc(0:nkgrid)                            !real
       real  chie_effgbs(0:nkgrid)                            !real
       real  chie_effgbw(0:nkgrid)                            !real
       real  tp1crt(0:nkgrid)                             !real
       real  tp3crt(0:nkgrid)                             !real
       real  freqst(0:ntime,0:20,nkgrid)                    !real
       real  gammast(0:ntime,0:20,nkgrid)                   !real
       real  timest(0:ntime)                                !real
       real  tp1_rec(0:ntime,2)                               !real
       real  tp3_rec(0:ntime,2)                               !real
       real  slbnorm(1:nlm,-ntml:ntm)                         !real
       real  vnewstare                                    !real
       real  vnewstari                                    !real
       real  stkx(0:20)                                   !real
       real  stkafreqs(0:nkgrid,0:20)                     !real
       real  stkagammas(0:nkgrid,0:20)                    !real
       real  stkgammast(0:ntime,0:20,nkgrid)                !real
       real  stkdfreqs(0:nkgrid,0:20)                     !real
       real  stkdgammas(0:nkgrid,0:20)                    !real
       real  stktwidths(0:nkgrid,0:20)                    !real
       real  stkdmxls(0:nkgrid,0:20)                      !real
       real  stkdmxlps(0:nkgrid,0:20)                      !real
       real  stkdi_effc(0:nkgrid,0:20)                    !real
       real  stkdi_effs(0:nkgrid,0:20)                    !real
       real  stkdi_effw(0:nkgrid,0:20)                    !real
       real  stkde_effc(0:nkgrid,0:20)                    !real
       real  stkde_effs(0:nkgrid,0:20)                    !real
       real  stkde_effw(0:nkgrid,0:20)                    !real
       real  stkchii_effc(0:nkgrid,0:20)                  !real
       real  stkchii_effs(0:nkgrid,0:20)                  !real
       real  stkchii_effw(0:nkgrid,0:20)                  !real
       real  stkchie_effc(0:nkgrid,0:20)                  !real
       real  stkchie_effs(0:nkgrid,0:20)                  !real
       real  stkchie_effw(0:nkgrid,0:20)                  !real
       real  stkdi_effgbc(0:nkgrid,0:20)                    !real
       real  stkdi_effgbs(0:nkgrid,0:20)                    !real
       real  stkdi_effgbw(0:nkgrid,0:20)                    !real
       real  stkde_effgbc(0:nkgrid,0:20)                    !real
       real  stkde_effgbs(0:nkgrid,0:20)                    !real
       real  stkde_effgbw(0:nkgrid,0:20)                    !real
       real  stkchii_effgbc(0:nkgrid,0:20)                  !real
       real  stkchii_effgbs(0:nkgrid,0:20)                  !real
       real  stkchii_effgbw(0:nkgrid,0:20)                  !real
       real  stkchie_effgbc(0:nkgrid,0:20)                  !real
       real  stkchie_effgbs(0:nkgrid,0:20)                  !real
       real  stkchie_effgbw(0:nkgrid,0:20)                  !real
       real  stkphinorm(-ntml:ntm,nkgrid,0:20)             !real
       real  stktp1crt(0:nkgrid,0:20)                     !real
       real  stktp3crt(0:nkgrid,0:20)                     !real
       real  park(-ntm:ntm)                          !real
       real  pflxe(-ntml:ntm,nspecmax)                           !real
       real  pflxm(-ntml:ntm,nspecmax)                           !real
       real  pflxte(-ntml:ntm,nspecmax)                          !real
       real  pflxtm(-ntml:ntm,nspecmax)                          !real
       real  eflxe(-ntml:ntm,nspecmax)                           !real
       real  eflxm(-ntml:ntm,nspecmax)                           !real
       real  eflxte(-ntml:ntm,nspecmax)                          !real
       real  eflxtm(-ntml:ntm,nspecmax)                          !real
       real  mflxe(-ntml:ntm,nspecmax)                          !real
       real  pflxea(nspecmax)                                !real
       real  pflxma(nspecmax)                                !real
       real  pflxtea(nspecmax)                               !real
       real  pflxtma(nspecmax)                               !real
       real  eflxea(nspecmax)                                !real
       real  eflxma(nspecmax)                                !real
       real  eflxtea(nspecmax)                               !real
       real  eflxtma(nspecmax)                               !real
       real  mflxea(nspecmax)                                !real
       real  peflx(-ntml:ntm)                                 !real
       real  eeflx(-ntml:ntm)                                 !real
       real  piflx(-ntml:ntm)                                 !real
       real  eiflx(-ntml:ntm)                                 !real
       real  peflxa                                       !real
       real  eeflxa                                       !real
       real  piflxa                                       !real
       real  eiflxa                                       !real
       real  cwave                                        !real
       real ominst           !real
       real  tol              !real
       real  gamma            !real
       real  fcv              !real
       real  fexp1            !real
       real  fexp3            !real
       real  thetamin         !real
       real  thetamax         !real
       real  alr              !real
       real  ali              !real
       real  al1r             !real
       real  al1i             !real
       real  alar             !real
       real  alai             !real
       real  al1ar            !real
       real  al1ai            !real
       real  aky1             !real
       real  aky2             !real
       real  aky3             !real
       real  aky4             !real
       real  aky5             !real
       real  aky6             !real
       real  aky7             !real
       real  aky8             !real
       real  aky9             !real
       real  aky10            !real
       real  aky11            !real
       real  aky12            !real
       real  aky13            !real
       real  aky14            !real
       real  aky15            !real
       real  aky16            !real
        real pk               !real
! rew add
       real  epsa             !real
       real  mach1            !real
       real  mach2            !real
       real  mach3            !real
       real  mach4            !real
       real  mach5            !real
       real  uprim1           !real
       real  uprim2           !real
       real  uprim3           !real
       real  uprim4           !real
       real  uprim5           !real
       real  egamma           !real       
       real  fprim1           !real
       real  fprim2           !real
       real  fprim3           !real
       real  fprim4           !real
       real  fprim5           !real
       real  tprim1           !real
       real  tprim2           !real
       real  tprim3           !real
       real  tprim4           !real
       real  tprim5           !real
       real  vnewk1           !real
       real  vnewk2           !real
       real  vnewk3           !real
       real  vnewk4           !real
       real  vnewk5           !real
       real  cnewk2           !real
       real  cnewk3           !real
       real  bakdif1          !real
       real  bakdif2          !real
       real  bakdif3          !real
       real  bakdif4          !real
       real  bakdif5          !real
       real  amass2           !real
       real  amass3           !real
       real  amass4           !real
       real  amass5           !real
       real  temp2            !real
       real  temp3            !real
       real  temp4            !real
       real  temp5            !real
       real  z2               !real
       real  z3               !real
       real  z4               !real
       real  z5               !real
********* Geo:
!Group
       real  aspectratio_loc  !real
       real  rmin_loc         !real
       real  shift_loc        !real
       real  delta_loc        !real
       real  kappa_loc        !real
       real  s_kappa_loc      !real
       real  s_delta_loc      !real
       real  q_loc            !real
       real  shat_loc         !real
       real  tiote_loc       !real
       real  nione_loc       !real
       real  dlntidr_loc     !real
       real  dlntedr_loc     !real
       real  dlnnidr_loc     !real
       real  dlnnedr_loc     !real
       real  beta_loc        !real
       real  xnu_loc         !real
       real  zeff_loc        !real
       real  dlnnimpdr_loc   !real
       real  fastionfrac_loc !real
! overrides
       real  ne_loc          !real
       real  te_loc          !real
       real  b00_loc         !real
       real  beta_loc_0      !real
       real  rmaj_mag_center !real
       real  dlnpdr_loc      !real
       real  q_prime_loc     !real
       real  p_prime_loc     !real
!outputs to gstotal from ifwritegeo
       real  b_unit                     !real
       real  a_unit                     !real
       real  b_norm                     !real
       real  b2_ave_geo                 !real
       real  epsl_geo(-ntml:ntm)        !real
       real  costheta_geo(-ntml:ntm)    !real
       real  costheta_p_geo(-ntml:ntm)    !real
       real  sintheta_geo(-ntml:ntm)    !real
       real  kxoky_geo(-ntml:ntm)       !real
       real  pk_geo(-ntml:ntm)          !real
       real  b_geo(-ntml:ntm)           !real
       real  qrat_geo(-ntml:ntm)        !real
       real  ppgeo(-ntml:ntm)          !real
       real  theta_c(-ntml:ntm)   !real
       real  dtheta_cdtheta(-ntml:ntm) !real
       real  epsl_geo_c(-ntml:ntm)        !real
       real  costheta_geo_c(-ntml:ntm)    !real
       real  costheta_p_geo_c(-ntml:ntm)    !real
       real  sintheta_geo_c(-ntml:ntm)    !real
       real  kxoky_geo_c(-ntml:ntm)       !real
       real  pk_geo_c(-ntml:ntm)          !real
       real  b_geo_c(-ntml:ntm)           !real
       real  qrat_geo_c(-ntml:ntm)        !real
       real  p_prime_zero               !real
!diagnostic
       real  beta_loc_out              !real
       real  alpha_mhd_loc_out     !real
       real  shat_loc_out              !real
       real  shat_mhd_loc_out     !real
       real  xmu_loc_out            !real
       real  delta_dor                !real
       real  d_prime_dor              !real
       real  k_prime_dor              !real
       real  volume_loc               !real
       real  volume_c                 !real
       real  rmaj0_loc                !real
       real  shat_mhd_loc             !real
       real  alpha_mhd_loc            !real
       real  alpha_c                  !real
       real  kxoky_0_geo(-ntml:ntm)       !real
       real  kxoky_p_geo(-ntml:ntm)       !real
       real  kxoky_s_geo(-ntml:ntm)       !real
       real  rmaj_theta(0:ntm)         !real
       real  z_theta(0:ntm)            !real
       real  r_curv(0:ntm)             !real
       real  z_l(0:ntm)                !real
       real  z_ll(0:ntm)               !real
       real  rmaj_l(0:ntm)             !real
       real  rmaj_ll(0:ntm)            !real
       real  grad_r_theta(0:ntm)       !real
       real  bt_theta(0:ntm)           !real
       real  bp_theta(0:ntm)           !real
       real  b_theta(0:ntm)            !real
       real  dl1dtheta(0:ntm)          !real
       real  d_0(0:ntm)                !real
       real  d_p(0:ntm)                !real
       real  d_ffp(0:ntm)              !real
       real  s1_check(0:ntm)           !real
       real  dbpdrho1(0:ntm)           !real
       real  dbpdrho2(0:ntm)           !real
       real  dbpdrho3(0:ntm)           !real
       real  theta_bar(0:ntm)           !real
       real  amhd(0:ntm)                !real
       real  bmhd(0:ntm)                !real
       real  cmhd(0:ntm)                !real
!1/28/98 rew add mhd_crit routine
       real  s_mhd               !real
       real  delt_mhd            !real
       real  eps_mhd             !real
       real  omega_mhd_ave(1:4000) !real
       real  omega_mhd_dev        !real
       real  gamma_mhd_ave(1:4000) !real
       real  gamma_mhd_ave_fin
       real  gamma_mhd_dev        !real
       real  damp_mhd             !real
       real  xmu_mhd              !real
       real y_mhd                  ! real
********* Share:
       real rmajor_s              ! real
       real epsilon_s             ! real
       real q_s                   ! real
       real shat_s                ! real
       real betae_s               ! real
       real xnu_s                 ! real
       real taui0_s               ! real
       real rlte_s                ! real
       real rlti_s                ! real
       real rln_s                 ! real
       real dt_s                  ! real
       real anrate_s              ! real
       real dnrate_s              ! real
       real dtnrate_s             ! real
       real anfreq_s              ! real
       real dnfreq_s              ! real
       real chie_s                ! real
       real chii_s                ! real
       real diff_s                ! real
       real diffmxl_s             ! real
       real deltheta_s            ! real
********* Sharegks:
       real bt_s                  ! real
       real zeff_s                ! real
       real rsurf_s               ! real
       real ne_s                  ! real
       real ni_s                  ! real
       real nz_s                  ! real
       real nfast_s               ! real
       real te_s                  ! real
       real ti_s                  ! real
       real qgks_s                ! real
       real lne_s                 ! real
       real lni_s                 ! real
       real lnz_s                 ! real
       real lte_s                 ! real
       real lti_s                 ! real
       real ls_s                  ! real
       real amass3_s              ! real
       real rmin_s                ! real
       real rmaj_s                ! real
       real aky1_s                ! real
       real shift_s               ! real
       real uprim1_s              ! real
       real uprim3_s              ! real
       real delt_s                ! real
       real agammas_s             ! real
       real dgammas_s             ! real
       real dtgammas_s            ! real
       real afreqs_s              ! real
       real dfreqs_s              ! real
       real chie_eff_s            ! real
       real chii_eff_s            ! real
       real diff_eff_s            ! real
******* end real variables
*
*******start integer variables
        integer igyro_fix                  !integer
        integer igyro_e                    !integer
! =0 old =1 fixes bessel function missing 1/B factor in original code
        integer  i_solve                           !integer
! NAG solver in place of gauss
!1/21/98 rew added i_pbar=1 option
! added bpar in analogy to apar,  and antotb like antota, gamtotb like gamtot
! i_bpar=0 old no B_par version
! i_bpar=1 new B_par version with x_bpar=0 should give bpar=0 result
!     x_bpar should be =1. for physical B_par not 0 result
       integer i_bpar                                      !integer
       integer  icontinue                             !integer
!rew _sav added for icontinue=1
        integer jend(-ntml:ntm)                                 !integer
        integer ngauss                                      !integer
        integer i_mhd                                            !integer
        integer ng2                                         !integer
        integer nlambda                                     !integer
        integer ntgrid                                      !integer
        integer ipar                                        !integer
        integer ioddpar                                     !integer
        integer ntgridl                                     !integer
        integer ntheta                                      !integer
        integer ntheth                                      !integer
        integer nthethl                                     !integer
        integer nperiod                                     !integer
! nperiod1 replaced with nperiod11 and nperiod10
       integer nperiod11                                   !integer
       integer nperiod10                                   !integer
       integer nperiod2                                    !integer
       integer ntgridr              !integer
       integer  nstep_mhd           !integer
       integer  iptot           !integer
       integer  igeo                     !integer
       integer  igeot                    !integer
       integer  igeo_print               !integer
       integer negrid                                      !integer
       integer icv                                         !integer
       integer nspec                                !integer
       integer ncspec(nspecmax)                            !integer
       integer nce(negmax)                                 !integer
       integer ivnewout                              !integer
       integer iaky                                         !integer
       integer igam                                         !integer
       integer istep                                        !integer
       integer ict                                          !integer
       integer isvar                                        !integer
       integer nwrite                                       !integer
       integer iphi                                         !integer
       integer nstep                                        !integer
!rew added nteststep, adamp
       integer nteststep                                    !integer
       integer icon                                         !integer
       integer ipmaxmax                                     !integer
       integer  ipmax                                        !integer
       integer  ipmin                                        !integer
       integer istepkprf(0:nkgrid,0:ipgrid)                     !integer
       integer iakymax(0:ipgrid)                                !integer
       integer ngksout                                      !integer
       integer  nsteps(0:nkgrid)                             !integer
       integer  noread                                       !integer
       integer  ikys                                         !integer
       integer  istepk(0:nkgrid)                             !integer
       integer  iistk                                         !integer
       integer  stknstep(0:nkgrid,0:20)                      !integer
!reads 10
!___________________________________________________
       integer  nscreen          !integer
       integer  nout             !integer
       integer  ngamstep       !integer
       integer ncspec1           !integer
       integer  ncspec2          !integer
       integer  ncspec3          !integer
       integer  ncspec4          !integer
       integer  ncspec5          !integer
       integer  nce1                !integer
       integer  nce2                !integer
       integer  nce3                !integer
       integer  nce4                !integer
       integer  nce5                !integer
       integer  naky                !integer
       integer  nt0              !integer
       integer  nt0min           !integer
       integer  nt0max           !integer
****** end of integer variables
*
******
c     external f04ade, f01aae
      contains
!
      include 'gks_gstotal.f'

      end module gks_var

