

   MODULE cer

     USE param,  ONLY : kj,ksplin,kbctim
     USE mhdpar, ONLY : kpsi
     USE nrtype, ONLY : DP ,I4B
! **********************************************************************
!                                                                      *
!   Common blocks for CER velocity data, magnetic components           *
!   and computed velocities                                            *
!                                                                      *
!   Fuctions computed in NCLBOOT on the full grid                      *
!                                                                      *
!   Epsi(kj)         Er/(R Bpol) V/M^2/T flux function                 *
!   Epsi_exp(kj)     Er/(RBpol) from cer_ion data                      *
!   ave_vpar_d(kj)   main ion flux averaged parallel velocity m/sec    *
!   Kpol_d(kj)       main ion poloidal velocity/Bpol                   *
!   Kpol_c(kj)       CER  ion poloidal velocity/Bpol                   *
!   angrot_d(kj)     main ion toroidal velocity/R   along Z=0 cord     *
!   udia_d(kj)       main ion diamagnetic velocity/(R*Bpol)            *
!   udia_c(kj)       CER ion  diamagnetic velocity/(R*Bpol)            *
!   ugrt(kj)         ion temperature gradient/(e*R*Bpol)               *
!   sqz_d(jk)        main ion orbit squeeze factor = 1 - dVexb/dr/wc   *
!                    where wc = ion cyclotron frequency                *
!                    and Vexb = ExB velocity                           *
!                                                                      *
!   Fuctions computed in FLUXAV on the full grid                       *
!                                                                      *
!   cer_bp  (kj)     Bp   along Z=0 cord                               *
!   cer_btdr(kj)     Bt/R along Z=0 cord                               *
!                                                                      *
!   DATA read from namelis1 input file                                 *
!                                                                      *
!   angrot_c(kj)     CER ion toroidal veloctiy/Rmajor along Z=0        *
!                    input as angrotin                                 *
!   Kpol_exp(kj)     CER ion poloidal velocity/cer_bp on full grid     *
!                                                                      *
!   kpolin(ksplin,kbctim)  spline knot values of Kpol_exp              *
!                                                                      *
!   rkpol(ksplin,kbctim)   spline knot location for kpolin             *
!                                                                      *
!   knotskpol(kbctim)      number of knots for kpolin                  *
!                                                                      *
!   Control variables                                                  *
!                                                                      *
!   j_cer            the ion index for the cer ion                     *
!                                                                      *
!   ikpol            flag to read kpolin (default = 0)                 *
!                                                                      *
!   Data Created  5/17/96  by G. M. Staebler                           *
!   Last Modified 8/16/96                                              *
!                                                                      *
! **********************************************************************
!
      INTEGER(I4B)         j_cer
      INTEGER(I4B)  knotskpol      
!
      REAL(DP)        cer_bp(kpsi), cer_btdr(kpsi), sqz_d(kj),     &
                      ave_vpar_d(kj), Kpol_d(kj), Kpol_c(kj),      &
                      Epsi(kj), angrot_d(kj), angrot_c(kj),        &
                      Epsi_exp(kj),Kpol_exp(kj), udia_d(kj),       &
                      udia_c(kj), ugrt(kj), kpolin(ksplin,kbctim), &
                      rkpol(ksplin,kbctim)

!
  END MODULE cer
