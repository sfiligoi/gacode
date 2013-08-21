!
      MODULE  ifs
        USE param, ONLY : kj,kjm1
        USE mhdpar,ONLY : kpsi
        IMPLICIT NONE
!
!     for use in the IFS (i.e., Dorland-Kotchenreuther) model ... HSJ
!
      REAL *8                                                            &
                  rhod_ifs(kj), rlt_ifs(kj), rln_ifs(kj),                &
                  rlne_ifs(kj), shat_ifs(kj), enb_ifs(kj),               &
                  tau_ifs(kj), eps_ifs(kj), zth_ifs(kj),                 &
                  rhod_psi_ifs(kpsi), chi_i_ifs(kj),                     &
                  chi_e_ifs(kj), RLTcrit_ifs(kj),dmassden(kjm1),         &
                  RLTcritz_ifs(kj), gnu_ifs(kj), d_ifs(kj),              &
                  vthi_ifs(kj), omci_ifs(kj), dorl_kotch,                &
                  dorl_kotche,dorl_kotchi,exbmult_ifs,                   &
                  xke_ifs(kj), xki_ifs(kj), xkang_ifs(kj),               &
                  g_ifs(kj), grho1_ifs(kj), grho2_ifs(kj),               &
                  rhod_max_ifs, time_ifs,                                &
                  wneo_elct_ifs
     INTEGER *4   include_ifs                             
!
!     this confinement model is used if
!              include_ifs = 1  (analysis and simulation)
!     and
!              dorl_kotch > 0.0 (multiplier in simulation only)
!     or
!              dorl_kotche,i > 0.0
!
!     rhod_ifs is the average horizontal minor radius in cm
!              it is interpolated onto the transport
!              rho grid  [ r(j=1,..nj) ] in subroutine RHOSET
!              after it is first calculated on the npsi mhd grid
!              in subroutine FLUXAV (rhod_psi_ifs).
!              rhod_ifs is normalized in rhoaset as well.
!    other parameters for the IFS model are described in subroutine
!    IP_CHI2. The routine IP2_CHI which evaluates this model is
!    called from DIFFUS.
!
     END MODULE IFS 
