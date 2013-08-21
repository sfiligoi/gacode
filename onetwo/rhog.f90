  MODULE rhog
      USE param,                             ONLY :kj
      USE nrtype,                            ONLY :I4B,DP 
      IMPLICIT NONE
!
! --- rhog saves some results (..old) from previous MHD solution,
! --- in case they must be restored,
! --- and holds some MHD-related quantitites (on the transport grid,
! --- i.e., RHOGrid)
! --- rmajorvec is rmajor at elevation of magnetic axis on rho grid
! --- (which corresponds also to psir grid).
! --- bprmaj is value of poloidal B field at rmajorvec(j)
!
      REAL(DP),public ::                                                &
                   curold(kj), ffprim(kj), ffpold(kj), pprim (kj),      &
                   ppold (kj), press (kj), pressb   (kj), psir  (kj),   &
                   psiold(kj), rold  (kj), rmajorvec(kj), bprmaj(kj),   &
                   qpsir (kj), press_thermal(kj),                       &
                   rmajor_r(kj),rminor_r(kj),ravg_r(kj),                &
                   psir_at_t0(kj),btotrmaj(kj)
!
   END  MODULE rhog
