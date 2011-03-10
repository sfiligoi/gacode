!COMM30-toroidal geometry - MHD (VMEC)
!  fiota(1,i)-iota-bar, 1/q (-)
!  fiota(2-4,i)-cubic spline coefficients for bar_iota (-)
!  rhm3d-radial grid normalized to a0*sqrt(tor flux) (m)
!  xm3d-radial grid normalized to sqrt(tor flux) (-)
!  elmna(1,i)-lambda coefficients normalized to xm3d**m3d (-)
!  elmna(2-4,i)-cubic spline coefficients for elmna (-)
!  rmna(1,i)-R coefficients normalized to xm3d**m3d (m)
!  rmna(2-4,i)-cubic spline coefficients for rmna (m)
!  zmna(1,i)-Z coefficients normalized to xm3d**m3d (m)
!  zmna(2-4,i)-cubic spline coefficients for zmna (m)
!  phitot-total toroidal flux (Wb)
!  km1n0-mode number for m=1,n=0 component (-)
!  nk3d-number of modes for lambda expansion (-)
!  nk3drz-number of modes for R and Z expansions (-)
!  nx3d-number of radial nodes (-)
!  m3d-poloidal mode numbers (-)
!  n3d-toroidal mode numbers (-)
      REAL           rhm3d,                   xm3d,
     &               fiota,                   elmna,
     &               rmna,                    zmna,
     &               phitot
      INTEGER        km1n0,                   nk3d,
     &               nk3drz,                  nx3d,
     &               m3d,                     n3d
      COMMON/COMM30/ rhm3d(mxnx),             xm3d(mxnx),
     &               fiota(4,mxnx),           elmna(4,mxnx,mxkm),
     &               rmna(4,mxnx,mxkm),       zmna(4,mxnx,mxkm),
     &               phitot,
     &               km1n0,                   nk3d,
     &               nk3drz,                  nx3d,
     &               m3d(mxkm),               n3d(mxkm)
     
