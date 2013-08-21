  MODULE bicube
! --- cspln is used to store the bicubic spline coefficents for psi
! --- wnoperm is a work array. it is assumed that wnoperm contains no
! --- permanent information. it may be used for arbitrary but temporary
! --- storage only. (see INCLUDE file storage.i for additional arrays.)
!
      USE mhdpar,   only : nh,nwh,nw


      IMPLICIT NONE
      INTEGER,PARAMETER :: nh2     = 2*nh
      INTEGER,PARAMETER :: nwork   = 2*nwh+nh2
!     parameter      (nwork   = 2*nwh+2*nw)        ! use this if nw > nh
      INTEGER,PARAMETER  ::  n2cspln = 2

      REAL *8 cspln(n2cspln,nw,nh2),wnoperm(nwork),pds(6)

  END MODULE bicube
