c
c --- parameters for fitting psi on rectangular boundary
c --- such that psi on plasma surface is constant
c
      logical    splinefitbd, fitboundary
      integer    nfitmax, mplasmax, nsymfit, nworkzx, mresidual,
     .           nfitpoints
      parameter (nfitmax  = 129)  ! max number of boundary values to fit
      parameter (mplasmax = 300)  ! max number of points on plasma
c                                   boundary to use in fit
      parameter (nsymfit  = ((nfitmax + 1)*nfitmax)/2)
      parameter (nworkzx  = 5*nfitmax + 2*mplasmax + nsymfit)
c
      common /fitparms/ mresidual, nfitpoints, splinefitbd, fitboundary
c
