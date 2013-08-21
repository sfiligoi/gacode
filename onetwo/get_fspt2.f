      subroutine TorGA_get_fspt2 (psi1d, farray, npsi1d, fspt2)
cProlog

c
      implicit none
c
      integer npsi1d
      doubleprecision psi1d(*), farray(*), fspt2(*)
c
c     uses subroutine CSPLINE
c
c     local variables:
c
      doubleprecision yp1, yp2
c use natural b.c.
      data yp1, yp2 /1.0d30, 1.0d30/ 
c
      call cspline (psi1d, farray, npsi1d, yp1, yp2, fspt2)
      return
c
      end
