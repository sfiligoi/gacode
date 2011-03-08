      subroutine norm_dv(normS,normDT,normT,kmax)
c
      implicit none
c
      include '../inc/vtrans.m'
c
      integer i,k,kmax
      real*8 normS, normDT, normT(5)
      real*8 sumnormS,sumnormT,sumnormDT
      sumnormS=0.D0
      sumnormT=0.D0
      sumnormDT=0.D0
      do i=1,nfields
        normS=0.D0
        normDT=0.D0
        normT(i) = 0.D0
        do k=1,kmax-1
          normS = normS + s(i,k)*s(i,k)
          normDT = normDT + st(i,k)*st(i,k)
          normT(i) = normT(i) + (Tstart(i,k)*vrho3(i,i,k))**2
        enddo
        normT(i) = DMAX1(1.D0,DSQRT(normT(i)))
        sumnormS = sumnormS + DSQRT(normS)/normT(i)
        sumnormDT = sumnormDT + DSQRT(normDT)/normT(i)
      enddo
      normDT = sumnormDT
      normS = sumnormS
      return
      end
