ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filter_dv
c
c     Use a sixth order filer to remove the 2dx mode from Tnew
c
      implicit none
c
      include '../inc/ptor.m'
      include '../inc/vtrans.m'
      real*8 f(1:ngrid+2)
      integer i,j
c
      do j=1,nfields
        do i=4,ngrid-4
         f(i) = Tnew(j,i)+(Tnew(j,i+3)-6.D0*Tnew(j,i+2)
     >   +15.D0*Tnew(j,i+1)-20.D0*Tnew(j,i)
     >   +15.D0*Tnew(j,i-1)-6.D0*Tnew(j,i-2)+Tnew(j,i-3))/64.D0
        enddo
c  extend grid to 0,-1,-2 with flat profile using i=1 value (zero gradient BC)
        i=3
        f(i) = Tnew(j,i)
     >  +(Tnew(j,i+3)-6.D0*Tnew(j,i+2)+15.D0*Tnew(j,i+1)-20.D0*Tnew(j,i)
     >      +15.D0*Tnew(j,i-1)-6.D0*Tnew(j,1)+Tnew(j,1))/64.D0
        i=2
        f(i) = Tnew(j,i)
     >  +(Tnew(j,i+3)-6.D0*Tnew(j,i+2)+15.D0*Tnew(j,i+1)-20.D0*Tnew(j,i)
     >    +15.D0*Tnew(j,1)-6.D0*Tnew(j,1)+Tnew(j,1))/64.D0
        i=1
        f(i) = Tnew(j,i)
     >  +(Tnew(j,i+3)-6.D0*Tnew(j,i+2)+15.D0*Tnew(j,i+1)-20.D0*Tnew(j,1)
     >    +15.D0*Tnew(j,1)-6.D0*Tnew(j,1)+Tnew(j,1))/64.D0
c  extend grid to ngrid,ngrid+1,ngrid+2 using Tnew(ngrid)=0 BC.
        i=ngrid-3
         f(i) = Tnew(j,i)
     >    + (-6.D0*Tnew(j,i+2)+15.D0*Tnew(j,i+1)-20.D0*Tnew(j,i)
     >      +15.D0*Tnew(j,i-1)-6.D0*Tnew(j,i-2)+Tnew(j,i-3))/64.D0
        i=ngrid-2
         f(i) = Tnew(j,i)
     >    + (+15.D0*Tnew(j,i+1)-20.D0*Tnew(j,i)
     >      +15.D0*Tnew(j,i-1)-6.D0*Tnew(j,i-2)+Tnew(j,i-3))/64.D0
        i=ngrid-1
         f(i) = Tnew(j,i)
     >    + (-20.D0*Tnew(j,i)
     >      +15.D0*Tnew(j,i-1)-6.D0*Tnew(j,i-2)+Tnew(j,i-3))/64.D0
        do i=1,ngrid-1
         Tnew(j,i) = f(i)
        enddo
      enddo
      return
      end
