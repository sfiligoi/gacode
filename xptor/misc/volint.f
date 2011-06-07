      real*8 function volint(Apro)
c
c Calculate the volume integral of some profile APro
c
      implicit none
c
c
      include '../inc/ptor.m'
c
      real*8 Apro(mxgrid)
      real*8 sum
      integer k
c
      sum=0.D0
      do k=1,ngrid
         sum=sum+Apro(k)*vprime(k,1)*dr(k,1)
c         write(*,*) k, Apro(k), sum
      enddo
c
      volint=sum
      return
      end
