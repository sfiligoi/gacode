ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function volintk(Apro,kgrid)
c
c Calculate the volume integral of some profile APro
c
      implicit none
c
      include '../inc/ptor.m'
c
      real*8 Apro(mxgrid)
      integer kgrid
      real*8 sum
      integer k
c
      sum=0.D0
      do k=1,kgrid
         sum=sum+Apro(k)*vprime(k,1)*dr(k,1)
      enddo
c
      volintk=sum
      return
      end
