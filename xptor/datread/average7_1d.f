      subroutine average7_1d(f,nk)
c*****************************************************
c
c performs a two pass seven point average 
c
c*****************************************************
      implicit none
c
      integer nk,k,i,m
      real*8 f(nk),g(nk)
c
c check if grid is too small
      if(nk.lt.7)return
c
        g(1)=f(1)
        g(nk)=f(nk)
        k=2
        g(k)=(f(k)+f(k+1)+f(k-1))/3.D0
        k=3
        g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.D0
        k=nk-1
        g(k)=(f(k)+f(k+1)+f(k-1))/3.D0
        k=nk-2
        g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.D0
        do k=4,nk-3
          g(k)=(f(k)+f(k-1)+f(k+1)+f(k-2)+f(k+2)
     >          + f(k+3)+f(k-3))/7.D0
        enddo
        k=2
        f(k)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=3
        f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        k=nk-1
        f(k)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=nk-2
        f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        do k=4,nk-3
          f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2)
     >            + g(k+3)+g(k-3))/7.D0
        enddo
      return
      end
   
