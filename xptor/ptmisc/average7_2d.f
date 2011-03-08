      subroutine average7_2d(f,na,nk,nf)
c*****************************************************
c
c performs a two pass seven point average 
c
c*****************************************************
      implicit none
c
      integer nf,nk,na,k,i,m
      real*8 f(na,nf),g(na)
c
      do k=1,nk
      enddo
c
c check if grid is too small
      if(nk.lt.7)return
c
      do i=1,nf
        g(1)=f(1,i)
        g(nk)=f(nk,i)
        k=2
        g(k)=(f(k,i)+f(k+1,i)+f(k-1,i))/3.D0
        k=3
        g(k)=(f(k,i)+f(k+1,i)+f(k-1,i)+f(k+2,i)+f(k-2,i))/5.D0
        k=nk-1
        g(k)=(f(k,i)+f(k+1,i)+f(k-1,i))/3.D0
        k=nk-2
        g(k)=(f(k,i)+f(k+1,i)+f(k-1,i)+f(k+2,i)+f(k-2,i))/5.D0
        do k=4,nk-3
          g(k)=(f(k,i)+f(k-1,i)+f(k+1,i)+f(k-2,i)+f(k+2,i)
     >          + f(k+3,i)+f(k-3,i))/7.D0
        enddo
        k=2
        f(k,i)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=3
        f(k,i)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        k=nk-1
        f(k,i)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=nk-2
        f(k,i)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        do k=4,nk-3
          f(k,i)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2)
     >            + g(k+3)+g(k-3))/7.D0
        enddo
      enddo
c
      return
      end
   
