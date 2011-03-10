cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine  lspolysmooth(korder,m,kmax,kx,fx)
c     Adapted by J. Kinsey from original lspolysmooth which 
c     used NAG libraries for implementation on the Cray
c
       implicit none
       integer ifail,ifail1,j,kplus1,nplus1,ndeg
       integer kmax,m,korder
       real*8 fx(0:kmax), kx(0:kmax), a(200)
       real*8 x(m),y(m),w(m),work1(3*m),work2(2*(korder+1))
       real*8 s(korder+1)
       real*8 xcap,p,ak,akmin,akmax,eps
c
       eps = 0.D0
       p = 0.D0
       akmin=kx(1)
       akmax=kx(m)
c
       do j=1,200
         a(j) = 0.D0
       enddo
c
       do j=1,m
        ak=kx(j)
        y(j)=fx(j)
        x(j)=ak
        w(j)=1.D0
        if(j.eq.1) w(j)=.1D0
        if(j.eq.m) w(j)=.1D0
        write(*,*) j,x(j),y(j)
       enddo
c
       kplus1=korder+1
       call dpolft(m,x,y,w,kplus1,ndeg,eps,s,ifail,a)
       write(*,*) 'call to dpolft okay ...'
       write(*,*) 'ndeg = ',ndeg
       do j=1,m
         write(*,*) j, a(j)
       enddo
c
       do j=1,m
        xcap=((x(j)-x(1))-(x(m)-x(j)))/(x(m)-x(1))
        write(*,*) j, xcap, a(j)
        call dp1vlu(ndeg,0,xcap,p,a)
       write(*,*) 'call to dp1vlu okay ...'
        fx(j)=p
        write(*,*) j,xcap,fx(j)
       enddo
c
      return
      end  
