      subroutine inter_cspl(n,r,datad,m,x,ds)
c
c... INTERPOLATE DATA TO NEW GRID
c
      implicit none
c
      integer j, k, n, m, ifail
      real*8 r(n), datad(1,n), x(m), ds(m), br(1,n)
      logical sk
c
      parameter (sk = .TRUE.)
c 
      do j=1,n
       br(1,j)=0.D0
      enddo
c
c      write(*,*) 'inside inter_cspl ...'
c      write(*,*) 'n,m = ',n,m
      ifail=0
c     call e01bee(n,r,data,br,ifail)
c     call e01bfe(n,r,data,br,m,x,ds,ifail)
c
c      do j=1,n
c       write(*,100)j,r(j),datad(1,j),br(1,j)
c      enddo
c
c      do j=1,m
c       write(*,100) j, x(j), ds(j)
c      enddo
c
      call dpchim(n,r,datad,br,1,ifail)
c      write(*,*) 'after call to dpchim ...',ifail,sk
c      do j=1,n
c       write(*,100) j,r(j),datad(1,j),br(1,j)
c      enddo
      call dpchfe(n,r,datad,br,1,sk,m,x,ds,ifail)
c      write(*,*) 'after call to dpchfe ...',ifail
c
 100  format(i2,2x,1p3e14.6)
      return
      end
