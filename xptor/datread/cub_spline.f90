! cub_spline.f90
!
! PURPOSE:
!  Map a raw dataset [xm(1:n),ym(1:n)] to an 
!  arbitrarily fine dataset [xi(1:ni),yi(1:ni)]
!  using cubic spline interpolation.
!
!   INPUT: xm(1:n),ym(1:n),n,xi(1:ni),ni
!
!  OUTPUT: yi(1:ni)
!
!  Note that the x-coordinates, xi(:), of the 
!  interpolates can be unequally-spaced.  
!
! REVISIONS
! 19 Aug 04: jc
!  Created.
! 25 May 05: jc
!  Moved to gyro/math subdir.
!-------------------------------------------------------

      subroutine cub_spline(XM,YM,n,xi,yi,ni)

      implicit none

      integer :: n,i,m,j,k,i_count
      integer :: ni,ii
  
      real, dimension(n) :: XM,YM,LM,UM,DM,EM
      real, dimension(0:n) :: CM
      real, dimension(ni) :: xi,yi
      real :: xv,fact,e,ff,g,h

      m = n - 1
      j = m - 1
      k = j - 1

  !  GENERATION OF TRIDIAGONAL SYSTEM FOR SECOND DERIVATIVE

  ! Added 06 Oct 04: jc
      CM(:) = 0.0

      do i = 1,j
       DM(i) = 2.* (XM(i+2) - XM(i))
       CM(i) = 6.*(YM(i+2)-YM(i+1)) / (XM(i+2)-XM(i+1)) &
         + 6.*(YM(i) - YM(i+1)) / (XM(i+1) - XM(i))
      enddo

      do i = 2,j
       LM(i) = XM(i+1) - XM(i)
      enddo
      do i = 1,k
       UM(i) = XM(i+2) - XM(i+1)
      enddo

!  SOLUTION OF TRIDIAGONAL SYSTEM

      call TRIDI(LM,DM,UM,CM,j)

!  EVALUATION AND PRINTING OF CUBIC SPLINES

      CM(0) = 0.0
      CM(n) = 0.0

      do ii=1,ni

       xv = xi(ii)

       do i = 1,m
        fact = XM(i+1) - XM(i)
        e = CM(i-1) / (6.0*fact)
        ff = CM(i) / (6.0*fact)
        g = (YM(i)/fact) - (CM(i-1)*fact/6.)
        h = (YM(i+1)/fact) - (CM(i)*fact/6.)
        EM(i) = e*(XM(i+1) - xv)**3 + ff*(xv - XM(i))**3 + &
             g*(XM(i+1)- xv) + h*(xv - XM(i))
       enddo

      i_count = 1

60    if (xv < XM(i_count+1)) go to 70

      i_count = i_count + 1

      go to 60

70    continue 

      yi(ii) = EM(i_count)

      enddo

      return

      end subroutine cub_spline

      subroutine TRIDI(L,D,U,COE,n)

      integer :: n,m,i
      real :: L(n),D(n),U(n),COE(n)

      m = n - 1 
      do i = 1,m
       L(i+1) = L(i+1) / D(i)
       D(i+1) = D(i+1) - L(i+1) *U(i)
       COE(i+1) = COE(i+1) - L(i+1) *COE(i)
      enddo

      COE(n) = COE(n)/D(n)
      do i = m,1,-1
       COE(i) = (COE(i) - U(i) *COE(i+1)) / D(i)
      enddo
 
      end subroutine TRIDI

