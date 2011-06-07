      SUBROUTINE W_BICUBIC(x,nx,y,ny,c,nc,xl,yl,npds,pds,iflag,message)
!***********************************************************************
!W_BICUBIC evaluates bicubic splines
!References:
!  W.A.Houlberg 3/2000
!Input:
!  x(i)-array of x ordinates
!  nx-number of x elements
!  y(j)-array of y ordinates
!  ny-number of y elements
!  c(2,nc,*)-array of bicubic spline coefficients
!  nc-second dimension of c
!  xl-target x value
!  yl-target y value
!  npds-number of values and partial derivatives
!Output:
!  pds(k)-array of values and partial derivatives
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        nc,                      npds,
     &               nx,                      ny
      REAL           c(2,nc,*),
     &               x(*),                    y(*),
     &               xl,                      yl
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           pds(6)
!Declaration of local variables
      INTEGER        i,                       j,
     &               k,                       km1,
     &               kp1,                     kp2,
     &               lxl,                     lx,
     &               ly,                      l,
     &               lx1
      REAL           d,                       dy,                       
     &               h,                       hx,                      
     &               s0,                      sh,
     &               sp0,                     sph,
     &               u,                       v
      REAL           su(2),                   sv(2),
     &               sux(2),                  suy(2),
     &               svx(2),                  sxy(2)
      REAL           spln0,                   spln1,
     &               spln2
      spln0(s0,sh,sp0,sph,h,d)=s0+d*(h*sp0+d*(3.0*(sh-s0)
     &                         -(sph+2.0*sp0)*h+d*(2.0*(s0-sh)
     &                         +(sph+sp0)*h)))
      spln1(s0,sh,sp0,sph,h,d)=sp0+d*(6.0*(sh-s0)/h-2.0*(sph+2.0*sp0)
     &                         +3.0*d*(2.0*(s0-sh)/h+(sph+sp0)))
      spln2(s0,sh,sp0,sph,h,d)=6.0*(sh-s0)/h**2-2.0*(sph+2.0*sp0)/h
     &                         +d*(2.0*(s0-sh)/h**2+(sph+sp0)/h)*6.0
!Initialization
      lx=0
      ly=0
      iflag=0
!Correlated table search for xl
      CALL RARRAY_SEARCH(x,nx,xl,lx)
      IF(lx.eq.0) THEN
        iflag=1
        message='W_BICUBIC(1)/ERROR:x<x(1)'
        GOTO 1000
      ELSEIF(lx.eq.nx) THEN
        iflag=1
        message='W_BICUBIC(2)/ERROR:x>x(n)'
        GOTO 1000
      ENDIF
!Correlated table search for yl
      CALL RARRAY_SEARCH(y,ny,yl,ly)
      IF(ly.eq.0) THEN
        iflag=1
        message='W_BICUBIC(1)/ERROR:y<y(1)'
        GOTO 1000
      ELSEIF(ly.eq.ny) THEN
        iflag=1
        message='W_BICUBIC(1)/ERROR:y>y(n)'
        GOTO 1000
      ENDIF
      lx1=lx+1
      hx=x(lx1)-x(lx)
      dy=y(ly+1)-y(ly)
      u=(xl-x(lx))/hx
      v=(yl-y(ly))/dy
      k=2*ly
      kp1=k+1
      kp2=k+2
      km1=k-1
      DO l=1,2
        lxl=lx-1+l
        i=2*(ly-1+l)
        j=i-1
        sv(l)=spln0(c(1,lxl,km1),c(1,lxl,kp1),c(1,lxl,k),c(1,lxl,kp2),
     &              dy,v)
        svx(l)=spln0(c(2,lxl,km1),c(2,lxl,kp1),c(2,lxl,k),c(2,lxl,kp2),
     &               dy,v)
        IF(npds.ge.3) THEN
          su(l)=spln0(c(1,lx,j),c(1,lx1,j),c(2,lx,j),c(2,lx1,j),hx,u)
          suy(l)=spln0(c(1,lx,i),c(1,lx1,i),c(2,lx,i),c(2,lx1,i),hx,u)
          IF(npds.ge.4) THEN
            sux(l)=spln1(c(1,lx,j),c(1,lx1,j),c(2,lx,j),c(2,lx1,j),hx,u)
            sxy(l)=spln1(c(1,lx,i),c(1,lx1,i),c(2,lx,i),c(2,lx1,i),hx,u)
          ENDIF
        ENDIF
      ENDDO
      pds(1)=spln0(sv(1),sv(2),svx(1),svx(2),hx,u)
      IF(npds.GT.1) THEN
        pds(2)=spln1(sv(1),sv(2),svx(1),svx(2),hx,u)
        IF(npds.GT.2) THEN
          pds(3)=spln1(su(1),su(2),suy(1),suy(2),dy,v)
          IF(npds.GT.3) THEN
            pds(4)=spln1(sux(1),sux(2),sxy(1),sxy(2),dy,v)
            IF(npds.GT.4) THEN
              pds(5)=spln2(sv(1),sv(2),svx(1),svx(2),hx,u)
              IF(npds.GT.5) THEN
                pds(6)=spln2(su(1),su(2),suy(1),suy(2),dy,v)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
 1000 RETURN
      END
