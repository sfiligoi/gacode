      SUBROUTINE V_BICUBIC(f,x,nx,y,ny,c,mx_nx,wk,iflag,message)
!***********************************************************************
!V_BICUBIC calculates bicubic spline coefficients
!References:
!  W.A.Houlberg 3/2000
!Input:
!  f(i,j)-array of values at (x(i),y(j))
!  x(i)-array of x values
!  nx-number of x values .ge. 4
!  y(i)-array of y values
!  ny-number of y values .ge. 4
!  mx_nx-first dimension of f and x arrays
!  wk(*)-work vector of length 2*nx*ny+2*max(nx,ny)
!Output:
!  c-array of spline coefficients
!   -c(1,i,j)=s
!   -c(2,i,j)=ds/dx
!   -c(1,i,j+ny)=ds/dy
!   -c(2,i,j+ny)=d(ds/dx)/dy
!    where s(x,y) is the spline approximation
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        nx,                      ny,
     &               mx_nx
      REAL           x(*),                    y(*),
     &               f(mx_nx,*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           c(2,mx_nx,*),            wk(*)
!Declaration of local variables
      INTEGER        iwk
      IF(nx.gt.mx_nx) THEN
        iflag=1
        message='V_BICUBIC/ERROR:nx greater than dimension'
        GOTO 1000
      ENDIF
      IF(nx.lt.4) THEN
        iflag=1
        message='V_BICUBIC/ERROR:nx less than 4'
        GOTO 1000
      ENDIF
      IF(ny.lt.4) THEN
        iflag=1
        message='V_BICUBIC/ERROR:ny less than 4'
        GOTO 1000
      ENDIF
      iwk=2*nx*ny
      CALL V_BICUBIC_NUC(x,f,nx,ny,wk(iwk+1),wk,mx_nx,ny,iflag,message)
      IF(iflag.gt.0) THEN
        iwk=LEN(message)
        message='V_BICUBIC(1)/'//message(1:iwk-13)
        GOTO 1000
      ENDIF
      CALL V_BICUBIC_NUC(y,wk,ny,2*nx,wk(iwk+1),c,ny,2*mx_nx,iflag,
     &                   message)
      IF(iflag.gt.0) THEN
        iwk=LEN(message)
        message='V_BICUBIC(2)/'//message(1:iwk-13)
        GOTO 1000
      ENDIF
 1000 RETURN
      END
      SUBROUTINE V_BICUBIC_NUC(tau,gtau,n,m,w,vs,ic1,ic2,iflag,message)
!***********************************************************************
!V_BICUBIC_NUC provides the nucleus of the bicubic spline routine
!  VBICUBIC
!References:
!  W.A.Houlberg 3/2000
!Input:
!tau-
!gtau-
!n-
!m-
!ic1-
!ic2-
!Output:
!w-()
!vs()-
!Output:
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        n,                       m,
     &               ic1,                     ic2
      REAL           tau(n),                  gtau(ic1,*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           w(n,2),                  vs(ic2,2,*)
!Declaration of local variables
      INTEGER        i,                       j,
     &               jj,                      jm1,
     &               jp1,                     k
      REAL           aa,                      bb,
     &               c1,                      c2,
     &               cc,                      dd,
     &               dtau,                    g,
     &               h,                       ratio,
     &               u,                       xilim
      w(2,1)=tau(3)-tau(1)
      IF(w(2,1).le.0.0) THEN
        iflag=1
        message='V_BICUBIC_NUC(1)/ERROR:grid not monotonic'
        GOTO 1000
      ENDIF
      DO k=1,m
        vs(k,1,1)=gtau(1,k)
      ENDDO
      xilim=tau(1)
      IF(n.ge.5) THEN
        xilim=tau(n-2)
        DO i=2,n-3
          j=i+1
          w(j,1)=tau(i+2)-tau(j)
          IF(w(j,1).le.0.0) THEN
            iflag=1
            message='V_BICUBIC_NUC(2)/ERROR:grid not monotonic'
            GOTO 1000
          ENDIF
          DO k=1,m
            vs(k,1,i)=gtau(j,k)
          ENDDO
        ENDDO
      ENDIF
      w(n-2,1)=tau(n)-xilim
      IF(w(n-2,1).le.0.0) THEN
        iflag=1
        message='V_BICUBIC_NUC(3)/ERROR:grid not monotonic'
        GOTO 1000
      ENDIF
      DO k=1,m
        vs(k,1,n-2)=gtau(n,k)
      ENDDO
      DO i=2,n-2
        DO k=1,m
          vs(k,2,i)=(vs(k,1,i)-vs(k,1,i-1))/w(i,1)
        ENDDO
      ENDDO
      dtau=tau(2)-tau(1)
      ratio=dtau/w(2,1)
      w(1,2)=(ratio-1.0)**2
      w(1,1)=ratio*(ratio-1.0)
      c1=ratio*(2.0*ratio-3.0)
      DO k=1,m
        vs(k,2,1)=(gtau(2,k)-gtau(1,k))/dtau+vs(k,2,2)*c1
      ENDDO
      IF(n.ge.5) THEN
        DO i=2,n-3
          j=i+1
          jj=i-1
          g=-w(j,1)/w(jj,2)
          c1=3.0*w(i,1)
          c2=3.0*w(j,1)
          DO k=1,m
            vs(k,2,i)=g*vs(k,2,jj)+c1*vs(k,2,j)+c2*vs(k,2,i)
          ENDDO
          w(i,2)=g*w(jj,1)+2.0*(w(i,1)+w(j,1))
        ENDDO
      ENDIF
      dtau=tau(n-1)-xilim
      ratio=dtau/w(n-2,1)
      g=-(ratio-1.0)**2/w(n-3,2)
      w(n-2,2)=ratio*(ratio-1.0)
      c1=ratio*(2.0*ratio-3.0)
      DO k=1,m
        vs(k,2,n-2)=(gtau(n-1,k)-vs(k,1,n-3))/dtau+vs(k,2,n-2)*c1
      ENDDO
      w(n-2,2)=g*w(n-3,1)+w(n-2,2)
      DO k=1,m
        vs(k,2,n-2)=(g*vs(k,2,n-3)+vs(k,2,n-2))/w(n-2,2)
      ENDDO
      DO j=n-3,1,-1
        DO k=1,m
          vs(k,2,j)=(vs(k,2,j)-w(j,1)*vs(k,2,j+1))/w(j,2)
        ENDDO
      ENDDO
      DO k=1,m
        DO jj=1,n
          j=n+1-jj
          IF(j.eq.1) THEN
            jm1=j
          ELSEIF(j.eq.n) THEN
            jm1=j-2
          ELSE
            jm1=j-1
          ENDIF
          DO i=1,2
            vs(k,i,j)=vs(k,i,jm1)
          ENDDO
        ENDDO
        DO j=2,n-1,n-3
          jm1=j-1
          jp1=j+1
          IF(jm1.eq.2) jm1=1
          IF(jp1.eq.n-1) jp1=n
          h=tau(jp1)-tau(jm1)
          u=tau(j)-tau(jm1)
          aa=vs(k,1,jm1)
          bb=vs(k,2,jm1)
          cc=(3.0*(vs(k,1,jp1)-vs(k,1,jm1))/h-(vs(k,2,jp1)
     &                                         +2.0*vs(k,2,jm1)))/h
          dd=(2.0*(vs(k,1,jm1)-vs(k,1,jp1))/h+(vs(k,2,jp1)
     &                                         +vs(k,2,jm1)))/h**2
          vs(k,1,j)=aa+u*(bb+u*(cc+dd*u))
          vs(k,2,j)=bb+u*(2.0*cc+3.0*dd*u)
        ENDDO
      ENDDO
 1000 RETURN
      END
