  MODULE smthspline
! smoothing cubic spline module.
    
    USE   nrtype,                               ONLY :DP,I4B

    IMPLICIT NONE
    
    CONTAINS


      SUBROUTINE smooth_values(ydata,xdata,ndata,smthdata,dsmthdatadx,d2smthdatadx2,  &
                               weight,no_weight,sdevfctr,ibc,derivl,derivr)
!------------------------------------------------------------------------------------------
!   INPUT
!   ydata(i)          noisy data , i =1,ndata
!   xdata(i)          abscissa associated with ydata
!   weights(i)        the weighting of each data point. 
!                     small weights lead to an interpolatory spline
!   no_weight = 0     means weights are given,eg. weight is input
!   no_weight = 1     means weights are internally 
!                     calculated in this routine (using the variance of the data)
!                     weight(i) is then defined in this routine  and becomes an
!                     output
!   sdevfctr          A multiplier for the weights,used if no_weight = 1.
!                     hence the actual weight becomes some (possibly fractional)
!                     multiple of the standard deviation
!   ibc               boundary condition specifier:
!                      ibc =1 input first deriv at ends,derivl at xdata(1)
!                                                       derivr at x(ndata)
!                          =2 means second deriv = 0 at ends, (no input required)
!                          =3 means periodic at ends,         (no input required)

!   OUPUT
!   smthdata(i)        the smoothed data
!   dsmthdatadx(I)     first deriv
!   d2smthdatadx2(I)   second derivative 
!   weight(i)          see above
! -------------------------------------------------------------------HSJ----------
       IMPLICIT NONE

       REAL(DP),DIMENSION(:) :: ydata,xdata,smthdata,weight,dsmthdatadx(ndata),           &
                                d2smthdatadx2(ndata)    
       !temporary arrays:
       REAL(DP),DIMENSION(:) :: a(ndata),b(ndata),c(ndata),d(ndata),e(ndata),g(ndata),    &  
                                zm(ndata),aa(ndata),dd(ndata),                            &
                                u(ndata+2),v(ndata+2),w(ndata+2),p(ndata+2),q(ndata+2),   &
                                r(ndata+2),s(ndata+2),t(ndata+2) 


       REAL(DP)              :: derivl,derivr,ave,sdev,adev,var,skew,curt,                &
                                deriv1,deriv2,sdevfctr
       INTEGER(I4B)          :: ndata, no_weight,ibc,ind,i


       a(:)=0.0_DP; b(:)=0.0_DP; c(:)=0.0_DP ;d(:)=0.0_DP ;e (:)=0.0_DP
       g(:)=0.0_DP ; zm(:) =0.0_DP ; aa(:) =0.0_DP ;  dd(:) =0.0_DP  
       u(:) =0.0_DP ;  v(:) =0.0_DP ;  w(:) =0.0_DP 
       p(:) =0.0_DP ;  q(:) =0.0_DP ;  r(:) =0.0_DP 
       s(:) =0.0_DP ;  t(:) =0.0_DP
       
       IF(no_weight == 1)THEN
          CALL moment(ydata,ndata,ave,adev,sdev,var,skew,curt)
          weight(:) = sdevfctr*sdev
       ENDIF


       ind = 0 ! determine spline coefficients 
       CALL Smspline(ndata,ibc,ind,xdata,ydata,a,b,c,d,e,g,weight,        &
                      derivl,derivl,u,v,w,p,q,r,s,t,aa,dd,zm,xdata(1),    &
                      smthdata(1),deriv1,deriv2)







       ind = 1  ! evaluate spline and first and seconde derivative point xx
       DO i = 1,ndata
            CALL Smspline(ndata,ibc,ind,xdata,ydata,a,b,c,d,e,g,weight,    &
                          derivr,derivl,u,v,w,p,q,r,s,t,aa,dd,zm,xdata(i),  &
                          smthdata(i),deriv1,deriv2)
       ENDDO




       RETURN

       END SUBROUTINE smooth_values







        SUBROUTINE Progon5(n,a,b,c,d,e,g,u,v,w, &
                            p,q,r,s,t,aa,dd,z)

!
!***********************************************************
!    Program Progon5 for solving linear system of the      *
!    form                                                  *
! Ú                                          ¿ Ú  ¿ Ú  ¿   *
! ³a1   b1   c1    0   ..........0   e1   d1 ³ ³z1³ ³g1³   *
! ³d2   a2   b2   c2    0 ........    0   e2 ³ ³z2³ ³g2³   *
! ³e3   d3   a3   b3   c3    0 ...    0    0 ³ ³z3³ ³g3³   *
! ³O    e4   d4   a4   b4   c4   .         0 ³ ³z4³=³g4³   *
! ³..........................................³ ³..³ ³  ³   *
! ³                                          ³ ³  ³ ³  ³   *
! ³bn   cn   0     0  .....  0   en  dn   an ³ ³zn³ ³gn³   *
! À                                          Ù À  Ù À  Ù   *
!  with pentadiagonal matrix by sweep method               *
!  n - number of equations                                 *
!  a,b,c,d,e - arrays (of size n) for elements of          *
!              the matrix's diagonal                       *
!  g - array (of size n) for elements of the right parts   *
!                                                          *
!  u,v,w,p,q,                                              *
!  r,s,t,aa,dd - work arrays (of size n+2)                 *
!                                                          *
!        Output:                                           *
!  z - array (of size n) containing solution of the        *
!      system                                              *
!***********************************************************
!
       REAL(DP)  a(1),b(1),c(1),d(1),e(1),g(1),u(1),aa(1),r(1)
       REAL(DP)  v(1),w(1),p(1),q(1),t(1),s(1),z(1),dd(1)
       REAL(DP)  a11,a12,a21,a22,b1,b2
       INTEGER(I4B) n,i
!
              v(1) = 0.
              v(2) = 0.
              w(1) = 0.
              w(2) = 0.
              u(1) = 0.
              u(2) = 0.
!
        DO i = 1,n
         dd(i)  =  d(i) + e(i)*v(i)
         aa(i)  =  a(i) + dd(i)*v(i+1) +e(i)*w(i)
         u(i+2) = (g(i) - dd(i)*u(i+1) - e(i)*u(i))/aa(i)
         v(i+2) = -(b(i)+dd(i)*w(i+1))/aa(i)
         w(i+2) = -c(i)/aa(i)
      ENDDO
!
            p(1) = 0.
            q(2) = 0.
            p(2) = 1.
            q(1) = 1.
!
           DO i = 1,n
         p(i+2) = -(dd(i)*p(i+1) + e(i)*p(i))/aa(i)
         q(i+2) = -(dd(i)*q(i+1) + e(i)*q(i))/aa(i)
      ENDDO
!
                t(n-1) = 0.
                s(n-1) = 0.
                t(n)   = 0.
                r(n)   = 0.
             r(n-1) = 1.
             s(n)   = 1.
!
       DO i = n-2,1,-1
          t(i) = v(i+2)*t(i+1) +w(i+2)*t(i+2) + u(i+2)
          s(i) = v(i+2)*s(i+1) +w(i+2)*s(i+2) + p(i+2)
          r(i) = v(i+2)*r(i+1) +w(i+2)*r(i+2) + q(i+2)
       ENDDO
!
       a11 =  1. - q(n+1) - w(n+1)*r(1)
       a12 = -(p(n+1) + v(n+1) + w(n+1)*s(1))
       a21 = -(v(n+2)*r(1) + w(n+2)*r(2) + q(n+2))
       a22 =  1. - p(n+2) - v(n+2)*s(1) - w(n+2)*s(2)
       b1  = w(n+1)*t(1) + u(n+1)
       b2  = v(n+2)*t(1) + w(n+2)*t(2) + u(n+2)
       z(n-1) = ( b1*a22 - b2*a12)/(a11*a22-a12*a21)
       z(n)   = (-b1*a21 + b2*a11)/(a11*a22-a12*a21)
!
       DO i = 1,n-2
          z(i) = t(i) + s(i)*z(n) + r(i)*z(n-1)
       ENDDO
!
       RETURN
       END SUBROUTINE progon5




       SUBROUTINE Smspline(n,ib,ind,x,y,a,b,c,d,e,g,rho, &
          ax,bx,u,v,w,p,q,r,s,t,aa,dd,zm,xx,sp,dsp,d2sp)
!
!************************************************************
!       Program Smspline constructs the smoothing spline    *
!                for given array of points                  *
!                                                           *
!  n   - number of points in array                          *
!  x   - array (of size n) containing x-coordinates         *
!  y   - array (of size n) containing y-coordinates         *
!  rho - array (of size n) containing weight coefficients   *
!                                                           *
!  ind - mode type:                                         *
!        ind = 0 - program calculates parameters of spline  *
!        ind = 1 - the spline parameters are known          *
!  ib  - type of boundary conditions                        *
!         ib = i (i=1-3) - conditions of ith type           *
!                                                           *
!  a,b,c,d,e,g,zm,aa,dd - work arrays (of size n)           *
!  u,v,w,p                                                  *
!  q,r,s,t - work arrays (of size n+2)                      *
!  ax, bx - parameter values in boundary conditions         *
!           of the first type                               *
!                                                           *
!         Output:                                           *
!  sp   - value of spline at point xx                       *
!  dsp  - value of its first derivative at point xx         *
!  d2sp - value of its second derivative at point xx        *
!                                                           *
!************************************************************
!
       REAL(DP)  a(1),b(1),c(1),d(1),e(1), &
            g(1),rho(1),aa(1),r(1)
       REAL(DP) u(1),v(1),w(1),p(1),q(1),t(1),s(1),zm(1),dd(1)
       REAL(DP) x(1), y(1),ax,bx,xx,sp,dsp,d2sp,h1,h2,h3,h4,di,h,tt
       INTEGER(I4B) ib,ind,n,nn,i,ni
!
        nn = n
        IF (ind.EQ.0) THEN
        SELECT CASE (ib)
!
        CASE (1)
         h1     = x(2)-x(1)
         a(1)   = h1/3.+(rho(1)+rho(2))/h1**2
         h2     = x(3) - x(2)
         b(1)   = h1/6.-(1./h1 +1./h2)*rho(2)/h1 &
                  - rho(1)/h1**2
         c(1)   = rho(2)/(h1*h2)
         g(1)   = (y(2) - y(1))/h1 - ax
         e(3)   = c(1)
         h1     = x(n) - x(n-1)
         h2     = x(n-1) - x(n-2)
         a(n)   = h1/3. +(rho(n-1)+rho(n))/h1**2
         b(n-1) = h1/6.-(1./h1 +1./h2)*rho(n-1)/h1 &
                 -rho(n)/h1**2
         c(n-2) = rho(n-1)/(h1*h2)
         e(1) = 0.
         e(2) = 0.
         d(1) = 0.
         d(2) = b(1)
         c(n-1) = 0.
         c(n) = 0.
         b(n) = 0.
         d(n) = b(n-1)
         g(n) = bx - (y(n) - y(n-1))/h1
         e(n) = c(n-2)
!
        CASE (2)
         a(1)   = 1.
         a(n)   = 1.
         b(1)   = 0.
         b(n-1) = 0.
         c(1)   = 0.
         c(n-2) = 0.
         g(1)   = 0.
         g(n)   = 0.
         e(1)   = 0.
         e(2)   = 0.
         e(3)   = 0.
         d(1)   = 0.
         b(n)   = 0.
         c(n)   = 0.
         c(n-1) = 0.
!
        CASE (3)
         h1 = x(2) - x(1)
         h2 = x(n) - x(n-1)
         h3 = x(3) - x(2)
         a(1) = (h1 + h2)/3. &
               + rho(n-1)/h2**2 + rho(2)/h1**2 + &
                (1./h1 + 1./h2)**2*rho(1)
         b(1) = h1/6. -((1./h2 + 1./h1)*rho(1)+ &
                        (1./h2 + 1./h3)*rho(2))/h2
         c(1) = rho(1)/(h1*h3)
         g(1) = (y(2) - y(1))/h1 - (y(1) - y(n-1))/h2
         d(2) = b(1)
         e(3) = c(1)
         h4   = x(n-1) - x(n-2)
         c(n-2) = rho(n-1)/(h2*h4)
         c(n-1) = rho(n)/(h2*h1)
         b(n-1) = h2/3. - ((1./h4 + 1./h2)*rho(n-1) + &
                           (1./h2 + 1./h1)*rho(n))/h2
         d(1) = b(n-1)
         e(1) = c(n-2)
         e(2) = c(n-1)
         g(n-1) = (y(n) - y(n-1))/h2 - (y(n-1) - y(n-2))/h4
         nn = n - 1
       END SELECT
!
      DO i = 2, n-1
         h1   = x(i)   - x(i-1)
       h2   = x(i+1) - x(i)
         a(i) = (h1+h2)/3. &
               + rho(i-1)/h1**2 + rho(i+1)/h2**2 + &
                 (1./h1 + 1./h2)**2*rho(i)
       g(i) = (y(i+1)-y(i))/h2 - (y(i)-y(i-1))/h1
            IF (i.LE.(n-2)) THEN
                h3 = (x(i+2)-x(i+1))
              b(i) = h2/6. - ((1./h1 + 1./h2)*rho(i) &
                               +(1./h2 + 1./h3)*rho(i+1))/h2
              d(i+1) = b(i)
            ENDIF
         IF (i.LE.(n-3)) THEN
               c(i) = rho(i+1)/(h2*h3)
               e(i+2) = c(i)
           ENDIF
      ENDDO
!
        CALL Progon5 &
             (nn,a,b,c,d,e,g,u,v,w,p,q,r,s,t,aa,dd,zm)
        aa(1) = y(1) - rho(1)*(zm(2)-zm(1))/(x(2)-x(1))
      aa(n) = y(n) + rho(n)*(zm(n)-zm(n-1))/(x(n)-x(n-1))

        IF (ib.EQ.3) THEN
            zm(n) = zm(1)
            di = (zm(2)-zm(1))/(x(2)-x(1))
            di = di - (zm(n)-zm(n-1))/(x(n)-x(n-1))
            aa(1) = y(1) - rho(1)*di
            aa(n) = aa(1)
        ENDIF
!
      DO i = 2, n-1
            di = (zm(i+1)-zm(i))/(x(i+1)-x(i))
            di = di - (zm(i)-zm(i-1))/(x(i)-x(i-1))
            aa(i) = y(i) - rho(i)*di
      ENDDO
!
      ELSE
!
        DO  i = 2,n
            ni= i
            IF(x(i).GT.xx) EXIT
        ENDDO
!
        i  = ni - 1
        h  = x(i+1) - x(i)
        tt = (xx - x(i))/h
        sp = aa(i)*(1.-tt) + aa(i+1)*tt &
            - h**2/6.*tt*(1.-tt)* &
             ((2.-tt)*zm(i) + (1.+tt)*zm(i+1))
        dsp= (aa(i+1) - aa(i))/h - &
              h/6.*((2. - 6.*tt + 3.*tt**2)*zm(i) &
            +(1. - 3.*tt**2)*zm(i+1))
        d2sp = (1. - tt)*zm(i) + tt*zm(i+1)
!
      ENDIF
      RETURN

      END       SUBROUTINE Smspline



        
      SUBROUTINE moment(DATA,n,ave,adev,sdev,var,skew,curt)
! ORIGINALlY FROM NUMERICAL RECEPIES
      IMPLICIT NONE
      INTEGER(I4B) n,j
      REAL(DP)  adev,ave,curt,sdev,skew,var,DATA(n)
      REAL(DP)  p,s,ep

      IF(n .LE. 1) THEN
         PRINT *,'Moments not calculated'
         AVE = DATA(1)
         ADEV =DATA(1)
         SDEV =DATA(1)
         VAR =SDEV*SDEV
         SKEW =0
         CURT =0
         RETURN
      ENDIF
      s=0.
      DO 11 j=1,n
        s=s+DATA(j)
11    CONTINUE
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      DO 12 j=1,n
        s=DATA(j)-ave
        ep=ep+s
        adev=adev+ABS(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    CONTINUE
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=SQRT(var)

      IF(var.NE.0.)THEN
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      ELSE
        PRINT *, 'no skew or kurtosis when zero variance in moment'
        SKEW =0.0
        CURT =0.0
      ENDIF

      RETURN

      END       SUBROUTINE moment


  END MODULE smthspline
