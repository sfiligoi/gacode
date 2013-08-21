 MODULE spline
   USE nrtype,                             ONLY : DP,I4B        
   IMPLICIT NONE
   REAL(DP),ALLOCATABLE, DIMENSION(:) :: s_coef
   REAL(DP) bc_0,bc_1,rho_w,sp,dsp,d2sp
   INTEGER(I4B) s_bc
   LOGICAL get_scoef


   CONTAINS
   SUBROUTINE Progon3 (a,b,c,d,x,n)
     !***********************************************************
     !        Program Progon3 for solving linear system         *
     !        of the form                                       *
     !                                                          *
     !        a(1)*x(1) + b(1)*x(2)      + c(1)*x(n) = d(1)     *
     !        ............................................      *
     !        c(i)*x(i-1) + a(i)*x(i)  + b(i)*x(i+1) = d(i)     *
     !        .............................................     *
     !        b(n)*x(1)     + a(n)*x(n-1) + b(n)*x(n)= d(n)     *
     !        with tridiagonal matrix by sweep method           *
     !    n - number of equations                               *
     !  a,b,                                                    *
     !  c,d - arrays (of size n) containing elements of the     *
     !        matrix' diagonals and of the right parts          *
     !  w,s,t,                                                  *
     !  u,v - work arrays (of size n+1)                         *
     !        Output:                                           *
     !    x - array (of size n) containing solution of the      *
     !        system                                            *
     !***********************************************************
     USE nrtype,                         ONLY : DP,I4B
     IMPLICIT NONE
     integer(I4B) i,n,i1
     REAL(DP)  a(n),b(n),c(n),d(n),x(n)
     !        tEMPORARY ARRAYS:
     REAL(DP) u(n+1),v(n+1),w(n+1),s(n+1),t(n+1)
     REAL(DP) z

     u(1) = 0.
     v(1) = 0.
     w(1) = 1.
     !
     DO  i = 1,n
        i1 = i + 1
        z  = 1._DP/(a(i) + c(i)*v(i))
        v(i1) = -b(i)*z
        u(i1) = (-c(i)*u(i) + d(i))*z
        w(i1) = - c(i)*w(i)*z
     ENDDO
     !
     s(n) = 1._DP
     t(n) = 0._DP
     DO  i = n-1,1, -1
        s(i) = v(i+1)*s(i+1) + w(i+1)
        t(i) = v(i+1)*t(i+1) + u(i+1)
     ENDDO
     !
     x(n) = (d(n) - b(n)*t(1) - c(n)*t(n-1))/ &
          (a(n) + b(n)*s(1) + c(n)*s(n-1))

     DO i = 1, n-1
        x(i) = s(i)*x(n) + t(i)
     ENDDO
     !
     RETURN
   END SUBROUTINE Progon3


   SUBROUTINE Spline_intrp (n,x,y,z,get_scoef,ib,ax,bx,xx,sp,dsp,d2sp)
     !***********************************************************
     !    Program Spline constructs the interpolating cubic     *
     !    spline for the function given in tabulated form       *
     !                                                          *
     !     n - number of knots                                  *
     !     x - array (of size n) containing the interpolation   *
     !         knots                                            *
     !     y - array (of size n) containing the function values *
     !         at knots                                         *
     !     z - array (of size n) containing the spline          *
     !         parameters                                       *
     !     get_scoef :                                          *
     !         true - program calculates the spline parameters  *
     !         false  - the spline parameters are known         *
     !    xx - point where values of splines and its            *
     !         derivatives are calculated                       *
     !    ib - type of boundary conditions                      *
     !               ind = i (i=1-4) - conditions of ith type   *
     !    ax, bx - parameter values in boundary conditions      *
     !         of the first or second type                      *
     !         Output:                                          *
     !    sp   - the spline value at point xx                   *
     !    dsp  - value of its first derivative at point xx      *
     !    d2sp - value of its second derivative at point xx     *

     !   INTERNAL (temporary arrays):
     !   a,b,                                                   *
     !   c,d - work arrays (of size n)                          *
     !***********************************************************
     USE nrtype,                        ONLY : DP,I4B
     IMPLICIT NONE
     INTEGER(I4B) n,ip
     REAL(DP)  x(n),y(n)
     REAL(DP)  a(n),b(n),c(n),d(n),z(n)
     REAL(DP) sp,dsp,d2sp,h1,xx,ax,bx,h2,am,al,g0,cc,gn,h,tt,rp,aa,bb
     INTEGER(I4b) ind,ib,ne,ns,nf,nj,j
     LOGICAL get_scoef
     !
     ind = 1
     IF(get_scoef) ind =0
     ip = 1
     IF (ind.EQ.0) THEN
        a(1) =  2.
        ne   =  n
        ns   =  2
        nf   =  n - 1
        SELECT CASE (ib)
           !
        CASE (1)
           b(1 ) = 0.
           c(1) =  0.
           d(1) =  2.*ax
           a(n) =  2.
           b(n) =  0.
           c(n) =  0.
           d(n) =  2.*bx
           !
        CASE (2)
           b(1) =  1.
           c(1) =  0.
           h1   =  x(2)-x(1)
           d(1) =  3.*(y(2)-y(1))/h1 - 0.5*h1*ax
           a(n) =  2.
           b(n) =  0.
           c(n) =  1.
           h1   =  x(n)-x(n-1)
           d(n) =  3.*(y(n)-y(n-1))/h1 + 0.5*h1*bx
           !
        CASE (3)
           h1   = x(2) - x(1)
           h2   = x(n) - x(n-1)
           am   = h2/(h1 + h2)
           al   = 1. - am
           b(1) = am
           c(1) = al
           d(1) = 3.*(am*(y(2) - y(1))/h1 + &
                al*(y(1) - y(n-1))/h2)
           h1   = x(n-1) - x(n-2)
           h2   = x(n) - x(n-1)
           am   = h1/(h1 + h2)
           al   = 1. - am
           a(n-1) = 2.
           b(n-1) = am
           c(n-1) = al
           d(n-1) = 3.*(am*(y(n) - y(n-1))/h2 + &
                al*(y(n-1) - y(n-2))/h1)
           nf = n - 2
           ne = n - 1
           !
        CASE (4)
           h1   = x(2)-x(1)
           h2   = x(3)-x(2)
           g0   = h1/h2
           a(2) = 1. + g0
           b(2) = g0
           c(2) = 0.
           am   = h1/(h1+h2)
           al   = 1.- am
           cc   = am*(y(3)-y(2))/h2 + al*(y(2)-y(1))/h1
           d(2) = cc + 2.*g0*(y(3)-y(2))/h2
           h2   = x(n)-x(n-1)
           h1   = x(n-1)-x(n-2)
           gn   = h1/h2
           a(n-1) = 1. + gn
           b(n-1) = 0.
           c(n-1) = gn
           am = h1/(h1+h2)
           al = 1. - am
           cc = am*(y(n)-y(n-1))/h2 +al*(y(n-1)-y(n-2))/h1
           d(n-1) = cc + 2.*gn*(y(n-1)-y(n-2))/h1
           ns = 3
           nf = n - 2
           ne = n - 2
           ip = 2
        END SELECT
        !
        DO  j  = ns ,nf
           h1   = x(j + 1) - x(j)
           h2  = x(j) - x(j-1)
           am   = h2/(h2 + h1)
           al   = 1. - am
           c(j) = al
           a(j) = 2.
           b(j) = am
           d(j) = 3.*(am*(y(j+1) - y(j))/h1 + &
                al*(y(j) - y(j-1))/h2)
        ENDDO
        !
        CALL Progon3 (a(ip),b(ip),c(ip),d(ip), &
             z(ip),ne)
        !
        SELECT CASE (ib)
           !
        CASE (3)
           z(n) = z(1)
           !
        CASE (4)
           z(1) = g0**2*z(3)+(g0**2-1.)*z(2)+ &
                2.*((y(2)-y(1))/(x(2)-x(1)) &
                -g0**2*(y(3)-y(2))/(x(3)-x(2)))

           z(n) = gn**2*z(n-2)+(gn**2-1.)*z(n-1)+ &
                2.*((y(n)-y(n-1))/(x(n)-x(n-1)) &
                -gn**2*(y(n-1)-y(n-2)) &
                /(x(n-1)-x(n-2)))
        END SELECT


        RETURN
     ENDIF

     
     !
     DO  j = 2,n
        nj= j
        IF(x(j).GT.xx) go to 1
     ENDDO
     !
1    j  = nj - 1
     h  = x(j+1) - x(j)
     tt = (xx - x(j))/h
     rp = (y(j+1) - y(j))/h
     aa = -2.*rp + z(j) + z(j+1)
     bb = -aa +rp - z(j)
     sp = y(j) + (xx - x(j))*(z(j) + tt*(bb + tt*aa))
     dsp= z(j) + tt*(bb + aa*tt) + tt*(bb + 2.*aa*tt)
     d2sp = (2.*bb + 6.*aa*tt)/h
     RETURN
   END SUBROUTINE Spline_intrp

 END MODULE spline
