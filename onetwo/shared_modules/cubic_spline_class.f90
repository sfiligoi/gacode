 MODULE cubic_spline
   USE nrtype,                             ONLY : DP,I4B  
   USE common_constants,                   ONLY : zeroc,izero

   IMPLICIT NONE
   PRIVATE


   TYPE spl1d_inst
      REAL(DP),DIMENSION(:) ,POINTER   :: s_coef,x,y
      REAL(DP) ax,bx,sp,dsp,d2sp,xx
      INTEGER(I4B) ib,n
   END TYPE spl1d_inst


   PUBLIC :: eval_1d_cubic_spline
   PUBLIC :: create_1d_spline_object
   PUBLIC :: done_1dspline
   PUBLIC :: spl1d_inst
   PUBLIC :: spline_interp
 

 
   CONTAINS







   SUBROUTINE  create_1d_spline_object(spline_inst,n,x,y,ax,bx,ib)

     !***********************************************************
     !    Subroutine  constructs the interpolating cubic        *
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
     !       ib = i (i=1-6) - conditions of ith type            *
     !                                                          *
     !       ib = 1:                                            *
     !          ax = derivative of spline at x(1)               *
     !          bx = derivative of spline at x(n)               *
     !       ib = 2:                                            *
     !          ax = second derivative at x(1)                  *
     !          bx = second derivative at x(n)                  *
     !       ib = 3:                                            *
     !          first deriv at x(1) = first deriv at x(n)       *
     !          second  deriv at x(1) = second  deriv at x(n)   *
     !       ib = 4:                                            *
     !          third derivative continuous at x(1) and x(n)    *
     !       ib = 5:                                            *
     !          third derivative continuous at x(2) and x(n-1)  *
     !          (this is a replacement for icsccu, icsevu)      *
     !       ib = 6:                                            *
     !          ax = derivative of spline at x(1)               *
     !          bx = second derivative at x(n)                  *
     !                                                          *
     !    ax, bx - parameter values in boundary conditions      *
     !         of the first or second type,                     *
     !         not used for ib=3,4,5                            *
     !                                                          *
     !   Output:                                                *
     !    sp   - the spline value at point xx                   *
     !    dsp  - value of its first derivative at point xx      *
     !    d2sp - value of its second derivative at point xx     *
     !                                                          *
     !   INTERNAL (temporary arrays):                           *
     !   a,b,                                                   *
     !   c,d - work arrays (of size n)                          *
     !***********************************************************
     USE nrtype,                        ONLY : DP,I4B
     IMPLICIT NONE
     TYPE(spl1d_inst) spline_inst  
     INTEGER(I4B) ip,n
     REAL(DP) sp,dsp,d2sp,h1,xx,ax,bx,h2,am,al,g0,cc,gn,h,tt,rp,aa,bb


     REAL(DP), DIMENSION(n) :: a,b,c,d,x,y,s_coef
     INTEGER(I4b) ind,ib,ne,ns,nf,nj,j,alloc_stat

 
 
    
     spline_inst%n = n ; spline_inst%ib = ib

     IF(ASSOCIATED(spline_inst%x))DEALLOCATE(spline_inst%x)
          ALLOCATE(spline_inst%x(spline_inst%n),stat=alloc_stat)
     if (alloc_stat .ne. 0) then 
       STOP  'Allocation unsuccessful in create_1d_spline_object'
     endif
     IF(ASSOCIATED(spline_inst%y))DEALLOCATE(spline_inst%y)
          ALLOCATE(spline_inst%y(spline_inst%n),stat=alloc_stat)
     if (alloc_stat .ne. 0) then 
       STOP  'Allocation unsuccessful in create_1d_spline_object'
     endif
     IF(ASSOCIATED(spline_inst%s_coef))DEALLOCATE(spline_inst%s_coef)
          ALLOCATE(spline_inst%s_coef(spline_inst%n),stat=alloc_stat)
     if (alloc_stat .ne. 0) then 
       STOP  'Allocation unsuccessful in create_1d_spline_object'
     endif
      
     spline_inst%ax   = ax      ;      spline_inst%bx   = bx
     spline_inst%x(:) = x(:)    ;      spline_inst%y(:) = y(:)

     ip = 1
        a(1) =  2._DP
        ne   =  spline_inst%n
        ns   =  2
        nf   =  spline_inst%n - 1
        SELECT CASE (ib)
           !
        CASE (1)
           b(1 ) = zeroc
           c(1) =  zeroc
           d(1) =  2._DP*spline_inst%ax
           a(n) =  2._DP
           b(n) =  zeroc
           c(n) =  zeroc
           d(n) =  2._DP*spline_inst%bx
           !
        CASE (2)
           b(1) =  1._DP
           c(1) =  zeroc
           h1   =  spline_inst%x(2)-spline_inst%x(1)
           d(1) =  3._DP*(spline_inst%y(2)-spline_inst%y(1))/h1 - 0.5_DP*h1*ax
           a(n) =  2._DP
           b(n) =  zeroc
           c(n) =  1._DP
           h1   =  spline_inst%x(n)-spline_inst%x(n-1)
           d(n) =  3._DP*(spline_inst%y(n)-spline_inst%y(n-1))/h1 + 0.5_DP*h1*bx
           !
        CASE (3)
           h1   = spline_inst%x(2) - spline_inst%x(1)
           h2   = spline_inst%x(n) - spline_inst%x(n-1)
           am   = h2/(h1 + h2)
           al   = 1._DP - am
           b(1) = am
           c(1) = al
           d(1) = 3._DP*(am*(spline_inst%y(2) - spline_inst%y(1))/h1 + &
                al*(spline_inst%y(1) - spline_inst%y(n-1))/h2)
           h1   = spline_inst%x(n-1) - spline_inst%x(n-2)
           h2   = spline_inst%x(n) - spline_inst%x(n-1)
           am   = h1/(h1 + h2)
           al   = 1._DP  - am
           a(n-1) = 2._DP
           b(n-1) = am
           c(n-1) = al
           d(n-1) = 3._DP *(am*(spline_inst%y(n) - spline_inst%y(n-1))/h2 + &
                al*(spline_inst%y(n-1) - spline_inst%y(n-2))/h1)
           nf = n - 2
           ne = n - 1
           !
        CASE (4)
           h1   = spline_inst%x(2)-spline_inst%x(1)
           h2   = spline_inst%x(3)-spline_inst%x(2)
           g0   = h1/h2
           a(2) = 1._DP + g0
           b(2) = g0
           c(2) = zeroc
           am   = h1/(h1+h2)
           al   = 1._DP- am
           cc   = am*(spline_inst%y(3)-spline_inst%y(2))/h2 + al*(spline_inst%y(2)-spline_inst%y(1))/h1
           d(2) = cc + 2._DP*g0*(spline_inst%y(3)-spline_inst%y(2))/h2
           h2   = spline_inst%x(n)-spline_inst%x(n-1)
           h1   = spline_inst%x(n-1)-spline_inst%x(n-2)
           gn   = h1/h2
           a(n-1) = 1._DP + gn
           b(n-1) = zeroc
           c(n-1) = gn
           am = h1/(h1+h2)
           al = 1._DP - am
           cc = am*(spline_inst%y(n)-spline_inst%y(n-1))/h2 +al*(spline_inst%y(n-1)-spline_inst%y(n-2))/h1
           d(n-1) = cc + 2._DP*gn*(spline_inst%y(n-1)-spline_inst%y(n-2))/h1
           ns = 3
           nf = n - 2
           ne = n - 2 
           ip = 2
           !
        CASE (5)
           h2 = x(2)-x(1)
           h1 = x(3)-x(2)
           a(1) = 1._DP+h2/h1
           b(1) = 2._DP*(h1+h2)*h2/h1**2+1._DP-h2**2/h1**2
           c(1) = zeroc
           d(1) = 2._DP*( (y(2)-y(1))/h2 - (y(3)-y(2))*h2**2/h1**3 ) + &
                  3._DP*( (y(3)-y(2))*h2**2/h1**3 + (y(2)-y(1))/h1 )
           h2 = x(n-1)-x(n-2)
           h1 = x(n) - x(n-1)
           a(n) = -(1._DP/h1**2+1._DP/(h1*h2))
           b(n) = zeroc
           c(n) = 1._DP/h2**2-1._DP/h1**2-2._DP*(h2+h1)/(h1*h2**2)
           d(n) = 2._DP*( (y(n-1)-y(n-2))/h2**3-(y(n)-y(n-1))/h1**3 ) - &
                    3._DP/h2**2*(h2/h1**2*(y(n)-y(n-1))+(y(n-1)-y(n-2))/h2)
           ns = 2
           nf = n-1
           ip = 1
           ne = n
           !
        CASE (6)
           !First derivative b.c. at x(1)
           b(1) =  zeroc
           c(1) =  zeroc
           d(1) =  2._DP*ax
           !Second derivative b.c. at x(n)
           a(n) =  2._DP
           b(n) =  zeroc
           c(n) =  1._DP
           h1   =  x(n)-x(n-1)
           d(n) =  3._DP*(y(n)-y(n-1))/h1 + 0.5_DP*h1*bx
           !
        END SELECT
        !
        DO  j  = ns ,nf
           h1   = spline_inst%x(j + 1) - spline_inst%x(j)
           h2   = spline_inst%x(j) - spline_inst%x(j-1)
           am   = h2/(h2 + h1)
           al   = 1._DP - am
           c(j) = al
           a(j) = 2._DP
           b(j) = am
           d(j) = 3._DP*(am*(spline_inst%y(j+1) - spline_inst%y(j))/h1 + &
                al*(spline_inst%y(j) - spline_inst%y(j-1))/h2)
        ENDDO
        !
        CALL Progon3 (a(ip),b(ip),c(ip),d(ip),s_coef(ip),ne)
        spline_inst%s_coef(:) = s_coef(:)

        SELECT CASE (ib)
           !
        CASE (3)
           spline_inst%s_coef(n) = spline_inst%s_coef(1)
           !
        CASE (4)
           spline_inst%s_coef(1) = g0**2*spline_inst%s_coef(3)+(g0**2-1.)*spline_inst%s_coef(2)+ &
                2._DP*((spline_inst%y(2)-spline_inst%y(1))/(spline_inst%x(2)-spline_inst%x(1)) &
                -g0**2*(spline_inst%y(3)-spline_inst%y(2))/(spline_inst%x(3)-spline_inst%x(2)))

           spline_inst%s_coef(n) = gn**2*spline_inst%s_coef(n-2)+(gn**2-1.)*spline_inst%s_coef(n-1)+ &
                2._DP*((spline_inst%y(n)-spline_inst%y(n-1))/(spline_inst%x(n)-spline_inst%x(n-1)) &
                -gn**2*(spline_inst%y(n-1)-spline_inst%y(n-2)) &
                /(spline_inst%x(n-1)-spline_inst%x(n-2)))
        END SELECT


        RETURN
     END SUBROUTINE create_1d_spline_object


!   SUBROUTINE  eval_1d_cubic_spline (n,x,y,get_scoef,ib,ax,bx,xx,sp,dsp,d2sp)
   SUBROUTINE  eval_1d_cubic_spline(spline_inst)

     !***********************************************************
     !    Subroutine  constructs the interpolating cubic     *
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
     !       {included here for reference,bc must be set when   *
     !       spline is created using                            *
     !       SUBROUTINE  create_1d_spline_object}               *
     !       ib = i (i=1-4) - conditions of ith type            *
     !                                                          *
     !       ib = 1:                                            *
     !          ax = derivative of spline at x(1)               *
     !          bx = derivative of spline at x(n)               *
     !       ib = 2:                                            *
     !          ax = second derivtive at x(1)                   *
     !          bx = second derivtive at x(n)                   *
     !       ib = 3:                                            *
     !          first deriv at x(1) = first deriv at x(n)       *
     !          second  deriv at x(1) = second  deriv at x(n)   *
     !       ib = 4:                                            *
     !          third derivative continuous at x(1) and x(n)    *
     !                                                          *
     !    ax, bx - parameter values in boundary conditions      *
     !         of the first or second type,                     *
     !         not used for ib=3,4                              *
     !                                                          *     !         
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
     TYPE(spl1d_inst) spline_inst  
     INTEGER(I4B) ip,n
     REAL(DP) sp,dsp,d2sp,h1,xx,ax,bx,h2,am,al,g0,cc,gn,h,tt,rp,aa,bb
     INTEGER(I4b) ind,ib,ne,ns,nf,nj,j
 
 
  
          
     ax = spline_inst%ax ; bx  =  spline_inst%bx ; n =  spline_inst%n
     ib = spline_inst%ib ; xx  =  spline_inst%xx

     DO  j = 2,spline_inst%n
        nj= j
        IF(spline_inst%x(j) .GT. xx) go to 1
     ENDDO
     !
1    j  = nj - 1
     h  = spline_inst%x(j+1) - spline_inst%x(j)
     tt = (xx - spline_inst%x(j))/h
     rp = (spline_inst%y(j+1) - spline_inst%y(j))/h
     aa = -2._DP*rp + spline_inst%s_coef(j) + spline_inst%s_coef(j+1)
     bb = -aa +rp - spline_inst%s_coef(j)
     IF (spline_inst%x(n).eq.xx) THEN
       spline_inst%sp = spline_inst%y(n) 
     ELSE
       spline_inst%sp = spline_inst%y(j) + (xx - spline_inst%x(j))*(spline_inst%s_coef(j) + tt*(bb + tt*aa))
     ENDIF
     spline_inst%dsp= spline_inst%s_coef(j) + tt*(bb + aa*tt) + tt*(bb + 2._DP*aa*tt)
     spline_inst%d2sp = (2._DP*bb + 6.*aa*tt)/h

     RETURN
   END SUBROUTINE eval_1d_cubic_spline


       SUBROUTINE done_1dspline(spline_inst)   ! destructor
!-----------------------------------------------------------
       TYPE(spl1d_inst) spline_inst 

            IF(ASSOCIATED(spline_inst%s_coef))DEALLOCATE(spline_inst%s_coef)
            IF(ASSOCIATED(spline_inst%x))DEALLOCATE(spline_inst%x)
            IF(ASSOCIATED(spline_inst%y))DEALLOCATE(spline_inst%y)
            NULLIFY(spline_inst%y,spline_inst%x,spline_inst%s_coef)
       RETURN
       END SUBROUTINE done_1dspline



      FUNCTION spline_interp(xdata,ydata,xint,spline_type_in,ax_in,bx_in,coef)
        !Given the data points at (xdata,ydata), return the spline 
        !interpolation at the points xint
        !OPTIONAL ARGUMENTS [default]:
        ! spline_type_in [5] - The type of spline (argument ib for 
        !                         create_1d_spline_object, see above)
        ! ax_in [0] - The left boundary condition (see create_1d_spline_object above)
        ! bx_in [0] - The right boundary condition (see create_1d_spline_object above)
        ! coef [output] - Length ndata array holding the spline coefficients as output
        ! 
        !The default boundary conditions are not-a-knot boundary conditions
        REAL(DP), INTENT(IN), DIMENSION(:) :: xdata, ydata, xint
        INTEGER, INTENT(IN), OPTIONAL :: spline_type_in
        REAL(DP), INTENT(IN), OPTIONAL :: ax_in, bx_in
        REAL(DP), INTENT(OUT), DIMENSION(SIZE(xdata)),OPTIONAL :: coef
        REAL(DP), DIMENSION(SIZE(xint)) :: spline_interp, yint2
        REAL(DP), DIMENSION(SIZE(xint),3) :: cspln
        TYPE(spl1d_inst) :: spline
        INTEGER :: nint, ndata, i, ier, spline_type
        REAL(DP) :: ax, bx
        REAL(DP), DIMENSION(SIZE(xdata)) :: xtmp,ytmp
        spline_type = 5
        ax = 0.
        bx = 0.
        IF (PRESENT(spline_type_in)) spline_type = spline_type_in
        IF (PRESENT(ax_in)) ax = ax_in
        IF (PRESENT(bx_in)) bx = bx_in
        nint = SIZE(xint)
        ndata = SIZE(xdata)
        xtmp = xdata
        ytmp = ydata
        if (xdata(2)<xdata(1)) then
          xtmp = xdata(ndata:1:-1)
          ytmp = ydata(ndata:1:-1)
        endif
!        call icsccu(xdata,ydata,ndata,cspln,nint,ier)
!        if (ier .ne. 0) then
!          print *,xdata,ydata
!          STOP 'Error in setting up spline_interp'
!        endif
!        call icsevu(xdata,ydata,ndata,cspln,nint,xint,yint2,nint,ier)
!        if (ier .ne. 0) then
!          print *,xint
!             STOP 'Error in evaluating spline_interp'
!        endif
        call create_1d_spline_object(spline,ndata,xtmp,ytmp,ax,bx,spline_type)
!        print *,spline%s_coef
        DO i=1,nint
          spline%xx = xint(i)
          call eval_1d_cubic_spline(spline)
          spline_interp(i) = spline%sp
!          print '(i4,3G20.8)',i,spline_interp(i),yint2(i),spline_interp(i)-yint2(i)
        ENDDO
        IF (xtmp(1).eq.xint(1)) spline_interp(1) = ytmp(1)
        IF (xtmp(ndata).eq.xint(nint)) spline_interp(nint) = ytmp(ndata)
        IF (present(coef))  coef = spline%s_coef
        call done_1dspline(spline)
        RETURN 
      END FUNCTION !spline_interp

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
     !        tEMPORARY ARRAYS:
     REAL(DP)  a(n),b(n),c(n),d(n),x(n)
     REAL(DP)  u(n+1),v(n+1),w(n+1),s(n+1),t(n+1)
     REAL(DP) z

     u(1) = zeroc
     v(1) = zeroc
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
     t(n) = zeroc
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


 END MODULE cubic_spline
