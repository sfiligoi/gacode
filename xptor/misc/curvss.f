      subroutine curvss (n,x,y,d,isw,s,eps,ys,ysp,sigma,td,
     *                   tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,
     *                   ierr)
c
      implicit none
c
      integer i,n,isw,nm1,nm3,ibak,ierr
      real*8 p,sum,sl,su,rdim1,yspim2,sigmap,delxi1,delyi1,
     *      dim1,delxi,delyi,di,betap,betapp,alpha,alphap,
     *      beta,hsd1p,hdim1,hdi,rsd1i,rsd2i,f,g,wim1,wim2,
     *      tui,wi,h,step
      real*8 x(n),y(n),d(n),s,eps,ys(n),ysp(n),sigma,td(n),
     *     tsd1(n),hd(n),hsd1(n),hsd2(n),rd(n),rsd1(n),
     *     rsd2(n),v(n)
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine determines the parameters necessary to
c compute a smoothing spline under tension. for a given
c increasing sequence of abscissae (x(i)), i = 1,..., n and
c associated ordinates (y(i)), i = 1,..., n, the function
c determined minimizes the summation from i = 1 to n-1 of
c the square of the second derivative of f plus sigma
c squared times the difference of the first derivative of f
c and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
c functions f with two continuous derivatives such that the
c summation of the square of (f(x(i))-y(i))/d(i) is less
c than or equal to a given constant s, where (d(i)), i = 1,
c ..., n are a given set of observation weights. the
c function determined is a spline under tension with third
c derivative discontinuities at (x(i)), i = 2,..., n-1. for
c actual computation of points on the curve it is necessary
c to call the function curv2.
c
c on input--
c
c   n is the number of values to be smoothed (n.ge.2).
c
c   x is an array of the n increasing abscissae of the
c   values to be smoothed.
c
c   y is an array of the n ordinates of the values to be
c   smoothed, (i. e. y(k) is the functional value
c   corresponding to x(k) ).
c
c   d is a parameter containing the observation weights.
c   this may either be an array of length n or a scalar
c   (interpreted as a constant). the value of d
c   corresponding to the observation (x(k),y(k)) should
c   be an approximation to the standard deviation of error.
c
c   isw contains a switch indicating whether the parameter
c   d is to be considered a vector or a scalar,
c          = 0 if d is an array of length n,
c          = 1 if d is a scalar.
c
c   s contains the value controlling the smoothing. this
c   must be non-negative. for s equal to zero, the
c   subroutine does interpolation, larger values lead to
c   smoother funtions. if parameter d contains standard
c   deviation estimates, a reasonable value for s is
c   float(n).
c
c   eps contains a tolerance on the relative precision to
c   which s is to be interpreted. this must be greater than
c   or equal to zero and less than equal or equal to one. a
c   reasonable value for eps is sqrt(2./float(n)).
c
c   ys is an array of length at least n.
c
c   ysp is an array of length at least n.
c
c   sigma contains the tension factor. this value indicates
c   the degree to which the first derivative part of the
c   smoothing functional is emphasized. if sigma is nearly
c   zero (e. g. .001) the resulting curve is approximately a
c   cubic spline. if sigma is large (e. g. 50.) the
c   resulting curve is nearly a polygonal line. if sigma
c   equals zero a cubic spline results. a standard value for
c   sigma is approximately 1.
c
c and
c
c   td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, and v are
c   arrays of length at least n which are used for scratch
c   storage.
c
c on output--
c
c   ys contains the smoothed ordinate values.
c
c   ysp contains the values of the second derivative of the
c   smoothed curve at the given nodes.
c
c   ierr contains an error flag,
c        = 0 for normal return,
c        = 1 if n is less than 2,
c        = 2 if s is negative,
c        = 3 if eps is negative or greater than one,
c        = 4 if x-values are not strictly increasing,
c        = 5 if a d-value is non-positive.
c
c and
c
c   n, x, y, d, isw, s, eps, and sigma are unaltered.
c
c this subroutine references package modules terms and
c snhcsh.
c
c-----------------------------------------------------------
c
      if (n .lt. 2) go to 16
      if (s .lt. 0.) go to 17
      if (eps .lt. 0. .or. eps .gt. 1.) go to 18
      ierr = 0
      p = 0.D0
      v(1) = 0.D0
      v(n) = 0.D0
      ysp(1) = 0.D0
      ysp(n) = 0.D0
      if (n .eq. 2) go to 14
      rsd1(1) = 0.D0
      rd(1) = 0.D0
      rsd2(n) = 0.D0
      rdim1 = 0.D0
      yspim2 = 0.D0
c
c denormalize tension factor
c
      sigmap = abs(sigma)*float(n-1)/(x(n)-x(1))
c
c form t matrix and second differences of y into ys
c
      nm1 = n-1
      nm3 = n-3
      delxi1 = 1.D0
      delyi1 = 0.D0
      dim1 = 0.D0
      do 1 i = 1,nm1
        delxi = x(i+1)-x(i)
        if (delxi .le. 0.) go to 19
        delyi = (y(i+1)-y(i))/delxi
        ys(i) = delyi-delyi1
        call terms (di,tsd1(i+1),sigmap,delxi)
        td(i) = di+dim1
        hd(i) = -(1.D0/delxi+1./delxi1)
        hsd1(i+1) = 1.D0/delxi
        delxi1 = delxi
        delyi1 = delyi
    1   dim1 = di
c
c calculate lower and upper tolerances
c
      sl = s*(1.D0-eps)
      su = s*(1.D0+eps)
      if (isw .eq. 1) go to 3
c
c form h matrix - d array
c
      if (d(1) .le. 0. .or. d(2) .le. 0.) go to 20
      betapp = 0.D0
      betap = 0.D0
      alphap = 0.D0
      do 2 i = 2,nm1
        alpha = hd(i)*d(i)*d(i)
        if (d(i+1) .le. 0.) go to 20
        beta = hsd1(i+1)*d(i+1)*d(i+1)
        hd(i) = (hsd1(i)*d(i-1))**2+alpha*hd(i)
     *                             +beta*hsd1(i+1)
        hsd2(i) = hsd1(i)*betapp
        hsd1(i) = hsd1(i)*(alpha+alphap)
        alphap = alpha
        betapp = betap
    2   betap = beta
      go to 5
c
c form h matrix - d constant
c
    3 if (d(1) .le. 0.) go to 20
      sl = d(1)*d(1)*sl
      su = d(1)*d(1)*su
      hsd1p = 0.D0
      hdim1 = 0.D0
      do 4 i = 2,nm1
        hdi = hd(i)
        hd(i) = hsd1(i)*hsd1(i)+hdi*hdi+hsd1(i+1)*hsd1(i+1)
        hsd2(i) = hsd1(i)*hsd1p
        hsd1p = hsd1(i)
        hsd1(i) = hsd1p*(hdi+hdim1)
    4   hdim1 = hdi
c
c top of iteration
c cholesky factorization of p*t+h into r
c
    5 do 6 i = 2,nm1
        rsd2i = hsd2(i)
        rsd1i = p*tsd1(i)+hsd1(i)-rsd2i*rsd1(i-1)
        rsd2(i) = rsd2i*rdim1
        rdim1 = rd(i-1)
        rsd1(i) = rsd1i*rdim1
        rd(i) = 1.D0/(p*td(i)+hd(i)-rsd1i*rsd1(i)
     *                           -rsd2i*rsd2(i))
        ysp(i) = ys(i)-rsd1(i)*ysp(i-1)-rsd2(i)*yspim2
    6   yspim2 = ysp(i-1)
c
c back solve of r(transpose)* r * ysp = ys
c
      ysp(nm1) = rd(nm1)*ysp(nm1)
      if (n .eq. 3) go to 8
      do 7 ibak = 1,nm3
        i = nm1-ibak
    7   ysp(i) = rd(i)*ysp(i)-rsd1(i+1)*ysp(i+1)
     *                       -rsd2(i+2)*ysp(i+2)
    8 sum = 0.D0
      delyi1 = 0.D0
      if (isw .eq. 1) go to 10
c
c calculation of residual norm
c  - d array
c
      do 9 i = 1,nm1
        delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
        v(i) = (delyi-delyi1)*d(i)*d(i)
        sum = sum+v(i)*(delyi-delyi1)
    9   delyi1 = delyi
      v(n) = -delyi1*d(n)*d(n)
      go to 12
c
c calculation of residual norm
c  - d constant
c
   10 do 11 i = 1,nm1
        delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
        v(i) = delyi-delyi1
        sum = sum+v(i)*(delyi-delyi1)
   11   delyi1 = delyi
      v(n) = -delyi1
   12 sum = sum-v(n)*delyi1
c
c test for convergence
c
      if (sum .le. su) go to 14
c
c calculation of newton correction
c
      f = 0.D0
      g = 0.D0
      wim2 = 0.D0
      wim1 = 0.D0
      do 13 i = 2,nm1
        tui = tsd1(i)*ysp(i-1)+td(i)*ysp(i)
     *                        +tsd1(i+1)*ysp(i+1)
        wi = tui-rsd1(i)*wim1-rsd2(i)*wim2
        f = f+tui*ysp(i)
        g = g+wi*wi*rd(i)
        wim2 = wim1
   13   wim1 = wi
      h = f-p*g
      if (h .le. 0.) go to 14
c
c update p - newton step
c
      step = (sum-dsqrt(sum*sl))/h
      if (sl .ne. 0.) step = step*dsqrt(sum/sl)
      p = p+step
      go to 5
c
c store smoothed y-values and second derivatives
c
   14 do 15 i = 1,n
        ys(i) = y(i)-v(i)
   15   ysp(i) = p*ysp(i)
      return
c
c n less than 2
c
   16 ierr = 1
      return
c
c s negative
c
   17 ierr = 2
      return
c
c eps negative or greater than 1
c
   18 ierr = 3
      return
c
c x-values not strictly increasing
c
   19 ierr = 4
      return
c
c weight non-positive
c
   20 ierr = 5
      return
      end
c------------------------------------------------------------------
      subroutine snhcsh (sinhm,coshm,x,isw)
c
      implicit none
c
      integer isw
      real*8 sinhm,coshm,x,ax
      real*8 sp10,sp11,sp12,sp13,sp20,sp21,sp22,sp23,sp24,
     *       sp31,sp32,sp33,sp41,sp42,sp43,sq30,sq31,sq32,
     *       sq40,sq41,sq42,cp0,cp1,cp2,cp3,cp4,xs,expx
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine returns approximations to
c       sinhm(x) = sinh(x)/x-1
c       coshm(x) = cosh(x)-1
c and
c       coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)
c with relative error less than 1.0e-6
c
c on input--
c
c   x contains the value of the independent variable.
c
c   isw indicates the function desired
c           = -1 if only sinhm is desired,
c           =  0 if both sinhm and coshm are desired,
c           =  1 if only coshm is desired,
c           =  2 if only coshmm is desired,
c           =  3 if both sinhm and coshmm are desired.
c
c on output--
c
c   sinhm contains the value of sinhm(x) if isw .le. 0 or
c   isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.
c   2).
c
c   coshm contains the value of coshm(x) if isw .eq. 0 or
c   isw .eq. 1 and contains the value of coshmm(x) if isw
c   .ge. 2 (coshm is unaltered if isw .eq. -1).
c
c and
c
c   x and isw are unaltered.
c
c-----------------------------------------------------------
c
      data sp13/.3029390D-5/,
     *     sp12/.1975135D-3/,
     *     sp11/.8334261D-2/,
     *     sp10/.1666665D0/
      data sp24/.3693467D-7/,
     *     sp23/.2459974D-5/,
     *     sp22/.2018107D-3/,
     *     sp21/.8315072D-2/,
     *     sp20/.1667035D0/
      data sp33/.6666558D-5/,
     *     sp32/.6646307D-3/,
     *     sp31/.4001477D-1/,
     *     sq32/.2037930D-3/,
     *     sq31/-.6372739D-1/,
     *     sq30/.6017497D1/
      data sp43/.2311816D-4/,
     *     sp42/.2729702D-3/,
     *     sp41/.9868757D-1/,
     *     sq42/.1776637D-3/,
     *     sq41/-.7549779D-1/,
     *     sq40/.9110034D1/
      data cp4/.2982628D-6/,
     *     cp3/.2472673D-4/,
     *     cp2/.1388967D-2/,
     *     cp1/.4166665D-1/,
     *     cp0/.5000000D0/
c
      ax = dabs(x)
      if (isw .ge. 0) go to 5
c
c sinhm approximation
c
      if (ax .gt. 4.45) go to 2
      xs = ax*ax
      if (ax .gt. 2.3) go to 1
c
c sinhm approximation on (0.,2.3)
c
      sinhm = xs*(((sp13*xs+sp12)*xs+sp11)*xs+sp10)
      return
c
c sinhm approximation on (2.3,4.45)
c
    1 sinhm = xs*((((sp24*xs+sp23)*xs+sp22)*xs+sp21)
     .               *xs+sp20)
      return
    2 if (ax .gt. 7.65) go to 3
c
c sinhm approximation on (4.45,7.65)
c
      xs = ax*ax
      sinhm = xs*(((sp33*xs+sp32)*xs+sp31)*xs+1.)/
     .             ((sq32*xs+sq31)*xs+sq30)
      return
    3 if (ax .gt. 10.1) go to 4
c
c sinhm approximation on (7.65,10.1)
c
      xs = ax*ax
      sinhm = xs*(((sp43*xs+sp42)*xs+sp41)*xs+1.)/
     .             ((sq42*xs+sq41)*xs+sq40)
      return
c
c sinhm approximation above 10.1
c
    4 sinhm = dexp(ax)/(ax+ax)-1.D0
      return
c
c coshm and (possibly) sinhm approximation
c
    5 if (isw .ge. 2) go to 7
      if (ax .gt. 2.3) go to 6
      xs = ax*ax
      coshm = xs*((((cp4*xs+cp3)*xs+cp2)*xs+cp1)*xs+cp0)
      if (isw .eq. 0) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)
     .                              *xs+sp10)
      return
    6 expx = dexp(ax)
      coshm = (expx+1.D0/expx)/2.D0-1.D0
      if (isw .eq. 0) sinhm = (expx-1.D0/expx)/(ax+ax)-1.D0
      return
c
c coshmm and (possibly) sinhm approximation
c
    7 xs = ax*ax
      if (ax .gt. 2.3) go to 8
      coshm = xs*(((cp4*xs+cp3)*xs+cp2)*xs+cp1)
      if (isw .eq. 3) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)
     .                              *xs+sp10)
      return
    8 expx = dexp(ax)
      coshm = ((expx+1.D0/expx-xs)/2.D0-1.D0)/xs
      if (isw .eq. 3) sinhm = (expx-1.D0/expx)/(ax+ax)-1.D0
      return
      end
c------------------------------------------------------------------
      subroutine terms (diag,sdiag,sigma,del)
c
      implicit none
c
      real*8 diag,sdiag,sigma,del,sigdel,denom,coshm,sinhm
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine computes the diagonal and superdiagonal
c terms of the tridiagonal linear system associated with
c spline under tension interpolation.
c
c on input--
c
c   sigma contains the tension factor.
c
c and
c
c   del contains the step size.
c
c on output--
c
c                sigma*del*cosh(sigma*del) - sinh(sigma*del)
c   diag = del*--------------------------------------------.
c                     (sigma*del)**2 * sinh(sigma*del)
c
c                   sinh(sigma*del) - sigma*del
c   sdiag = del*----------------------------------.
c                (sigma*del)**2 * sinh(sigma*del)
c
c and
c
c   sigma and del are unaltered.
c
c this subroutine references package module snhcsh.
c
c-----------------------------------------------------------
c
      if (sigma .ne. 0.) go to 1
      diag = del/3.D0
      sdiag = del/6.D0
      return
    1 sigdel = sigma*del
      call snhcsh (sinhm,coshm,sigdel,0)
      denom = sigma*sigdel*(1.D0+sinhm)
      diag = (coshm-sinhm)/denom
      sdiag = sinhm/denom
      return
      end
