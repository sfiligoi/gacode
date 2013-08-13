

      subroutine cntour (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin,
     .                   yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,
     .                   xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,
     .                   nh,cspln,n2cspln,nh2,itty,iptsm,ierr,
     .                   bpmag,iconvg,delta_psi)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- INTERFACE TO CONTOUR ROUTINES ---------------------- HSJ ---------
c
      include 'shape.i'
c
      logical   fl1, fl2
      dimension xc(*), yc(*), x(nw), y(nh),
     .          cspln(n2cspln,nw,nh2), bpmag(*)
c
      fl1        = use_cnt1
      fl2        = use_cnt2
      iauto_save = iauto
c
   10 ierr       = 0
      iauto      = iauto_save
      iautoc     = 0
c
      if (use_cnt1) then
        call cntour1 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin,
     .                yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,
     .                xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,
     .                nh,cspln,n2cspln,nh2,itty,iptsm,ierr,bpmag,iconvg,
     .                delta_psi)
c
      else if (use_cnt2) then
        call cntour2 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin,
     .                yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,
     .                xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,
     .                nh,cspln,n2cspln,nh2,itty,iptsm,ierr,bpmag,iconvg,
     .                delta_psi)
      else
        call STOP ('subroutine CNTOUR: use_cnt problem', 164)
      end if
c
      if (ierr .ne. 0) then
        use_cnt1 = .not. use_cnt1
        use_cnt2 = .not. use_cnt2
        if (fl1 .neqv. use_cnt1)  go to 10
        call interface_dump_psi (psivl)     ! dump the psi values
        print *,'IERR =',IERR
        call STOP ('subroutine CNTOUR: IERR is not zero', 261)
      end if
c
      use_cnt1 = fl1         ! restore to state at entry
      use_cnt2 = fl2
      return
c
      end

      subroutine coefcalc (c, all, alt, alr, alb, deltar, deltarzsq,
     .                     deltarsq, rl, interior)
c
      implicit none
c
c ----------------------------------------------------------------------
c
c --- subroutine COEFCALC returns the coefficients of the unknowns
c     surrounding
c
c --- unknown psil:
c    coefficent of psill is c(1)
c                  psilt    c(2)
c                  psilr    c(3)
c                  psilb    c(4)
c                  rhs(l)   c(5)
c    here psil is the central pt.psill is point to left,psilt
c    is point on top,psilr is point to right,psilb is point at
c    bottom. c(5) is coefficient of source term.
c    see notes,pg.144,mhdcalc notebook.
c
c --- input ---
c
c  deltar                 r grid spacing
c  deltarsq               deltar*deltar
c  deltarzsq               (deltar/deltaz)**2
c  all                    fraction (0 .lt. all .le. 1) of
c                         deltar that point psill is to the left
c  alr                    same for point to the right
c  alt                    fraction (0 .lt. alt .le. 1) of
c                         deltaz that point psilt is on top
c  alb                                   on bottom
c  rl                     r coordinate of point psil
c  zl                     z
c
c --- output ---
c
c  c(1-9)                 as explained above
c
c ---------------------------------------------------------- HSJ 3/26/93
c
      logical  interior
      real*8   c(*), alzz, alrr, ac, rl, all, alt, alr, alb,
     .         deltar, deltarsq, deltarzsq
c
      alzz = alb*alt*(alb+alt)
      alrr = alr*all*(all+alr)
      ac   =    (-2.0*(all+alr)+(deltar/rl)*(all*all-alr*alr)) / alrr
      ac   = ac - 2.0*(alt+alb)*deltarzsq/alzz
c
c --- coefficient of point to left, c(1)
c
      c(1) = (alr/alrr) * (2.0+alr*deltar/rl) / ac
c
c --- coefficient of point on top, c(2)
c
      c(2) = 2.0 * alb*deltarzsq/alzz/ac
c
c --- coefficient of point to right, c(3)
c
      c(3) = (all/alrr)*(2.0-deltar*all/rl)/ac
c
c --- coefficient of point on bottom, c(4)
c
      c(4) = 2.0 * alt*deltarzsq/alzz/ac
c
      if (.not. interior)  return
c
c --- coefficient of source, c(5)
c
      c(5) = deltarsq/ac
      c(6) = all
      c(7) = alt
      c(8) = alr
      c(9) = alb
c
      return
c
      end

      subroutine dump_psi_values (psi, nw, nh, rmhdgrid, zmhdgrid,
     .                            iounit, isignpsi, zero, isym,
     .                            map, mapset,psif,mf,psitrace)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- DUMP_PSI_VALUES dumps out the rmhdgrid, zmhdgrid and psi values
c --- for examination by another program if a fatal error occurred. HSJ
c
      logical  mapset
      integer  map(*)
      real*8   psi(nw,nh), rmhdgrid(nw), zmhdgrid(nh), zero(*),
     .         psif(*)
c
      call getioun(iounit,iounit)
      open (unit = iounit, file = 'dump_psi.dat', status = 'UNKNOWN')
c
      write  (iounit, 10)  nw, nh, isignpsi, isym, mapset,mf
   10 format (2x, 4(2x, i5), l8,i4)
      write  (iounit, 20)  (rmhdgrid(i), i=1,nw)
      write  (iounit, 20)  (zmhdgrid(i), i=1,nh)
      write  (iounit, 20)  ((psi(i,j),   i=1,nw), j=1,nh)
      write  (iounit, 20)  (zero(j), j=1,nw*nh)
   20 format (6(2x, 1pe20.12))
      if (mapset)
     .  write  (iounit, 30)  (map(j), j=1,nw*nh)
   30   format (6(2x, i5))
      if (mf .gt. 0)  write (iounit, 20) (psif(i), i=1,mf)
      write (iounit, 20) psitrace
      call giveupus(iounit)
      close (unit = iounit)
c
      write (6, *)  'file "dump_psi.dat" has current psi values in it'
      write (6, *)  'you may use program PSIPLOT to examine this file'
      return
c
      end

      subroutine findintsc (rl, zl, rp, zp, np, deltar, deltaz,
     .                      ll, lt, lr, lb, all, alt, alr, alb,
     .                      updownsym, jstart)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- input ---
c
c  rp(j)       r values of plasma contour
c  zp(j)       z
c  np          number of points in rp,zp
c              note the contour must be closed (i.e., rp(1) = rp(np)
c                                                 and zp(1) = zp(np))
c  rl           r coordinate of central point
c  zl           z coordinate of central point
c  ll            if ll = 0 then find the intersection at the left
c                in this case the output will be all
c  lr            if lr = 0 then find the intersection at the right
c                in this case the output will be alr
c  lt            if lt = 0 then find the intersection at the top
c                in this case the output will be alt
c  lb            if lb = 0 then find the intersection at the bottom
c                in this case the output will be alb
c
c updownsym          set to true if up/down symmetric solution
c                    is to be generated.
c
c  jstart         set to true if updownsym is true and this is
c                 the first row of unknowns. otherwise set to false.
c --- output ---
c
c  all             fraction,(0 .lt. all .le. 1.0)of distance deltar
c                  to boundary at the left of point l.
c                  all  = 1.0 if point ll  is inside
c  alr             fraction,(0 .lt. all .le. 1.0)of distance deltar
c                  to boundary at the right of point l.
c                  alr  = 1.0 if point lr  is inside
c  alt             fraction,(0 .lt. all .le. 1.0)of distance deltaz
c                  to boundary on top  of point l.
c                  alt  = 1.0 if point lt  is inside
c  alb             fraction,(0 .lt. all .le. 1.0)of distance deltaz
c                  to boundary on bottom  of point l.
c                  alb  = 1.0 if point lb  is inside
c ---------------------------------------------------------- HSJ 3/26/92
c
      logical  updownsym, jstart
      integer  j, np, ll, lt, lr, lb
      real*8   all, alt, alr, alb, rl, zl, rp(*), zp(*), deltaz, deltar,
     .         rf, rb, zf, zb, dz, dr, slope, binc, rinctp, zinctp
c
c find intersection to the left
c
      if (ll .eq. 0) then
        all = 0.0
            do j=2,np
                rf = rl-rp(j)
                rb = rl-rp(j-1)
                if (rf .lt. 0.0 .and. rb .lt. 0.0)  go to 10
                zf = zl-zp(j)
                zb = zl-zp(j-1)
                if (zb*zf .le. 0.0) then
                    dz = zp(j)-zp(j-1)
                    dr = rp(j)-rp(j-1)
                    if (dr .eq. 0.0) then
                      all = (rl-rp(j))/deltar
                      if (all .ge. 0.0 .and.  all .le. 1.0)  go to 100
                      all = 0.0
                      go to 10
                    else
                      slope = dz/dr
                      if (slope .eq. 0.0)
     .            call STOP ('subroutine FINDINTSC: zero slope #1', 109)
                      binc   = zp(j)-slope*rp(j)
                      rinctp = (zl-binc)/slope
                      all    = (rl-rinctp)/deltar
                      if (all .ge. 0.0 .and. all .le. 1.0)  go to 100
                      all = 0.0
                      go to 10
                    end if
                end if
   10           continue
            end do
        call STOP ('subroutine FINDINTSC: no intersection #1', 110)
      else
        all = 1.0    ! point to left, ll, is inside
      end if
c
c find intersection to the right
c
  100 if (lr .eq. 0) then
        alr = 0.0
            do j=2,np
                rf = rl-rp(j)
                rb = rl-rp(j-1)
                if (rf .gt. 0.0 .and. rb .gt. 0.0)  go to 20
                zf = zl-zp(j)
                zb = zl-zp(j-1)
                if (zb*zf .le. 0.0) then
                    dz = zp(j)-zp(j-1)
                    dr = rp(j)-rp(j-1)
                    if (dr .eq. 0.0) then
                      alr = (rp(j)-rl)/deltar
                      if (alr .ge. 0.0 .and.  alr .le. 1.0)  go to 110
                      alr = 0.0
                      go to 20
                    else
                      slope = dz/dr
                      if (slope .eq. 0.0)
     .            call STOP ('subroutine FINDINTSC: zero slope #2', 111)
                      binc = zp(j)-slope*rp(j)
                      rinctp = (zl-binc)/slope
                      alr = (rinctp-rl)/deltar
                      if (alr .ge. 0.0 .and. alr .le. 1.0)  go to 110
                      alr = 0.0
                      go to 20
                    end if
                end if
   20           continue
            end do
        call STOP ('subroutine FINDINTSC: no intersection #2', 112)
      else
        alr = 1.0     ! point to right, lr, is inside
      end if
c
c --- find intersection at top
c
  110 if (lt .eq. 0) then
        alt = 0.0
            do j=2,np
                zf = zl-zp(j)
                zb = zl-zp(j-1)
                if (zf .gt. 0.0 .and. zb .gt. 0.0)  go to 30
                rf = rl-rp(j)
                rb = rl-rp(j-1)
                if (rb*rf .le. 0.0) then
                    dz = zp(j)-zp(j-1)
                    if (dz .eq. 0.0) then
                      alt = (zp(j)-zl)/deltaz
                      if (alt .ge. 0.0 .and.  alt .le. 1.0)  go to 120
                      alt = 0.0
                      go to 30
                    else
                      dr = rp(j)-rp(j-1)
                      if (dr .eq. 0.0)
     .        call STOP ('subroutine FINDINTSC: infinite slope #1', 113)
                      slope = dz/dr
                      binc = zp(j)-slope*rp(j)
                      zinctp = slope*rl+binc
                      alt = (zinctp-zl)/deltaz
                      if (alt .ge. 0.0 .and. alt .le. 1.0)  go to 120
                      alt = 0.0
                      go to 30
                    end if
                end if
   30           continue
            end do
        call STOP ('subroutine FINDINTSC: no intersection #3', 114)
      else
        alt = 1.0     ! point on top, lt, is inside
      end if
c
c --- find intersection at bottom
c
  120 if (updownsym .and. jstart)  go to 200
      if (lb .eq. 0) then
        alb = 0.0
            do j=2,np
                zf = zl-zp(j)
                zb = zl-zp(j-1)
                if (zf .lt. 0.0 .and. zb .lt. 0.0)  go to 40
                rf = rl-rp(j)
                rb = rl-rp(j-1)
                if (rb*rf .le. 0.0) then
                    dz = zp(j)-zp(j-1)
                    if (dz .eq. 0.0) then
                      alb = (zl-zp(j))/deltaz
                      if (alb .ge. 0.0 .and.  alb .le. 1.0)  go to 250
                      alb = 0.0
                      go to 40
                    else
                      dr = rp(j)-rp(j-1)
                      if (dr .eq. 0.0)
     .        call STOP ('subroutine FINDINTSC: infinite slope #2', 115)
                      slope = dz/dr
                      binc = zp(j)-slope*rp(j)
                      zinctp = slope*rl+binc
                      alb = (zl-zinctp)/deltaz
                      if (alb .ge. 0.0 .and. alb .le. 1.0)  go to 250
                      alb = 0.0
                      go to 40
                    end if
                end if
   40           continue
            end do
        call STOP ('subroutine FINDINTSC: no intersection #4', 116)
      else
        alb = 1.0     ! point on bottom, lb, is inside
      end if
c
      go to 250
c
c --- up/down symmetric part
c
  200 lb  = lt
      alb = alt
c
  250 if (ll .eq. 0 .and. all .le. 0.0  .or.
     .    lr .eq. 0 .and. alr .le. 0.0  .or.
     .    lt .eq. 0 .and. alt .le. 0.0  .or.
     .    lb .eq. 0 .and. alb .le. 0.0)
     .    call STOP ('subroutine FINDINTSC: bad things happening', 117)
      return
c
      end

      subroutine fitdriver (iounit, rplasbdry, zplasbdry, nplasbdry,
     .                      residin, maxinit, iterations, updownsym,
     .                      jsymetric, nw, nh, map, psi1d, coeff,deltar,
     .                      deltarzsq, deltarsq, interior, omega, psi2d,
     .                      rmhdgrid, zmhdgrid, psibdry, residmax,
     .                      psidifbdext, iextiters)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- subroutine is driver for least squares fitting routine.
c   This subroutine repeatedly adjusts a subset of the psi values
c   on the rectangular boundary of the plasma. Points on the rectangular
c   boundary which are not fit are updated using interpolation of the
c   fitted boundary values. The psi values are adjusted until the
c   difference in interpolated values on the plasma surface
c   match as closely as possible the given psi value on the plasma
c   surface.
c
c --- input
c
c  nfit      no of points on rectangular boundary of computational
c             grid which will be fit( for the typical nw by nh grid
c             there are a maximum of 2(nw+nh-2) such points. Typically
c             nfit is taken as a fraction of this maximum to limit
c             the size of the fitting problem.
c
c  mresid     no of residuals . mresid is the no of points on the
c             plasma boundary which are used in the residual calculations
c             (typically there are many more points on the plasma boundary
c             than are used in the minimization process (to limit the size
c             of matrix xjac)
c
c  nw,nh      # grid    points in r and z direction
c  jsymetric  set to 1 for non up down symmetric case,
c             set to nh/2+1,(nh must be odd) for up/down
c             symmetric case
c  iounit     unit number for diagnostic output. set to 0 if diagnostics
c             are not to be printed.
c
c  updownsym  logical variable,set to true if this is an up/down
c             symmetric case,false otherwise.
c
c  psi1d(j)   j = 1,2...nw*nh. One-dimensional version of psi(i,j)
c             NOTE THAT DUE TO HISTORICAL REASONS PSI1D IS NOT
c             THE FORTRAN DEFAULT STORAGE SCHEME. INSTEAD,
c             PSI1D(k) = PSI(I,J) FOR K = (I-1)*NH+J !!!!!!!!!!!
c  splinefitbd
c
c ------------------------------------------------------------------ HSJ
c
      include 'fitparms.i'
c
      integer     isetbdry(mplasmax),isetrect(nfitmax),
     .            jsymetric,k1,k2,k3,k4,k5,nplasbdry,
     .            k1min,k2min,k3min,k4min,k5min,nsig,maxfn,iopt,
     .            k,kk,kkk,ierr,ier,nfit,iounit,nw,nh,ksave(2),
     .            icross,infer,iterations,lastpt,iextiters,
     .            k1save,k2save,k3save,k4save,k5save,map(*),
     .            maxinit,mresid,nbdrypts,nplasbdrysym,j,i,
     .            msave,nsave
c
      logical     updownsym, interior
      real*8      xjac(mplasmax,nfitmax),xjtj(nsymfit),
     .            residual(mplasmax),parm(4),xfit(nfitmax),
     .            work(nworkzx),psi1d(*),rjsymetric,psibdry,
     .            prod,deltar,deltarsq,deltarzsq,eps,hold,omega,
     .            psi2d(nw,nh),coeff(9,*),ssq,zplasbdry(*),
     .            residmax,rincr,rmhdgrid(*),residin,delta,
     .            rplasbdry(*),zmhdgrid(*),psidifbdext(*),
     .            fitresmax,xjacmax
c
c      integer     kstore, kstorest          ! for INCLUDE file storage.i
c      real*8      sdum, tdum, udum, vdum, wdum, xdum, ydum, zdum, rdum
c
      include 'storage.i'
c
      equivalence (work(1)    , sdum(1))
      equivalence (xjtj(1)    , tdum(1))
      equivalence (residual(1), udum(1))
      equivalence (isetbdry(1), udum(mplasmax+1))
      equivalence (isetrect(1), udum(mplasmax+1+mplasmax+1))
c
c --- note that subroutine RESIDUALFIT also uses some storage in storage.i
c
      external residualfit           ! called by least squares fitter
c
c --- first some error checking
c
      mresid = mresidual
      nfit   = nfitpoints
      ierr   = 0
      if (nsymfit .gt. kstorest .or. nworkzx .gt. kstorest .or.
     .      nfitmax+mplasmax .gt. kstore  .or.
     .      2*mplasmax+nfitmax+2 .gt. kstore) then
          ierr = 1
          if (iounit .ne. 0)  write (iounit, 1) kstore,nsymfit,nfitmax,
     .                                          mplasmax,nsymfit,nworkzx
    1     format ('  subroutine FITDRIVER reports:'   /
     .            '  incompatible parameter settings' /
     .            '  kstore =',i5,' nsymfit = ', i5,
     .            '   nfitmax =',i5                   /
     .            '   mplasmax =',i5,' nsymfit =',i5,' nworkzx =',i5)
      end if
      if (nfit .gt. nfitmax) then
          ierr = 1
          if (iounit .ne. 0)  write (iounit, 2)  nfit, nfitmax
    2     format ('  subroutine FITDRIVER reports:' /
     .            '  nfit must be .le. nfitmax '    /
     .            '  nfit,nfitmax =',2(2x,i7)       /
     .            '  program must stop')
      end if
      if (mresid .gt. mplasmax) then
          ierr = 1
          if (iounit .ne. 0)  write (iounit, 3)  mresid, mplasmax
    3     format ('  subroutine FITDRIVER reports:' /
     .            '  mresid must be .le. mplasmax ' /
     .            '  mresid,mplasmax =',2(2x,i7)    /
     .            '  program must stop')
      end if
      if (mresid .le. nfit) then
          ierr = 1
          if (iounit .ne. 0)  write (iounit, 4)  mresid, nfit
    4     format ('  subroutine FITDRIVER reports:'          /
     .            '  problem is not overdetermined'          /
     .            '  # residuals =',i5, '  # fit parms =',i5 /
     .            '  program must stop')
      end if
c
c --- check symmetric or non-symmetric input case
c
      if (jsymetric .eq. 1 .and. .not. updownsym) then
          nbdrypts = 2 * (nw+nh-2)    ! # pts on rectangular boundary
          if (nfit .lt. 6) then
              if (iounit .gt. 0)  write (iounit, 6) nfit
    6         format ('  subroutine FITDRIVER reports:'        /
     .                '  numbr of fitted boundary points must' /
     .                '  be ge to 6. # points selected by user =',2x,i5)
              ierr = 1
          end if
      else if (jsymetric .eq. nh/2 + 1 .and. updownsym) then
          nbdrypts = nh+nw-1 ! # pts on rect. bdry with up/down symmetry
          if (nfit .lt. 4 .and. .not. splinefitbd) then
              if (iounit .gt. 0)  write (iounit, 7) nfit
    7         format ('  subroutine FITDRIVER reports:'         /
     .                '  number of fitted boundary points must' /
     .                '  be ge to 4. # points selected by user =',2x,i5)
              ierr = 1
          end if
          if (nfit .lt. 7 .and. splinefitbd) then
              if (iounit .gt. 0)  write (iounit, 20) nfit
   20         format ('  subroutine FITDRIVER reports:'         /
     .                '  number of fitted boundary points must' /
     .                '  be ge to 7. # points selected by user =',2x,i5)
              ierr = 1
          end if
      else
          if (iounit .gt. 0)
     .    write  (iounit, 8)  jsymetric, nw, nh, updownsym
    8     format ('  subroutine FITDRIVER reports:'  /
     .            '  jsymetric is not set correctly' /
     .            '  jsymetric,nw,nh =',3(2x,i7)     /
     .            '  updownsym =',l8)
          ierr = 1
      end if
      if (mresid .gt. nplasbdry) then
          ierr = 1
          if (iounit .ne. 0)  write (iounit, 5)  mresid, nplasbdry
    5     format ('  subroutine FITDRIVER reports:'            /
     .            '  # requested residuals is greater than '   /
     .            '  # available ,mresid,nplasbdry =',2(2x,i5) /
     .            '  program must stop')
      end if
      if (ierr .ne. 0)
     .  call STOP ('subroutine FITDRIVER: problem #1', 166)
c
c --- next set up for call to least squares fitter
c
      nsig    = 3         ! number of digits for convergence
      eps     = 1.0e-6    ! consider converged if rel diff in ssq < this
      delta   = 1.0e-6    ! consider converged if norm grad (ssq) < this
      maxfn   = 10*nfit   ! max function evaluations, 10 * # unknowns
      iopt    = 2         ! user supplies initial values of parm..
c                         ..and strict descent is desired
      parm(1) = 1.0       ! initial value of marq. parameter
      parm(2) = 2.0       ! scaling factor for parm(1)
c
      parm(3) = 220.0     ! upper bound on parm(1)
      parm(4) = 0.01      ! use central diffs if parm(1) less than this
c
c --- set up initial xfit vector,and index vector to rectangular
c --- boundary, isetrect
c
      k          = 0
      rjsymetric = FLOAT (jsymetric)
      if (updownsym) then
          xfit(1) = psi1d(jsymetric)! first point for symmetric boundary
          xfit(nfit) = psi1d((nw-1)*nh+nh/2+1) ! last point for sym bdry
          rincr      = FLOAT (nbdrypts) / FLOAT ((nfit-1))
          k          = 1
          isetrect(k   ) = jsymetric
          isetrect(nfit) = (nw-1)*nh+nh/2+1
          do while (k .le. nfit-2)
              k  = k + 1
              kk = rjsymetric + (k-1)*rincr
              if (kk .le. nh) then             ! up the lhs
                   xfit(k)     = psi1d(kk)
                   isetrect(k) = kk
              else if (kk .le. nh+nw-1) then   ! across the top
                   kkk = kk-nh
                   kkk = (kkk+1)*nh
                   xfit(k) = psi1d(kkk)
                   isetrect(k) = kkk
              else if (kk .lt. nh+nw-1+nh/2) then  ! down rhs
                   kkk = kk-nh-nw+1
                   kkk = nh-kkk
                   kkk = (nw-1)*nh+kkk
                   xfit(k) = psi1d(kkk)
                   isetrect(k) = kkk
              else
                   if (iounit .ne. 0)
     .             write  (iounit, 9)  kk, nw, nh, nfit
    9              format ('  subroutine FITDRIVER reports:'     /
     .                     '  error in symmetric boundary setup' /
     .                     '  kk,nw,nh,nfit = ',4(2x,i5))
                   call STOP ('subroutine FITDRIVER: sym bdry', 127)
              end if
          end do
            if (splinefitbd) then
c
c             include the two corner points by changing the
c             closest points to the corner to the corner points:
c
              k1min = nw*nh
              k2min = nw*nh
              do k=2,nfit-1
                  k1    = IABS (isetrect(k)-nh)        ! upper lh corner
                  k2    = IABS (isetrect(k)-nw*nh)     ! upper rh corner
                  k1min = MIN (k1, k1min)
                  k2min = MIN (k2, k2min)
                  if (k1min .eq. k1)  k1save = k
                  if (k2min .eq. k2)  k2save = k
              end do
              if (k1min .ne. 0) then
                  xfit(k1save) = psi1d(nh)
                  isetrect(k1save) = nh
              end if
              if (k2min .ne. 0) then
                  xfit(k2save) = psi1d(nw*nh)
                  isetrect(k2save) = nw*nh
              end if
            else  ! if not spline fit then do not include the two corner
c                   points. the reason for this is that these corner points
c                   do not enter into the five point finite difference form and
c                   their value is irrelevant and the fit would be singular
c                   if we tried to include them
c
                do k=2,nfit-1
                   if (isetrect(k) .eq. nh) then
                        if (isetrect(k-1) .lt. nh-1 ) then
                               isetrect(k) = nh-1
                        else if (isetrect(k+1) .gt. 2*nh) then
                               isetrect(k) = 2*nh
                        else      ! just delete the point
                            do j=k+1,nfit
                                isetrect(j-1) = isetrect(j)
                            end do
                            nfit = nfit-1
                        end if
                   end if
                   if (isetrect(k) .eq. nw*nh) then
                        if (isetrect(k-1) .lt. (nw-1)*nh) then
                               isetrect(k) = (nw-1)*nh
                        else if (isetrect(k+1) .lt. nw*nh-1) then
                               isetrect(k) = nw*nh-1
                        else      ! just delete the point
                            do j=k+1,nfit
                                isetrect(j-1) = isetrect(j)
                            end do
                            nfit = nfit-1
                        end if
                   end if
                end do
            end if
      else                                          ! non-symmetric case
          rincr  = FLOAT ( nbdrypts) / FLOAT (nfit-1)
          k = 1
          xfit(k) = psi1d(1)
          isetrect(k) = 1
          do while (k .lt. nfit)
              k = k+1
              kk = rjsymetric+(k-1)*rincr
              if (kk .le. nh) then                   ! up the rhs
                   xfit(k) = psi1d(kk)
                   isetrect(k) = kk
              else if (kk .le. nh+nw-1) then         ! across the top
                   kkk = kk-nh
                   kkk = (kkk+1)*nh
                   xfit(k) = psi1d(kkk)
                   isetrect(k) = kkk
              else if (kk .le. 2*nh+nw-2) then       ! down rhs
                   kkk = kk-nh-nw+1
                   kkk = nh-kkk
                   kkk = (nw-1)*nh+kkk
                   xfit(k) = psi1d(kkk)
                   isetrect(k) = kkk
              else if (kk .le. 2*(nw+nh-2)) then     ! across the bottom
                   kkk = kk-2*nh-nw+2
                   kkk = nw-kkk
                   kkk = (kkk-1)*nh+1
                   xfit(k) = psi1d(kkk)
                   isetrect(k) = k
              else
                   if (iounit .ne. 0)
     .             write  (iounit, 10)  kk, nw, nh, nfit
   10              format ('  subroutine FITDRIVER reports:'         /
     .                     '  error in non-symmetric boundary setup' /
     .                     '  kk,nw,nh,nfit = ',4(2x,i5))
                   call STOP ('subroutine FITDRIVER: non sym', 128)
              end if
          end do
c
c             include the other three corner points and the two
c             points half way up the vertical sides by changing the
c             closest points to these points:
c
              k1min = nw*nh
              k2min = nw*nh
              k3min = nw*nh
              k4min = nw*nh
              k5min = nw*nh
              do k=2,nfit-1
                  k1 = IABS (isetrect(k)-nh/2-1)       ! middle lh side
                  k2 = IABS (isetrect(k)-nh)           ! upper lh corner
                  k3 = IABS (isetrect(k)-nw*nh)        ! upper rh corner
                  k4 = IABS (isetrect(k)-(nw-1)*nh-nh/2-1) ! mid rh side
                  k5 = IABS (isetrect(k)-(nw-1)*nh-1)  ! lower rh corner
                  k1min = MIN (k1,k1min)
                  k2min = MIN (k2,k2min)
                  k3min = MIN (k3,k3min)
                  k4min = MIN (k4,k4min)
                  k5min = MIN (k5,k5min)
                  if (k1min .eq. k1)  k1save = k
                  if (k2min .eq. k2)  k2save = k
                  if (k3min .eq. k3)  k3save = k
                  if (k4min .eq. k4)  k4save = k
                  if (k5min .eq. k5)  k5save = k
              end do
              if (k1min .ne. 0) then
                  xfit(k1save) = psi1d(nh/2+1)         ! middle lh side
                  isetrect(k1save) = nh/2+1
              end if
              if (k2min .ne. 0) then
                  xfit(k2save) = psi1d(nh)             ! upper lh corner
                  isetrect(k2save) = nh
              end if
              if (k3min .ne. 0) then
                  xfit(k3save) = psi1d(nh*nw)          ! upper rh corner
                  isetrect(k3save) = nh*nw
              end if
              if (k4min .ne. 0) then
                  xfit(k4save) = psi1d((nw-1)*nh+nh/2+1)! middle rh side
                  isetrect(k4save) = (nw-1)*nh+nh/2+1
              end if
              if (k5min .ne. 0) then
                  xfit(k5save) = psi1d((nw-1)*nh+1)    ! lower rh corner
                  isetrect(k5save) = (nw-1)*nh+1
              end if
      end if
c
c --- end of xfit setup
c
c --- set up index vector for plasma boundary, isetbdry:
c
      if (updownsym) then              ! symmetry about z = 0 is assumed
c
c        first get the start and end plasma contour points in the
c        upper half of the plasma (it is assumed that exactly
c        two points with z = 0 exist on the plasma boundary
c
         k = 0
         icross = 0
         do while (k .le. nplasbdry-1 .and. icross .lt. 2)
              k = k+1
              prod = zplasbdry(k)*zplasbdry(k+1)
              if (prod .le. 0.0  ) then
                   icross = icross+1
                   if (zplasbdry(k) .ge. 0.0) then
                       ksave(icross) = k
                   else
                       ksave(icross) = k+1
                   end if
               end if
         end do
         if (icross .ne. 2) then
             if (iounit .gt. 0)  write (iounit, 11)
   11        format ('  subroutine FITDRIVER reports:'           /
     .               '  unable to determine symmetric start/end' /
     .               '  of plasma boundary')
             call STOP ('subroutine FITDRIVER: problem #2', 129)
         end if
         if (ksave(2) .lt. ksave(1)) then
             hold = ksave(2)
             ksave(2) = ksave(1)
             ksave(1) = hold
         end if
         nplasbdrysym = ksave(2)-ksave(1)+1
         do j=ksave(1),ksave(2),1
              if (zplasbdry(j) .lt. 0.0) then
                  if (iounit .ne. 0)
     .            write  (iounit, 12)  ksave(1), ksave(2)
   12             format (' subroutine FITDRIVER reports:'          /
     .                    '   symmetric plasma boundary points not' /
     .                    '   properly ordered, ksave = ', 2(2x, i5))
                  call STOP ('subroutine FITDRIVER: ksave problem', 130)
               end if
         end do
c
         isetbdry(1) = ksave(1)
         isetbdry(mresid) = ksave(2)
         rincr = FLOAT (nplasbdrysym) / FLOAT (mresid-1)
         k = 1
         do while (k .lt. mresid -1)
            k = k+1
            kk = (k-1)*rincr
            if (kk .ge. isetbdry(mresid)) then
               if (iounit .ne. 0)  write (iounit, 13)  kk, mresid
   13          format ('  subroutine FITDRIVER reports:' /
     .                 '  kk,mresid =',2(2x,i5)          /
     .                 '  program must stop')
               call STOP ('subroutine FITDRIVER: sm pl bd', 131)
            end if
            isetbdry(k) = kk+ksave(1)
         end do
c
c        check if any boundary points are duplicates
c
         kkk = 0
         do k=1,mresid-1
             if (isetbdry(k) .eq. isetbdry(k+1)) then ! delete duplicate
                 kkk = kkk+1
                 do kk=k+1,mresid-1
                     isetbdry(kk) = isetbdry(kk+1)
                 end do
             end if
         end do
         mresid = mresid-kkk
c
c        check if any boundary points have negative z values
c
         kkk = 0
         do k=1,mresid-1
            if (zplasbdry(isetbdry(k)) .lt. 0.0) then ! delete duplicate
                kkk = kkk+1
                do kk=k+1,mresid-1
                    isetbdry(kk) = isetbdry(kk+1)
                end do
            end if
         end do
         mresid = mresid-kkk
      else                    ! non-symmetric case
          rincr  = FLOAT (nplasbdry) / FLOAT (mresid-1)
          lastpt = nplasbdry-rincr
          rincr  = FLOAT (lastpt) / FLOAT (mresid-2)
          k      = 1
          isetbdry(k) = 1
          do while (k .le. mresid)
              k = k+1
              kk = (k-1)*rincr
              if (kk .gt. lastpt -1) then
                   if (iounit .ne. 0)
     .             write  (iounit, 14)  kk, lastpt, nplasbdry
   14              format ('  subroutine FITDRIVER reports:'  /
     .                     '  kk,lastpt,nplasbdry =',3(2x,i5) /
     .                     '  program must stop')
                    call STOP ('subroutine FITDRIVER: nsm pl bd', 132)
              end if
              isetbdry(k) = kk
           end do
      end if
c
c --- call the fitter. The fitter calls subroutine residualfit
c --- which adjusts the boundary values of psi in fit vector xfit.
c --- xfit is copied into psi1d and the exterior boundary value problem
c --- delstar psi =0  is solved in the exterior region with boundary
c --- conditons given by psi1d on the rectangular boundry and with
c --- psi = const value on the plasma boundary.
c
      write (6, *)  'final values of mresid,nfit used =',mresid,nfit
      call zxssq1 (residualfit,mresid,nfit,nsig,eps,delta,maxfn,
     .             iopt,parm,xfit,ssq,residual,xjac,mplasmax,xjtj,work,
     .             infer,ier,zmhdgrid,psidifbdext,iextiters,
     .             residin,maxinit,iterations,updownsym,
     .             jsymetric,nw,nh,map,psi1d,coeff,deltar,
     .             deltarzsq,deltarsq,interior,omega,psi2d,
     .             rmhdgrid,psibdry,residmax,
     .             isetbdry,isetrect,rplasbdry,zplasbdry)
      if (iounit .gt. 0) then
          fitresmax = 0.0
          do j=1,mresid
              fitresmax = MAX (fitresmax, ABS (residual(j)))
              if (fitresmax .eq. ABS (residual(j)))  i = j
          end do
          xjacmax = 0.0
          do k=1,mresid
             do kk=1,nfit
                  xjacmax = MAX (ABS (xjac(k,kk)), xjacmax)
                  if (xjacmax .eq. ABS (xjac(k,kk))) then
                       msave = k
                       nsave = kk
                  end if
             end do
          end do
          xjacmax = xjac(msave,nsave)
          write  (iounit, 15)  infer, ier
   15     format ('  subroutine ZXSSQ1 reports:' /
     .            '  infer,ier =',2(2x,i5))
          write  (iounit, 16)  (work(i), i=1,5)
   16     format ('  norm of ssq gradient =',1pe12.2      /
     .            '  # funct evaluations  =',1pe12.2      /
     .            '  est sign digits in xfit  = ',1pe12.2 /
     .            '  marq param final value   = ',1pe12.4 /
     .            '  # iterations             = ',1pe12.2)
          write  (iounit, 17)  fitresmax, i, rplasbdry(i), zplasbdry(i)
   17     format ('  maximum residual = ',1pe12.4         /
     .            '  occurs at boundary point number ',i5 /
     .            '  with r,z location = ',2(2x,1pe12.4))
          write  (iounit, 18)  xjacmax, msave, nsave,
     .                         rplasbdry(msave), zplasbdry(msave),
     .                         isetrect(nsave)
   18     format ('  max d(residual(i))/d(xfit(j) =',1pe12.4     /
     .            '  occurs at i,j =',2(2x,i5)                   /
     .            '  rplasbdry(i),zplasbdry(i), = ,',2(2x,1pe12.4) /
     .            '  index of boundary psi point corresponding ' /
     .            '  to xfit(j) is ', i5)
      end if
c
      return
c
      end

      subroutine interface_dump_psi (psitrace)
c
      USE param
      USE mhdpar
      USE mhdgrid
      USE nub 
      USE io,only : n66,n77 
      USE mhdcom
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c    bring in the psi values from cntour so they can be dumped
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'nub.i'
c      include 'mhdpar.i'
c      include 'mhdcom.i'
c      include 'small.i'
c      include 'mhdgrid.i'
c      include 'zerocom.i'
c
      cconst = -1.0
      call multpl1 (p,nwh,cconst)
c
c     use this call to generate the bicubic spline:
c
****  call ibcccu (p, rmhdgrid, nw, zmhdgrid, nh, cspln,
**** .             nw, wnoperm, ier)
c
      call dump_psi_values (p, nw, nh, rmhdgrid, zmhdgrid,
     .                      n77, 1, zero, 0,
     .                      map, .false., psif, mf, psitrace)
      return
c
      end

      subroutine mapvar (nw, nh, wzero, map, nplas, rp, zp, np,
     .                   rmhdgrid, zmhdgrid, coeff, updownsym)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- THIS MODULE DOES NOT USE ANY INCLUDE FILES AND DOES NOT CALL
c --- ROUTINES OUTSIDE OF THIS MODULE, SO IT IS SELF-CONTAINED.
c --- ITS PRIMARY PURPOSE IS TO CALCULATE THE COEFFICIENTS OF THE
c --- VARIABLE FINITE DIFFERENCE FROM THE GRAD-SHAFRANOV EQUATION
c ------------------------------------------------------HSJ-----------
c --- subroutine creates the vector map(k),k = 1,2,...nw*nh ,where the grid
c --- points are numbered by columns as k  =  (i-1)*nh+j,i=1,2..nw,j=1,2...nh
c --- map(k)is set to zero if grid point k is outside or on the boundary
c --- of the plasma. map(k) is set to j for points interior to the plasma.
c --- Where the unknowns, j, are inside the plasma and are numbered along
c --- rows, starting at the bottom of the (r,z) grid and working upward
c --- across the rows.
c --- hence map(k) = l means global unknown at grid point k = (i-1)*nh+j
c --- is mapped into local (to the plasma interior) unknown grid point l.
c --- nplas is the total number of points found inside the plasma
c --- boundary given by (rp(j),zp(j)),j = 1,2..np
c
c --- next the subroutine calculates the coefficients of the unknown
c --- psi values. using the 5 point diamond difference we have
c                            psi lt
c
c              psi ll        psi l         psi lr
c
c                            psi lb
c any combination of the four points (ll,lt,lr,lb) can be on the
c plasma boundary.
c  there are 16 cases to be accounted for:
c  ll  lt  lr  lb
c   1   1   1   1       case 15 all points inside
c   1   1   1   0            14 lb outside
c   1   1   0   1            13 lr outside
c   1   1   0   0            12 lr,lb outside
c   1   0   1   1            11 lt outside
c   1   0   1   0            10 lt,lb outside
c   1   0   0   1             9 lt,lr outside
c   1   0   0   0             8 lt,lr,lb outside
c   0   1   1   1             7 ll outside
c   0   1   1   0             6 ll,lb outside
c   0   1   0   1             5 ll,lr outside
c   0   1   0   0             4 ll,lr,lb ouside
c   0   0   1   1             3 ll,lt outside
c   0   0   1   0             2 ll,lt,lb outside
c   0   0   0   1             1 ll,lt,lr outside
c   0   0   0   0             0 ,ll,lt,lr,lb outside
c
c  case 10 has no vertical extent and is pathological
c        5        horizontal
c        0        vertical and no horizontal
c  if any of these three cases occur grid is too coarse. set icoarse
c  and return.
c
c --- input
c
c   nw          number of grid points in radial direction
c   nh                                   axial
c   wzero(k)    k = 1,2...nw*nh . wzero is defined so that
c               wzero(k) = 1. if global point k is inside the
c               plasma and wzero(k) = 0 if k is outside.
c rmhdgrid(i)   i = 1,2..nw horizontal grid
c zmhdgrid(j)   j = 1,2..nh vertical grid
c  rp(i)          i = 1,2..np
c  zp(i)          i = 1,2..np
c  np             (rp,zp) are the points on the plasma boundary
c
c --- output ---
c
c   map(k)      k = 1,2,...nw*nh
c               see explanation above
c   nplas       number of points inside the plasma
c               (and thus the number of unknown psi values
c                to be found)
c    coeff(1,l)        coefficient of psi ll
c         (2,l)                           lt
c         (3,l)                           lr
c         (4,l)                           lb
c         (5,l)                           source term
c         (6,l)
c         (7,l)
c         (8,l)
c         (9,l)
c
c ---------------------------------------------------------- HSJ 3/26/92
c
      logical  icoarse, interior, updownsym, jstart
      integer  k, i, j, nw, nh, ll, lt, lr, lb, l, map(*), ii, np,
     .         nplas, lll, llr, llt, llb, jsymetric
      real*8   all, alt, alr, alb, wzero(*), rmhdgrid(*), zmhdgrid(*),
     .         coeff(9,*), deltar, deltaz, deltarsq, deltarzsq,
     .         rl, zl, rp(*), zp(*)
c
      interior  = .true.
      jsymetric = 1
      if (updownsym)  jsymetric = nh/2+1 ! jsymetric is starting j value
      do j=1,nw*nh
        map(j) = 0
      end do
      l = 0
      do j=jsymetric,nh
        do i=1,nw
          k = (i-1)*nh+j
          if (wzero(k) .gt. 0.0) then
            l      = l + 1
            map(k) = l
          end if
        end do
      end do
      nplas = l
c
      deltar    = rmhdgrid(2)-rmhdgrid(1)
      deltaz    = zmhdgrid(2)-zmhdgrid(1)
      deltarsq  = deltar*deltar
      deltarzsq = deltarsq/(deltaz*deltaz)
      icoarse   = .false.
      do j=jsymetric,nh
          jstart = .false.
          if (updownsym .and. j .eq. jsymetric)  jstart = .true.
          do i=1,nw
              k = (i-1)*nh+j
              if (map(k) .ne. 0) then        ! do only if pt l is inside
                  ll  = map(k-nh)            ! ll is to left of l
                  lt  = map(k+1)             ! lt is above l
                  lr  = map(k+nh)            ! lr is to right of l
                  lb  = map(k-1)             ! lb is below l
                  l   = map(k)               ! l  is the current unknown
                  lll = MIN (ll, 1)
                  llr = MIN (lr, 1)
                  llt = MIN (lt, 1)
                  llb = MIN (lb, 1)
                  ii  = lll*8+llt*4+llr*2+llb        ! ii is case number
                  if (ii .eq. 0 .or. ii .eq. 5 .or. ii .eq. 10) then
                     write (6, *)  'pathological grid encountered'
                     write (6, *)  'nw, nh =', nw, nh
                     write (6, *)  'grid will be refined'
                     icoarse = .true.
                     return
                  end if
                  rl = rmhdgrid(i)
                  zl = zmhdgrid(j)
                  call findintsc (rl,zl,rp,zp,np,deltar,deltaz,
     .                            ll,lt,lr,lb,all,alt,alr,alb,
     .                            updownsym,jstart)
                  call coefcalc (coeff(1,l),all,alt,alr,alb,
     .                           deltar,deltarzsq,deltarsq,rl,interior)
              end if    ! map(k) inside
          end do        ! end on row index i
      end do            ! end on column index j
c
      return
c
      end

      subroutine residualfit (xfit,mresid,nfit,fresid,zmhdgrid,
     .               psidifbdext,icall,residin,maxinit,iterations,
     .               updownsym,jsymetric,nw1,nh1,map,psi1d,coeff,
     .               deltar,deltarzsq,deltarsq,interior,omega,psi2d,
     .               rmhdgrid,psibdry,residmax,
     .               isetbdry,isetrect,rplasbdry,zplasbdry)
c
c
c ----------------------------------------------------------------------
c --- subroutine called by least squares fitter. Arguments after zmhdgrid
c --- in the above calling sequence are passed through zxssq1 without
c --- modfication.
c ----------------------------------------------------------------------
c --- input
c  xfit(j)       j = 1,2..nfit, the current value of the fit parameters
c
c  mresid         the number of residuals . (this is the number of
c                 points on the plasma boundary which are used in
c                      the fit)
c  nfit           the number of unknown fit values
c
c  splinefitbd
c
c --- output
c
c  fresid(j)      j = 1,2...mresid   the residuals
c  psidifbdext(j)    j = 1,2..?  sum squares residual as function
c                    of call no.
c
c --- remaining arguments are explained in subroutine OUTER_SOLUTION
c ----------------------------------------------------------------------
c
c --- arguments passed on to outer_solution:
c
      USE mhdpar
      USE bicube
      USE replace_imsl,                ONLY : my_ibcccu,my_dbcevl1, 
     .                                        my_icsicu,my_icsevu

      implicit none
      integer            maxinit,nw1,nh1,jsymetric,map(*),iterations,
     .                   isetbdry(*),isetrect(*)
      logical            updownsym,interior
      real*8             psi1d(*),residin,omega,psi2d(nw1,nh1),psibdry,
     .                   rmhdgrid(*),coeff(9,*),deltar,deltarzsq,
     .                   deltarsq,residmax,psidifbdext(*)
c
      integer    ndum
      parameter (ndum = 2)   ! ndum is for psi1dr which is equivalenced
c
      include 'fitparms.i'
      include 'imsl.i'

c
c --- additional arguments for this subroutine
c
      integer    j,nfit,mresid,jsave,ier,nknots,npts,i,k,nhh,
     .           icall,ksave
c
      real*8     fresid(*),xfit(*),bpar(4),psiknot(nfitmax),
     .           zknot(nfitmax),rknot(nfitmax),c1dspl(nfitmax,3),
     .           zmhdgrid(*),psi1dr(ndum),rloc,zloc,sumsq,
     .           rplasbdry(*),zplasbdry(*)
c
c --- include file storage.i variables:
c
c      integer    kstore
c      real*8     sdum,tdum,udum,vdum,wdum,xdum,ydum,zdum,rdum
c
      include 'storage.i'   ! note that some of storage is in use by..
c                           ..the driver that calls this routine
c
      equivalence (vdum(1),zknot(1))
      equivalence (vdum(nfitmax+1),rknot(1))
      equivalence (vdum(2*nfitmax+1),psiknot(1))
      equivalence (wdum(1),c1dspl(1,1))
      equivalence (xdum(1),psi1dr(1))
c
      if (3*nfitmax .gt. kstore)
     .  call STOP ('subroutine RESIDUALFIT: kstore problem', 133)
      if (nw1 .ne. nw  .or.  nh1 .ne. nh)
     .  call STOP ('subroutine RESIDUALFIT: nw1 problem', 134)
c
c --- copy xfit into the proper spaces in psi1d:
c
      do j=1,nfit
          psi1d(isetrect(j)) = xfit(j)
      end do
c
      if (splinefitbd) then
c
c --- use spline for points on bdry not fit directly
c
c --- interpolate to get values of psi1d on the boundary in between
c --- the fitted points:
c --- left vertical boundary:
c
      if (updownsym) then
c
c        find index pointing to nh:
c
         j = 1
         do while (isetrect(j) .lt. nh .and. j .le. nfit)
           j = j+1
         end do
         nknots = j ! # knots to use in fit (# of unknowns on this side)
         if (isetrect(nknots) .ne. nh)
     .     call STOP ('subroutine RESIDUALFIT: lv1 problem', 135)
         do k=1,nknots
           psiknot(k) = psi1d(isetrect(k))
           zknot(k) = zmhdgrid(isetrect(k))  ! ok for first side
         end do
c
c        boundary condition on spline. At the starting point we have
c        zero derivative in psi (in the z direction) due to symmetry.
c        at the ending point we use the natural spline condition:
c
         bpar(1) = 1.0
         bpar(2) = (6.0/deltar)*(psiknot(2)-psiknot(1))/deltar
         bpar(3) = 0.0
         bpar(4) = 0.0
c
c        fit a one d spline to psi,as a function of z:
c
         call my_icsicu(zknot,psiknot,nknots,bpar,c1dspl,nfitmax,ier)
         npts = nh/2+1                    ! # values to be determined
         call my_icsevu (zknot,psiknot,nknots,c1dspl,nfitmax,
     .                zmhdgrid(npts),psi1d(npts),npts,ier)
         if (ier .ne. 0)
     .     call STOP ('subroutine RESIDUALFIT: lv2 problem', 136)
c
      else                             ! non symmetric case
c
c         find index pointing to nh:
c
          j = 1
          do while (isetrect(j) .lt. nh .and. j .le. nfit)
            j = j + 1
          end do
          nknots = j
          if (isetrect(nknots) .ne. nh)
     .      call STOP ('subroutine RESIDUALFIT: rnonsym1 problem', 137)
          do k=1,nknots
            psiknot(k) = psi1d(isetrect(k))
            zknot(k) = zmhdgrid(isetrect(k))
          end do
c
c         boundary condition on spline (here we use a natural spline
c         on both ends)
c
          bpar(1) = 0.0
          bpar(2) = 0.0
          bpar(3) = 0.0
          bpar(4) = 0.0
          call my_icsicu (zknot,psiknot,nknots,bpar,c1dspl,nfitmax,ier)
          npts = nh                  ! number of values to be determined
          call my_icsevu (zknot,psiknot,nknots,c1dspl,nfitmax,
     .                 zmhdgrid(1),psi1d(1),npts,ier)
          if (ier .ne. 0)
     .    call STOP ('subroutine RESIDUALFIT: rnonsym2 problem #1', 138)
c
      end if
      jsave = nknots ! jsave is start index for next segment of boundary
c
c --- top horizontal boundary:
c
c         starting index is nknots from above
c         find index pointing to nw*nh
c
          j      = jsave
          do while (isetrect(j) .lt. nw*nh .and. j .le. nfit)
            j    = j + 1
          end do
          ksave  = j
          nknots = ksave-nknots+1
          if (isetrect(j) .ne. nw*nh)
     .      call STOP ('subroutine RESIDUALFIT: th1 problem', 139)
          do k=1,nknots
               j          = isetrect(jsave+k-1)
               psiknot(k) = psi1d(j)
c
c              find i correspondiong to j:
c              nh must divide j exactly
c
               if (MOD (j, nh) .ne. 0)
     .           call STOP ('subroutine RESIDUALFIT: thnh problem', 140)
               i        = j/nh
               rknot(k) = rmhdgrid(i)
          end do
c
c         boundary condition on spline (here we use a natural spline
c         on both ends)
c
          bpar(1) = 0.0
          bpar(2) = 0.0
          bpar(3) = 0.0
          bpar(4) = 0.0
          call my_icsicu (rknot,psiknot,nknots,bpar,c1dspl,nfitmax,ier)
          npts    = nw                       ! # values to be determined
          call my_icsevu (rknot,psiknot,nknots,c1dspl,nfitmax,
     .                 rmhdgrid,psi1dr,npts,ier)
          do i=1,nw
            j        = i*nh
            psi1d(j) = psi1dr(i)
          end do
          if (ier .ne. 0)
     .      call STOP ('subroutine RESIDUALFIT: th2 problem', 141)
          jsave = ksave
c
c --- right vertical boundary. the spline routines require that
c --- zknots increases so we have some additional manipulations
c
      if (updownsym) then
c
c         index pointing to (nw-1)*nh+nh/2+1 is nfit
c         (constructed that way in subroutine FITDRIVER). Note that
c         the number of knots may not work out to be the same on
c         the left and right vertical boundaries due to integer
c         truncation in the math
c
          nknots = nfit-jsave+1
          i      = 0
          do k=nfit,jsave,-1
               i          = i+1
               psiknot(i) = psi1d(isetrect(k))
               zknot(i)   = zmhdgrid(isetrect(k)-(nw-1)*nh)
          end do
          if (zknot(nknots) .ne. zmhdgrid(nh))
     .      call STOP ('subroutine RESIDUALFIT: rvm1 problem', 142)
c
c         boundary condition on spline. At the starting point we have
c         zero derivative in psi (in the z direction) due to symmetry.
c         at the ending point we use the natural spline condition:
c
          bpar(1) = 1.0
          bpar(2) = (6.0/deltar)*(psiknot(2)-psiknot(1))/deltar
          bpar(3) = 0.0
          bpar(4) = 0.0
          call my_icsicu(zknot,psiknot,nknots,bpar,c1dspl,nfitmax,ier)
          npts    = nh / 2 + 1               ! # values to be determined
          call my_icsevu (zknot,psiknot,nknots,c1dspl,nfitmax,
     .                 zmhdgrid(npts),psi1dr,npts,ier)
          if (ier .ne. 0)
     .      call STOP ('subroutine RESIDUALFIT: rv2 problem', 143)
          do i=npts,nh
            j        = (nw-1)*nh+i
            psi1d(j) = psi1dr(i-npts+1)
          end do
c
      else             ! non symmetric case
c
c         find index pointing to nh+nw-1 and to 2*nh+nw-2
c
          j = jsave
          k = 0
          do while (j .lt. nfit .and. k .ne. j)
            j = j+1
            if (isetrect(j) .eq. (nw-1)*nh+1)  k = j
          end do
          nknots = k-jsave+1
          jsave  = k
          do k=1,nknots
            psiknot(k) = psi1d(isetrect(k))
            zknot(k) = zmhdgrid(isetrect(k)-(nw-1)*nh)
          end do
          if (psiknot(1     ) .ne. psi1d((nw-1)*nh+1) .or.
     .        psiknot(nknots) .ne. psi1d( nw*nh     ))
     .        call STOP ('subroutine RESIDUALFIT: rvs2 problem', 144)
          if (zknot(1     ) .ne. zmhdgrid(1 ) .or.
     .        zknot(nknots) .ne. zmhdgrid(nh))
     .        call STOP ('subroutine RESIDUALFIT: rvs3 problem', 145)
c
c         boundary condition on spline (here we use a natural spline
c         on both ends)
c
          bpar(1) = 0.0
          bpar(2) = 0.0
          bpar(3) = 0.0
          bpar(4) = 0.0
          call my_icsicu (zknot,psiknot,nknots,bpar,c1dspl,nfitmax,ier)
          npts = nh                  ! number of values to be determined
          call my_icsevu (zknot,psiknot,nknots,c1dspl,nfitmax,
     .                 zmhdgrid,psi1dr,npts,ier)
          if (ier .ne. 0)
     .    call STOP ('subroutine RESIDUALFIT: rnonsym2 problem #2', 146)
          do j=1,nh
            k        = (nw - 1) * nh + j
            psi1d(k) = psi1dr(j)
          end do
c
      end if
c
c --- across the bottom (non-symmetric case only)
c --- we need to set up the points in increasing r for spline routines
c
      if (.not. updownsym) then
            j = 1
            rknot(j)   = rmhdgrid(1)
            psiknot(1) = psi1d(1)
            do  k=nfit,jsave,-1
                i          = (isetrect(k)-1)/nh +1
                j          = j+1
                rknot(j)   = rmhdgrid(i)
                psiknot(j) = psi1d(isetrect(k))
            end do
            nknots = nfit-jsave+1
            if (  rknot(nknots) .ne. rmhdgrid(nw) .or.
     .          psiknot(nknots) .ne. psi1d((nw-1)*nh+1))
     .            call STOP ('subroutine RESIDUALFIT: lhb problem', 147)
c
c         boundary condition on spline
c         here we use natural spline on both ends
c
          bpar(1) = 0.0
          bpar(2) = 0.0
          bpar(3) = 0.0
          bpar(4) = 0.0
          call my_icsicu(rknot,psiknot,nknots,bpar,c1dspl,nfitmax,ier)
          npts = nw                        ! # values to be determined
          call my_icsevu (rknot,psiknot,nknots,c1dspl,nfitmax,
     .                 rmhdgrid,psi1dr,npts,ier)
          if (ier .ne. 0)
     .      call STOP ('subroutine RESIDUALFIT: rnonsym4 problem', 148)
          do i=1,nw
              k = (i-1)*nh+1
              psi1d(k) = psi1dr(i)
          end do
      end if
      end if    ! end spline fit boundary
c
c --- we are done with loading the boundary values. now
c --- calculate the exterior solution wioth the given values of
c --- psi on the rectangular boundary stored in psi1d:
c
      residin = 1.0e-11   ! TEMPORARY
c
      call outer_solution (residin,maxinit,iterations,
     .                     updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .                     deltar,deltarzsq,deltarsq,interior,omega,pds,
     .                     rmhdgrid,psibdry,residmax)
      write (6, *)  'residin ,residmax =',residin,residmax
      write (6, *)  'maxinit,iterations =',maxinit,iterations
c
c --- we are done with the inner iterations of the exterior solution.
c --- we now have psi over the entire rectangular grid.
c --- (or over the upper half if updownsym  = .true.)
c
       do j=jsymetric,nh
         do i=1,nw
           k          = (i - 1) * nh + j
           psi2d(i,j) = psi1d(k)
         end do
       end do
c
c --- fit a bicubic spline over the upper half (or entire if not
c     symmetric) rectangular grid.
c     use this two-dimensional spline fit to calculate psi on the plasma
c     boundary and obtain the residulas:
c
      nhh = nh
      if (jsymetric .ne. 1)  nhh = jsymetric
      imslmd ='1494c203'
      call my_ibcccu (psi2d(1,jsymetric),rmhdgrid,nw,
     .      zmhdgrid(jsymetric),nhh,cspln,nw,wnoperm,ier)
      sumsq = 0.0
      do j=1,mresid
        rloc = rplasbdry(isetbdry(j))
        zloc = zplasbdry(isetbdry(j))
        call my_dbcevl1(rmhdgrid,nw,zmhdgrid(jsymetric),nhh,
     .               cspln,nw,rloc,zloc,pds,ier,1)
        if (ier .ne. 0)
     .      call STOP ('subroutine RESIDUALFIT: resdb1 problem', 149)
        fresid(j) = psibdry-pds(1)
        sumsq     = sumsq+fresid(j)**2
      end do
      icall = icall+1
      write (6, *)  'iter sumsq (residualfit)',icall,sumsq
      if (icall .gt. maxinit)  icall = maxinit-100   ! save last 100
      psidifbdext(icall) = sumsq
      return
c
      end

      subroutine surfarea (rp, zp, mp, sarea)
c
      USE constnts
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     subroutine calculates the surface area of the flux contour,
c     given its the contour points in one plane.
c
c --- input
c     rp(i)     i = 1...mp  the major radius of point i
c     zp(i)                         elevation
c     mp        number of points
c
c --- output
c     sarea     the surface area in units of rp*zp
c
c --------------------------------------------------------- 11/30/94 HSJ
c
c      include 'constnts.i'    ! pick up pi here
c
      dimension  rp(*), zp(*)
c
      sarea    = 0.0
      do j=1,mp-1
        dl    = SQRT ((rp(j+1)-rp(j))**2+(zp(j+1)-zp(j))**2)
        rj    = rp(j+1)+rp(j)            ! take care of 0.5 below
        sarea = sarea+rj*dl
      end do
      sarea   = sarea*pi                 ! includes 0.5 from above
      return
c
      end

      subroutine torflux2 (psi, nw, nh, rmhdgrid, zmhdgrid, wzero,
     .                     rbd, zbd, nbd, torflux, fpsi, psival,
     .                     npsi, isetw)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- calculate the toroidal flux without using any flux surface
c --- average quantitites (because they become ill-defined near a
c --- separatrix).
c
c     INPUT:
c     argument list
c     rbd(j)            j = 1,2,...nbd
c     zbd(j)
c     nbd               rbd,zbd,nbd decribes the contour
c     fpsi(j)           f(psi)
c     psival(j)         psi grid
c     isetw             if 0 set up wzero
c     psi(i,j)          i = 1,2..nw,j=1,2...nh  psi values
c     rmhdgrid(i)       i = 1,2..nw
c     zmhdgrid(j)       j = 1,2..nh   mhdgrid
c     nw,nh             size of MHD grid
c
c     OUTPUT:
c     wzero(k)          k = 1,2,....nwh
c     torflux           toroidal flux inside surface described by
c                       rbd,zbd
c --------------------------------------------------------- 11/31/94 HSJ
c
      integer npsi, i, j, k, iflag, nw, nh, isetw, nbd
      real*8  torflux, da, rbd(nbd), zbd(nbd), fp, fpsi(*), psival(*),
     .        psi(nw,nh), wzero(*), rmhdgrid(*), zmhdgrid(*)
c
c --- set up wzero - it may not have been set previously
c
      iflag = 2
      if (isetw .eq. 0)
     .  call zlim (wzero, nw, nh, nbd, rbd, zbd,
     .             rmhdgrid, zmhdgrid, iflag)
c
c --- get the toroidal flux inside the given surface
c
      do j=1,nh
        do i=1,nw
          k = (i-1)*nh+j
          if (wzero(k) .ne. 0.0 ) then                    ! get f(psi)
            call interp(psi(i,j),psival,npsi,fpsi,fp)
            torflux = torflux + fp / rmhdgrid(i)
          end if
        end do
      end do
c
      da      = (rmhdgrid(2)-rmhdgrid(1))*(zmhdgrid(2)-zmhdgrid(1))
      torflux = torflux*da
      return
c
      end

      subroutine zxssq1 (func,m,n,nsig,eps,delta,maxfn,iopt,parm,
     .                   x,ssq,f,xjac,ixjac,xjtj,work,infer,ier,
     .                   zmhdgrid,psidifbdext,iextiters,
     .                   residin,maxinit,iterations,updownsym,
     .                   jsymetric,nw,nh,map,psi1d,coeff,deltar,
     .                   deltarzsq,deltarsq,interior,omega,psi2d,
     .                   rmhdgrid,psibdry,residmax,
     .                   isetbdry,isetrect,rplasbdry,zplasbdry)

c
c ----------------------------------------------------------------------
c   IMSL routine name   - zxssq
c ----------------------------------------------------------------------
c
c   computer            - cray/single
c
c   latest revision     - june 1, 1982
c
c   purpose             - minimum of the sum of squares of m functions
c                           in n variables using a finite difference
c                           levenberg-marquardt algorithm
c
c   usage               - call zxssq(func,m,n,nsig,eps,delta,maxfn,iopt,
c                           parm,x,ssq,f,xjac,ixjac,xjtj,work,infer,ier)
c
c   arguments    func   - a user supplied subroutine which calculates
c                           the residual vector f(1),f(2),...,f(m) for
c                           given parameter values x(1),x(2),...,x(n).
c                           the calling sequence has the following form
c                             call func(x,m,n,f)
c                             where x is a vector of length n and f is
c                               a vector of length m.
c                             func must appear in an external statement
c                               in the calling program.
c                             func must not alter the values of
c                               x(i),i = 1,...,n, m, or n.
c                m      - the number of residuals or observations
c                           (input)
c                n      - the number of unknown parameters (input).
c                nsig   - first convergence criterion. (input)
c                           convergence condition satisfied if on two
c                           successive iterations, the parameter
c                           estimates agree, component by component,
c                           to nsig digits.
c                eps    - second convergence criterion. (input)
c                           convergence condition satisfied if, on two
c                           successive iterations the residual sum
c                           of squares estimates have relative
c                           difference less than or equal to eps. eps
c                           may be set to zero.
c                delta  - third convergence criterion. (input)
c                           convergence condition satisfied if the
c                           (euclidean) norm of the approximate
c                           gradient is less than or equal to delta.
c                           delta may be set to zero.
c                             note, the iteration is terminated, and
c                             convergence is considered achieved, if
c                             any one of the three conditions is
c                             satisfied.
c                maxfn  - input maximum number of function evaluations
c                           (i.e., calls to subroutine func) allowed.
c                           the actual number of calls to func may
c                           exceed maxfn slightly.
c                iopt   - input options parameter.
c                         iopt = 0 implies brown's algorithm without
c                           strict descent is desired.
c                         iopt = 1 implies strict descent and default
c                           values for input vector parm are desired.
c                         iopt = 2 implies strict descent is desired with
c                           user parameter choices in input vector parm.
c                parm   - input vector of length 4 used only for
c                         iopt equal two.  parm(i) contains, when
c                           i = 1, the initial value of the marquardt
c                             parameter used to scale the diagonal of
c                             the approximate hessian matrix, xjtj,
c                             by the factor (1.0 + parm(1)).  a small
c                             value gives a newton step, while a large
c                             value gives a steepest descent step.
c                             the default value for parm(1) is 0.01.
c                           i = 2, the scaling factor used to modify the
c                             marquardt parameter, which is decreased
c                             by parm(2) after an immediately successful
c                             descent direction, and increased by the
c                             square of parm(2) if not.  parm(2) must
c                             be greater than one, and two is default.
c                           i = 3, an upper bound for increasing the
c                             marquardt parameter.  the search for a
c                             descent point is abandoned if parm(3) is
c                             exceeded.  parm(3) greater than 100.0 is
c                             recommended.  default is 120.0.
c                           i = 4, value for indicating when central
c                             rather than forward differencing is to be
c                             used for calculating the jacobian.  the
c                             switch is made when the norm of the
c                             gradient of the sum of squares function
c                             becomes smaller than parm(4).  central
c                             differencing is good in the vicinity
c                             of the solution, so parm(4) should be
c                             small.  the default value is 0.10.
c                x      - vector of length n containing parameter
c                           values.
c                         on input, x should contain the initial
c                           estimate of the location of the minimum.
c                         on output, x contains the final estimate
c                           of the location of the minimum.
c                ssq    - output scalar which is set to the residual
c                           sums of squares, f(1)**2+...+f(m)**2, for
c                           the final parameter estimates.
c                f      - output vector of length m containing the
c                           residuals for the final parameter estimates.
c                xjac   - output m by n matrix containing the
c                           approximate jacobian at the output vector x.
c                ixjac  - input row dimension of matrix xjac exactly
c                           as specified in the dimension statement
c                           in the calling program.
c                xjtj   - output vector of length (n+1)*n/2 containing
c                           the n by n matrix (xjac-transposed) * (xjac)
c                           in symmetric storage mode.
c                work   - work vector of length 5*n + 2*m + (n+1)*n/2.
c                         on output, work(i) contains for
c                           i = 1, the norm of the gradient described
c                             under input parameters delta and parm(4).
c                           i = 2, the number of function evaluations
c                             required during the work(5) iterations.
c                           i = 3, the estimated number of significant
c                             digits in output vector x.
c                           i = 4, the final value of the marquardt
c                             scaling parameter described under parm(1).
c                           i = 5, the number of iterations (i.e., changes
c                             to the x vector) performed.
c                           see programming notes for description of
c                             the latter elements of work.
c                infer  - an integer that is set, on output, to
c                           indicate which convergence criterion was
c                           satisfied.
c                         infer = 0 indicates that convergence failed.
c                           ier gives further explanation.
c                         infer = 1 indicates that the first criterion
c                           was satisfied.
c                         infer = 2 indicates that the second criterion
c                           was satisfied.
c                         infer = 4 indicates that the third criterion
c                           was satisfied.
c                         if more than one of the convergence criteria
c                           were satisfied on the final iteration,
c                           infer contains the corresponding sum.
c                           (e.g., infer = 3 implies first and second
c                           criteria satisfied simultaneously).
c                ier    - error parameter (output)
c                         terminal error
c                           ier = 129 implies a singularity was detected
c                             in the jacobian and recovery failed.
c                           ier = 130 implies at least one of m, n, iopt,
c                             parm(1), or parm(2) was specified
c                             incorrectly.
c                           ier = 132 implies that after a successful
c                             recovery from a singular jacobian, the
c                             vector x has cycled back to the
c                             first singularity.
c                           ier = 133 implies that maxfn was exceeded.
c                         warning error
c                           ier = 38 implies that the jacobian is zero.
c                             the solution x is a stationary point.
c                           ier = 39 implies that the marquardt
c                             parameter exceeded parm(3).  this
c                             usually means that the requested
c                             accuracy was not achieved.
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - leqt1p,ludecp,luelmp,uerset,uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   copyright           - 1982 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c --- original IMSL call:
c
c     SUBROUTINE ZXSSQ  (FUNC,M,N,NSIG,EPS,DELTA,MAXFN,IOPT,PARM,
c    *                   X,SSQ,F,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
c
c --- modified version, all additional arguments are simply passed
c --- on to subroutine func. None of them is used in zxssq1 itself.
c --- the additional arguments are explained in subroutine outer_solution
c
c     specifications for arguments
c
       USE replace_imsl,  ONLY : my_uertst,my_uerset,my_leqt1p
      implicit  integer (i-n), real*8 (a-h, o-z)
      external           func
      integer            m,n,nsig,maxfn,iopt,ixjac,infer,ier
      real*8             eps,delta,parm(*),x(n),ssq,f(m),xjac(*),
     .                   xjtj(*),work(*)
c
c     xjac used internally in packed form
c     specifications for local variables
c
      integer            imjc,igrad1,igradl,igradu,idelx1,idelxl,
     .                   idelxu,iscal1,iscall,iscalu,ixnew1,ixnewl,
     .                   ixbad1,ifpl1,ifpl,ifpu,ifml1,ifml,ieval,
     .                   ibad,isw,iter,j,ijac,i,k,l,is,js,li,lj,icount,
     .                   izero,level,levold
      real*8             al,cons2,dnorm,dsq,
     .                   erl2,erl2x,f0,f0sq,f0sqs4,g,half,
     .                   hh,one,onep10,onep5,onesf0,ax,
     .                   prec,rel,rhh,sig,sqdif,ssqold,sum,ten,
     .                   tenth,xdif,xhold,up,zero,
     .                   xdabs,relcon,p01,two,huntw,delta2
c
c --- additional declarations required for modified calling arguments:
c
      integer            maxinit,nw,nh,jsymetric,map(*),iterations,
     .                   isetbdry(*),isetrect(*),iextiters
      logical            updownsym,interior
      real*8             psi1d(*),residin,omega,psi2d(nw,nh),psibdry,
     .                   rmhdgrid(*),coeff(9,*),deltar,deltarzsq,
     .                   deltarsq,residmax,zmhdgrid(*),
     .                   rplasbdry(*),zplasbdry(*),psidifbdext(*)
c
c --- end of additional declarations
c
      data               sig/14/    ! this is questionable
      data               ax/0.1/
      data               p01,tenth,half,zero,one,onep5,two,
     .                   ten,huntw,onep10/0.01,0.1,0.5,0.0,
     .                   1.0,1.5,2.0,10.0,1.2e2,1.0e10/
c
c     error checks
c
      ier   = 0
      level = 0
      call my_uerset (level, levold)
      if (m .le. 0 .or. m .gt. ixjac .or. n .le. 0 .or. iopt .lt. 0
     .                                             .or. iopt. gt. 2)
     .  go to 305
      imjc = ixjac-m
      if (iopt .ne. 2)  go to 5
      if (parm(2) .le. one .or. parm(1) .le. zero)  go to 305
c
c     machine-dependent constants
c
    5 prec   = ten**(-sig - one )
      rel    = ten**(-sig * half)
      relcon = ten**(-nsig      )
c
c     work vector is concatenation of scaled hessian,gradient,delx,scale,
c     xnew,xbad,f(x+del),f(x-del)
c
      igrad1 = ((n+1)*n)/2
      igradl = igrad1+1
      igradu = igrad1+n
      idelx1 = igradu
      idelxl = idelx1+1
      idelxu = idelx1+n
      iscal1 = idelxu
      iscall = iscal1+1
      iscalu = iscal1+n
      ixnew1 = iscalu
      ixnewl = ixnew1+1
      ixbad1 = ixnew1+n
      ifpl1  = ixbad1+n
      ifpl   = ifpl1+1
      ifpu   = ifpl1+m
      ifml1  = ifpu
      ifml   = ifml1+1
      imjc   = ixjac - m
c
c     initialize variables
c
      al = one
      cons2 = tenth
      if (iopt .eq. 0)  go to 20
      if (iopt .eq. 1)  go to 10
      al = parm(1)
      f0 = parm(2)
      up = parm(3)
      cons2 = parm(4)*half
      go to 15
   10 al = p01
      f0 = two
      up = huntw
   15 onesf0 = one/f0
      f0sq = f0*f0
      f0sqs4 = f0sq**4
   20 ieval = 0
      delta2 = delta*half
      erl2 = onep10
      ibad = -99
      isw = 1
      iter = -1
      infer = 0
      ier = 0
      do j=idelxl,idelxu
        work(j) = zero
      end do
      go to 165
c
c     main loop
c
   30 ssqold = ssq
c
c     calculate jacobian
c
      if (infer .gt. 0 .or. ijac .ge. n .or. iopt .eq. 0 .or.
     .    icount .gt. 0)  go to 55
c
c     rank one update to jacobian
c
      ijac = ijac+1
      dsq = zero
      do j=idelxl,idelxu
        dsq = dsq + work(j) * work(j)
      end do
      if (dsq .le. zero)  go to 55
      do i=1,m
         g = f(i)-work(ifml1+i)
         k = i
         do j=idelxl,idelxu
           g = g +  xjac(k)*work(j)
           k = k + ixjac
         end do
         g = g/dsq
         k = i
         do j=idelxl,idelxu
           xjac(k) = xjac(k)-g*work(j)
           k = k+ixjac
         end do
      end do
      go to 80
c
c     jacobian by incrementing x
c
   55 ijac = 0
      k = -imjc
      do 75 j=1,n
         k = k+imjc
         xdabs = ABS (x(j))
         hh    = rel * MAX (xdabs, ax)
         xhold = x(j)
         x(j) = x(j)+hh
****     CALL FUNC (X,M,N,WORK(IFPL))
         call func (x,m,n,work(ifpl),
     .              zmhdgrid,psidifbdext,iextiters,
     .              residin,maxinit,iterations,updownsym,
     .              jsymetric,nw,nh,map,psi1d,coeff,deltar,
     .              deltarzsq,deltarsq,interior,omega,psi2d,
     .              rmhdgrid,psibdry,residmax,
     .              isetbdry,isetrect,rplasbdry,zplasbdry)
         ieval = ieval+1
         x(j)  = xhold
         if (isw .eq. 1)  go to 65
c
c        central differences
c
         x(j) = xhold-hh
****     CALL FUNC (X,M,N,WORK(IFML))
         call func (x,m,n,work(ifml),
     .                   zmhdgrid,psidifbdext,iextiters,
     .                   residin,maxinit,iterations,updownsym,
     .                   jsymetric,nw,nh,map,psi1d,coeff,deltar,
     .                   deltarzsq,deltarsq,interior,omega,psi2d,
     .                   rmhdgrid,psibdry,residmax,
     .                   isetbdry,isetrect,rplasbdry,zplasbdry)
         ieval = ieval+1
         x(j) = xhold
         rhh = half/hh
         do i=ifpl,ifpu
           k = k+1
           xjac(k) = (work(i)-work(i+m))*rhh
         end do
         go to 75
c
c        forward differences
c
   65    rhh = one/hh
         do i=1,m
           k = k+1
           xjac(k) = (work(ifpl1+i)-f(i))*rhh
         end do
   75 continue
c
c     calculate gradient
c
   80 erl2x     = erl2
      erl2      = zero
      k         = -imjc
      do j=igradl,igradu
        k       = k + imjc
        sum     = zero
        do i=1,m
          k     = k + 1
          sum   = sum+xjac(k)*f(i)
        end do
        work(j) = sum
        erl2    = erl2+sum*sum
      end do
      erl2      = SQRT (erl2)
c
c     convergence test for norm of gradient
c
      if (ijac .gt. 0     )  go to 95
      if (erl2 .le. delta2)  infer = infer + 4
      if (erl2 .le. cons2 )  isw   = 2
c
c     calculate the lower super triangle of jacobian (transposed) * jacobian
c
   95 l  = 0
      is = -ixjac
      do i=1,n
         is = is+ixjac
         js = -ixjac
         do j=1,i
           js  = js+ixjac
           l   = l+1
           sum = zero
           do k=1,m
             li  = is+k
             lj  = js+k
             sum = sum+xjac(li)*xjac(lj)
           end do
           xjtj(l) = sum
        end do
      end do
c
c     convergence checks
c
      if (infer .gt. 0    )  go to 315
      if (ieval .ge. maxfn)  go to 290
c
c     compute scaling vector
c
      if (iopt .eq. 0)  go to 120
      k = 0
      do j=1,n
        k              = k+j
        work(iscal1+j) = xjtj(k)
      end do
      go to 135
c
c     compute scaling vector and norm
c
  120 dnorm = zero
      k     = 0
      do j=1,n
        k              = k + j
        work(iscal1+j) = SQRT (xjtj(k))
        dnorm          = dnorm+xjtj(k)*xjtj(k)
      end do
      dnorm = one / SQRT (dnorm)
c
c     normalize scaling vector
c
      do j=iscall,iscalu
        work(j) = work(j)*dnorm*erl2
      end do
c
c     add l-m factor to diagonal
c
  135 icount = 0
  140 k = 0
      do i=1,n
        do j=1,i
          k       = k+1
          work(k) = xjtj(k)
        end do
        work(k)        = work(k)+work(iscal1+i)*al
        work(idelx1+i) = work(igrad1+i)
      end do
c
c     cholesky decomposition
c
  155 call my_leqt1p (work,1,n,work(idelxl),n,0,g,xhold,ier)
      if (ier .eq. 0)  go to 160
      ier = 0
      if (ijac .gt.  0 )  go to 55
      if (ibad .le.  0 )  go to 240
      if (ibad .ge.  2 )  go to 310
      go to 190
  160 if (ibad .ne. -99)  ibad = 0
c
c     calculate sum of squares
c
  165 do j=1,n
        work(ixnew1+j) = x(j)-work(idelx1+j)
      end do
c
****  CALL FUNC (WORK(IXNEWL),M,N,WORK(IFPL))
      call func (work(ixnewl),m,n,work(ifpl),
     .           zmhdgrid,psidifbdext,iextiters,
     .           residin,maxinit,iterations,updownsym,
     .           jsymetric,nw,nh,map,psi1d,coeff,deltar,
     .           deltarzsq,deltarsq,interior,omega,psi2d,
     .           rmhdgrid,psibdry,residmax,
     .           isetbdry,isetrect,rplasbdry,zplasbdry)
      ieval = ieval+1
      ssq   = zero
      do i=ifpl,ifpu
        ssq = ssq+work(i)*work(i)
      end do
      if (iter .ge. 0)  go to 185
c
c     ssq for initial estimates of x
c
      iter = 0
      ssqold = ssq
      do i=1,m
        f(i) = work(ifpl1+i)
      end do
      go to 55
  185 if (iopt .eq. 0)  go to 215
c
c     check descent property
c
      if (ssq .le. ssqold)  go to 205
c
c     increase parameter and try again
c
  190 icount = icount+1
      al     = al*f0sq
      if (ijac .eq. 0)  go to 195
      if (icount .ge. 4 .or. al .gt. up)  go to 200
  195 if (  al .le. up)  go to 140
      if (ibad .eq. 1 )  go to 310
      ier = 39
      go to 315
  200 al = al / f0sqs4
      go to 55
c
c     adjust marquardt parameter
c
  205 if (icount .eq. 0) al = al/f0
      if (erl2x .le. zero)  go to 210
      g = erl2/erl2x
      if (erl2 .lt. erl2x)  al = al * MAX (onesf0, g)
      if (erl2 .gt. erl2x)  al = al * MIN (f0    , g)
  210 al = MAX (al, prec)
c
c     one iteration cycle completed
c
  215 iter = iter+1
      do j=1,n
        x(j) = work(ixnew1+j)
      end do
      do i=1,m
        work(ifml1+i) = f(i)
        f(i)          = work(ifpl1+i)
      end do
c
c     relative convergence test for x
c
      if (al .gt. 5.0)  go to 30
      do j=1,n
        xdif = ABS (work(idelx1+j)) / MAX (ABS (x(j)), ax)
        if (xdif .gt. relcon)  go to 235
      end do
      infer = 1
c
c     relative convergence test for ssq
c
  235 sqdif = ABS (ssq-ssqold) / MAX (ssqold, ax)
      if (sqdif .le. eps)  infer = infer+2
      go to 30
c
c     singular decomposition
c
  240 if (ibad)  255, 245, 265
c
c     check to see if current iterate has cycled back to the last singular point
c
  245 do j=1,n
        xhold = work(ixbad1+j)
        if (ABS (x(j)-xhold) .gt. relcon * MAX (ax, ABS (xhold)))
     .  go to 255
      end do
      go to 295
c
c     update the bad x values
c
  255 do j=1,n
        work(ixbad1+j) = x(j)
      end do
      ibad = 1
c
c     increase diagonal of hessian
c
  265 if (iopt .ne. 0)  go to 280
      k = 0
      do i=1,n
        do j=1,i
          k       = k + 1
          work(k) = xjtj(k)
        end do
        work(k) = onep5*(xjtj(k)+al*erl2*work(iscal1+i))+rel
      end do
      ibad = 2
      go to 155
c
c     replace zeroes on hessian diagonal
c
  280 izero = 0
      do j=iscall,iscalu
        if (work(j) .le. zero) then
          izero   = izero+1
          work(j) = one
        end if
      end do
      if (izero .lt. n)  go to 140
      ier = 38
      go to 315
c
c     terminal error
c
  290 ier = ier + 1
  295 ier = ier + 1
      ier = ier + 1
  305 ier = ier + 1
  310 ier = ier + 129
      if (ier .eq. 130)  go to 335
c
c     output erl2, ieval, nsig, al, and iter
c
  315 g = sig
      do j=1,n
        xhold = ABS (work(idelx1+j))
        if (xhold .gt. zero) then
          g   = MIN (g, -LOG10 (xhold) + LOG10 (MAX (ax, ABS (x(j)))))
        end if
      end do
      if (n .gt. 2)  go to 330
      do j=1,n
        work(j+5) = work(j+igrad1)
      end do
  330 work(1) = erl2+erl2
      work(2) = ieval
      ssq     = ssqold
      work(3) = g
      work(4) = al
      work(5) = iter
  335 call my_uerset (levold, levold)
      if (ier .eq. 0)  go to 9005
 9000 continue
      call my_uertst (ier, 'zxssq1')
 9005 return
c
      end
