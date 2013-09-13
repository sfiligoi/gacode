  MODULE contour 
! ---------------------------------------------------------------------------
! -- supplies  plasma contour related routines
! -----------------------------------------------------------------HSJ-------

      USE nrtype,                                  ONLY : DP,I4b
      USE error_handler

      USE io_gcnmp,                                ONLY : nlog,ncrt

      USE MPI_data,                                ONLY : myid,master

      USE common_constants,                        ONLY : izero,zeroc,twopi

      USE replace_imsl,                            ONLY : my_dbcevl1

      LOGICAL  use_cnt1, use_cnt2
      DATA use_cnt1,use_cnt2      /.TRUE.,.FALSE. /

      CONTAINS 

      SUBROUTINE cntour (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin, &
                        yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,   &
                        xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,  &
                        nh,cspln,n2cspln,nh2,itty,iptsm,ierr,           &
                        bpmag,iconvg,delta_psi )
! ------------------------------------------------------------------------
! -- INTERFACE TO CONTOUR ROUTINES 
! ---------------------------------------------------------HSJ------------

      IMPLICIT  NONE

      LOGICAL        use_cnt1, use_cnt2

      REAL(DP)   xc(*), yc(*), x(nw), y(nh),cspln(n2cspln,nw,nh2),     &
                 bpmag(*),xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,     &
                 yxmin,yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,   &
                 xmax,ymin,ymax,delta_psi

      INTEGER(I4B) nw,nh,nh2,n2cspln,ipts,itty,iptsm,ierr,iconvg,      &
                  iauto,iauto_save,iautoc

      iauto_save = iauto
      use_cnt1     = .TRUE.
      use_cnt2     = .FALSE.

   10 ierr       = izero
      iauto      = iauto_save
      iautoc     = izero

      IF (use_cnt1) THEN
        CALL cntour1 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin,     &
                     yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,       &
                     xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,      &
                     nh,cspln,n2cspln,nh2,itty,iptsm,ierr,bpmag,iconvg,  &
                     delta_psi)

      ELSE IF (use_cnt2) THEN
        CALL cntour2 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin,     &
                     yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,       &
                     xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,      &
                     nh,cspln,n2cspln,nh2,itty,iptsm,ierr,bpmag,iconvg,  &
                     delta_psi)
      ELSE
        IF(myid == master)     &
           WRITE (ncrt, 4)
    4      FORMAT (/ ' subroutine CNTOUR detects an error:'   / &
                  ' must select cntour1 or cntour2 with use_cnt'        )
        lerrno = iomaxerr + 210_I4B
        CALL terminate(lerrno,nlog)
      END IF

      IF (ierr .NE. izero) THEN
!       CALL interface_dump_psi (psivl)     ! dump the psi values
        IF(myid == master) THEN
           IF (use_cnt1) &
           WRITE (ncrt, 6)
     6     FORMAT(/, 'Sub cntour1  has returned ierr ne 0')
           IF (use_cnt2) &
           WRITE (ncrt, 7)
     7     FORMAT(/, 'Sub cntour2  has returned ierr ne 0')
           lerrno = iomaxerr + 211_I4B
           CALL terminate(lerrno,nlog)
         ENDIF
      ENDIF

      RETURN

      END       SUBROUTINE cntour


      SUBROUTINE cntour1 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin, &
                          yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin, &
                          xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y, &
                          nh,cspln,n2cspln,nh2,itty,iptsm,ierr, &
                          bpmag,iconvg,delta_psi)
!
! ----------------------------------------------------------------------
!
!                     version 5/17/94
!                     version 7/07/94
!             (same as 5/17/94 but uses IMSL spline routines)
!
!          CNTOUR1 AND ITS OLDER VERSION, CNTOUR2, ARE NOT INTENDED TO
!          TO BE GENERAL-PURPOSE CONTOURING ROUTINES. THESE ROUTINES
!          ARE DESIGNED SPECIFICALLY FOR FINDING PLASMA FLUX
!          SURFACES AND ASSUME THAT:
!
!          1) THE PLASMA SURFACE IS CONVEX (BEAN SHAPES MAY BE
!             PROBLEMATIC)
!          2) THE FUNCTION TO BE CONTOURED IS MONOTONICALLY
!             DECREASING IN THE REGION OF INTEREST
!
!     this routine uses the bicubic spline representation of psi to get
!     psi values at an arbitrary point on the MHD grid. the spline
!     coefficients must exist (in cspln) before this subroutine is called.
!     given (xaxd,yaxd) and a psi value, psivl, generate a contour of
!     ordered points,(xc(i),yc(i)), i = 1,ipts, which has (xaxd,yaxd) as
!     as an interior point. the contour must fully encircle (xaxd,yaxd).
!     the search is limited to a rectangle specified by (xmin,xmax,ymin
!     ymax).  dx and dy determine the basic increment in x and y for
!     the coarse grid search.  they should be picked so that over the
!     distance dx or dy, psi is single valued.  a sliding interval
!     search rather than a binary search is used for the coarse grid so
!     that non-monotonic psi can be handled.
!     arcl is the arclength in meters taken from the current point to get
!     to the next point.
!     bperr is relative change in
!     poloidal b field between (xc(i),yc(i)) and (xc(i+1),yc(i+1)).
!     if the number of elements (xc,yc) exceeds the limit,iptsm,
!     bperr is relaxed (up to twice) by a constant increment
!     dbperr before an error exit is taken.
!     bperr = 0.03 is suggested.
!     (xemin,xemax,yemin,
!     yemax) are min and max values of contour.  (yxmin,yxmax), are y
!     values at x = xemin and x = xemax respectively.  (xymin,xymax) are x
!     values at y = yemin and y = yemax respectively.  iauto is used as an
!     "error recovery switch".  that is, if a flux surface is found to
!     pass outside the search box an error exit is taken by way of label
!     1000, if iauto = 0.  if iauto = 1, the input value of psi, psivl, is
!     increased toward the value of psi at (xaxd,yaxd) until a closed
!     surface inside the search box is found.  obviously this option
!     makes sense only for the limiter flux surface.  it could be used
!     to find the plasma boundary if an appopriate search box were
!     defined.  however its primary use here is to slightly modify
!     psilim so that a closed limiter flux surface is found.  the option
!     is necessary due to the fact some codes generate eqdsks
!     using linear interpolation, which differs from the bicubic inter-
!     polation used here.  the user is informed that psivl was changed
!     by returning iautoc = 1.  if no change occurred iautoc = 0.
!     if no error ocurred ierr = 0 is returned. if an error ocurred then
!     ierr = 1 is returned and the results of this calculation are not valid.
!
! --- input (see description above)
!
!  xaxd
!  yaxd
!  psivl   psi value for which contour is desired
!  arcl    arc length of step (in same units of length as xaxd for example)
!  bperr   arcl,bperr control point density on contour returned
!  dx
!  dy      defines step length taken
!          input dx and dy as 0.0 to set step length in this subroutine
!  xmin
!  xmax
!  ymin
!  ymax    define the search box
!             note: the search box could be defined as the entire
!                   MHD grid (i.e., xmin = x(1),xmax = x(nw),ymin = y(1),
!                   ymax = y(nh)) however this will normally not work
!                   satisfactorily because there are local minima and
!                   maxima around the f coils which would confuse the
!                   contour tracing and require many more calculations.
!                   (these additional calculations are not done here!)
!                   the best way to define the search box is to set it
!                   equal to a rectangle,within which the plasma must
!                   reside under all circumstances. the only such
!                   rectangle to suggest itself is the one defined by
!                   the extremes of the limiter. hence it is suggested
!                   that xmin be set to the min x value of the limiter
!                   points vector,xmax the max x value of limiter points
!                   vector and similarly for ymin,ymax.
!  iaouto   = 1      automatic error recovery, = 0 no recovery
!                  (iauto may be set to 0 in this subroutine
!                  if a failure occurs)
!                   There are three ways to run
!                     a)iauto=0    should be ok for all interior plasma
!                                  surfaces. The value of psi to be traced,
!                                  psivl,is not altered. Either the routine
!                                  returns with the desired contour or
!                                  a failure is indicated with ierr
!                                  (iconvg is not used in this mode)
!
!                        the other two modes are intended for use on the
!                        plasma surface:
!                     b)iauto=1,iconvg=0 . May change the value of psivl!!
!                       (hence use call by reference). If the contour
!                       with the value of psivl could not be traced or
!                       passes outside the "box" (see below) the value
!                       of psivl will be REPEATEDLY  brought closer
!                       to the value
!                       on the magnetic axis. Either by some fraction
!                       of delta_psi if it is set or some fraction of
!                       psiaxis-psilim if it is not. The first contour
!                       succesfully traced will be returned.
!                     c)iauto=1,iconvg=1 as b except that a binary search
!                       is done to find  the plasma surface. FOR THIS
!                       PURPOSE THE "PLASMA SURFACE" IS DEFINED AS THE
!                       LARGEST CLOSED FLUX SURFACE THAT CAN BE FOUND
!                       WHICH DOES NOT PASS OUTSIDE THE SEARCH BOX.
!                       The actual search box used is not critical for
!                       diverted plasmas but is very important for
!                       limited plasmas. It is recommended that the
!                       boundary finder "bound" be used in those cases
!                       where an appropriate psivl for the plasma
!                       boundary is not known a priori.
!
!  x(1..nw)
!  y(1..nh)        the MHD grid vectors
!  cspln(nw,nh,4)  bicubic spline coefficient array
!  cspln(n2cspln,nw,nh2)  for IMSL routines
!  itty            iounit for diagnostic messages. set to 0 to suppress.
!  iptsm           max allowed number of points to be returned in
!                  vectors xc,yc
!  iconvg          used only if iauto = 1 is set . if iauto = 1 and iconvg = 1
!                  the contour is converged to the outermost cntour
!                  (meaning the contour with the smallest psi) that never
!                  leaves the search box given by xmin,xmax,ymin,ymax.
! delta_psi        this is the psi grid spacing. Define delta_psi as
!                  psi(penultimate value)-psi(last(i.e., edge) value)
!                  If the value of psi to be traced leads to an error
!                  and iauto=1 then delta_psi is used to control the
!                  amount by which the value of psi is adjusted. Otherwise
!                  delta_psi is not used,so you set its value accordingly.
!
! --- output:
!
!  xemin
!  xemax
!  yemin
!  yemax
!  yxmin
!  yxmax
!  xymax
!  xymin             limits of contour (see above)
!  xc(1,..ipts)
!  yc(1,..ipts)      the ipts contour points found
!  bpmag(1,...,ipts) magnitude of poloidal bfield at xc, yc
!  ipts
!  ierr     error flag, = 0 no error, = 1 error, output not usable.
!           if ierr > 1 then ierr equals the number of points found
!           on the contour up to this point. since this number will
!           be exceeded if the code continues,the subroutine must
!           stop. the error condition on bperr has allready been relaxed
!           before this condition causes the subroutine to terminate.
!           hence in order to continue this subroutine could be called
!           again with a larger value of arcl.
!           (note: in practice we have not found this to be a problem,
!           provided that iptsm is large enough to begin with to give
!           a reasonable representation of the contour).
!
! ------------------------------------------------------------------ HSJ

!

      IMPLICIT  NONE
!
      REAL(DP)   pds(6), xc(*), yc(*), x(nw), y(nh),             &
                 cspln(n2cspln,nw,nh2), bpmag(*)
      REAL(DP)   fpiov4,spiov4, twopi

      REAL(DP)   xemin ,xemax,yemin ,yemax ,xymin,xymax,         &
                 yxmin , yxmax ,psiaxd,xmin,xmax,ymin,ymax,x2,   &
                 y2,ymult,sint,cost,xns,yns,thetnew,psival,dang, &
                 bperr,dx,dy,psivl,delta_psi,bpmintol,psiin,     &
                 bperrsav,dbperr,dthetmin,dthet,dxx,dyy,dxmin,   &
                 dymin,serrt,derrt,x1,y1,psi1,psi2,cmult,        &
                 dpsi,ysave,arclmin,bp1,bp2,sintsq,costsq,       &
                 dtarcl,xsave,psisave,xmult,dpsids,psivl0,dapsi, &
                 serr,dmult,tstart,psivln,psiout,delx,dely,step, &
                 estmrad


      REAL(DP)   xaxd, yaxd, thet, a, bincp,tend,xn,yn,arcl

      INTEGER(I4B) nitr8,nstepc,itry,irestart,ierr,iautoc,ilsave,&
                   jlsave,ier,nothing,nw,nh,n2cspln,nh2, iflg,   &
                   isgn,iauto,ipts,iptsm,iconvg,nsteps,j,        &
                   maxtry,newti,iqpt,isgnsave,itty


      DATA       fpiov4     , spiov4     , twopi      , nitr8 &
               / 3.926990818, 5.497787145, 6.283185307, 10 /

      LOGICAL    bptol, reverse
!
      nstepc = 2  ! # steps along ray before tangent line method is used
      bpmintol = 0.05 ! if min bp drops below bpmintol and iauto = 1..
!                     ..the contour is searched for an x point
      itry     = 0
      irestart = 0
      psiin    = 0.0
      ierr     = 0
      iautoc   = 0
      bperrsav = bperr
      dbperr   = bperr/2.0
      arclmin  = 0.10 * arcl ! used to set minimum angle increment below
      ilsave   = 1
      jlsave   = 1
!
! --- get psi at (xaxd,yaxd)
!
!***  call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,xaxd,yaxd,pds,ier,
!***.                           1,ilsave,jlsave)
!
!      CALL dbcevl1 (x, nw, y, nh, cspln, nw, xaxd, yaxd, pds, ier, 1)
      CALL my_dbcevl1 (x, nw, y, nh, cspln, nw, xaxd, yaxd, pds, ier, 1)
      IF (ier .NE. 0) THEN
        IF(myid == master)     &
           WRITE (itty, 4) ier
    4      FORMAT (/ ' subroutine CNTOUR1,contour_module,f90, detects an error:'   / &
                  ' EVAL_BICUBIC_SPLINE returned ier =', i5 / &
                  ' while getting psi on axis'               )
        lerrno = iomaxerr + 206_I4B
        CALL terminate(lerrno,nlog)
      END IF
      psiaxd = pds(1)
!
      IF (psiaxd .LT. psivl) THEN
        IF(myid == master)                                        &
        WRITE (itty, 3)  psiaxd, psivl
    3   FORMAT (/ ' subroutine CNTOUR1 detects an ERROR:'       / &
                  ' psi on magnetic axis    ', 1pe12.6          / &
                  ' psi contour to be found ', 1pe12.6          / &
                  ' must have max psi on axis, non-recoverable' / &
                  ' execution must stop'                         )
        lerrno = iomaxerr + 207_I4B
        CALL terminate(lerrno,nlog)
      END IF
!
    5 xemin    =  1.0e10     ! entry point for generating contour..
      xemax    = -1.0e10     ! ..come back to here if a correctable..
      yemin    =  1.0e10     ! ..failure occured.
      yemax    = -1.0e10
      nothing  = -1
      xymin    =  0.0
      xymax    =  0.0
      yxmin    =  0.0
      yxmax    =  0.0
      ipts     =  0
      irestart = irestart + 1
      tstart   = fpiov4
      tend     = tstart + twopi
      thet     = tstart
      dthetmin = 1.0e-5   ! min. increment in theta from one ray to next
!                           note that this is increased below,
!                           after the first point is found
      dthet = 0.0
      dxx   = dx
      dyy   = dy
      IF (dxx .EQ. 0.0)  dxx = x(3)-x(1)  ! 2.0 times grid spacing is OK
      IF (dyy .EQ. 0.0)  dyy = y(3)-y(1)
      dxmin = 0.05*(x(2)-x(1)) ! subroutine stops if dxx becomes < dxmin
      dymin = 0.05*(y(2)-y(1)) ! subroutine stops if dyy becomes < dymin
      dx    = dxx
      dy    = dyy
      IF (dxmin .GE. dxx)  dxmin = 0.5 * dxx
      serrt = 3.5e-06
      derrt = 0.5e-07
!
! --- serrt is absolute error convergence criteria for Newton's method below
!
! --- loop over theta from tstart to tend ( = tstart+twopi)
!
   10 thet   = thet + dthet
      IF (thet .GT. twopi .AND. tstart .NE. 0.0) THEN
        thet = thet - twopi
        tend = tend - twopi
      END IF
      thet   = MIN (thet, tend)
      IF (thet .EQ. tend)  go to 200   ! normal exit from routine
      dxx    = dx
      dyy    = dy
      iqpt   = 0
!
! --- get equation of ray emanating from (xaxd,yaxd)
!
      CALL findeqlin (xaxd, yaxd, thet, a, bincp, iflg, isgn)
!
! --- now have y = a*x+bincp    (iflg = 0)   or
! ---          x = a*y+bincp    (iflg = 1)
!
!     start search from axis point
!
      x1   = xaxd
      y1   = yaxd
      bp1  = 0.0D0
      psi1 = psiaxd
      cost = COS (thet)
      sint = SIN (thet)
!
! --- sliding interval search. max width of interval ~1.41*(dx or dy)
!
      nsteps = 0
      cmult  = 1.0
      xsave  = -1.0e20
   40 IF (iflg .EQ. 1)  go to 50
!
! --- search in x
!
      xmult = 1.0
   41 x2    = x1+isgn*dxx*xmult*cmult    ! increment x
      IF (x2 .EQ. x1)  go to 1000        ! psivl is not within reach
      IF (ABS(x2- x1)  .LT. 1.e-12 ) go to 1000
      x2    = MAX (x2,xmin)
      x2    = MIN (x2,xmax)              ! limit the search
      IF (x2 .EQ. x1)  go to 1000        ! psivl is not within reach
      y2    = a*x2+bincp                 ! get corresponding value of y
      IF (y2 .LT. ymin .OR. y2 .GT. ymax) THEN
        xmult = 0.5 * xmult
        go to 41
      END IF
      go to 60
!
! --- search in y
!
   50 ymult = 1.0_DP
   51 y2    = y1+isgn*dyy*ymult*cmult    ! increment
      IF (y2 .EQ. y1)  go to 1000        ! psivl is not within reach
      IF (ABS(y2- y1)  .LT. 1.e-12)  go to 1000  ! psivl is not within r
      y2    = MAX (y2,ymin)
      y2    = MIN (y2,ymax)              ! limit the search
      IF (y2 .EQ. y1)  go to 1000        ! psivl is not within reach
      x2    = a*y2+bincp                 ! get corresponding value of x
      IF (x2 .LT. xmin .OR. x2 .GT. xmax) THEN
          ymult = ymult * 0.5
          go to 51
      END IF
!
!* 60 call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,x2,y2,pds,ier,1,
!*** .                          ilsave,jlsave)
!   60 CALL dbcevl1 (x,nw,y,nh,cspln,nw,x2,y2,pds,ier,1)
   60 CALL my_dbcevl1 (x,nw,y,nh,cspln,nw,x2,y2,pds,ier,1)
      IF (ier .NE. 0) THEN
        IF(myid == master)THEN
           WRITE (itty, *)  'CNTOUR1 reports, call to EVAL_BI.., ier =', ier
           WRITE (itty, *)  'arguments    x2, y2 =', x2  , y2
           WRITE (itty, *)  'xgrid(1), xgrid(nw) =', x(1), x(nw)
           WRITE (itty, *)  'ygrid(1), ygrid(nh) =', y(1), y(nh)
        ENDIF
!        CALL STOP ('subroutine CNTOUR1: problem #3', 159)
        lerrno = 216 + iomaxerr
        call terminate(lerrno,nlog)
      END IF
      psi2 =  pds(1)
      dpsi = (psivl-psi1)*(psivl-psi2)
      IF (dpsi .LE. 0.0)  go to 70
!
! ----------------------------------------------------------------------
!
      IF (cmult .LT. 0.99) THEN
        cmult = cmult * 0.5
        go to 40
      END IF
!
! ----------------------------------------------------------------------
! --- if psi starts to increase we are in the vicinity of the saddle
! --- (i.e., x point). in this case we want to be sure that we stay
! --- on the contour which envelopes the plasma. the only robust way
! --- to do this is to first search for the minimum in psi along
! --- the ray. onece this minimum is found we guarantee that the
! --- point we find is on the plasma contour by not allowing the ray
! --- to extend past this minimum.
! --- to find the minimum we use Newton's method to search for the
! --- point along the ray where the directional derivative of pis
! --- along the ray is zero.
! ----------------------------------------------------------------------
!
      reverse = .FALSE.
      IF (psi2-psi1 .GT. 0.0) THEN    ! psi is increasing
         reverse  = .TRUE.
         xsave    = x1   ! save current guess for possible later restore
         ysave    = y1
         psisave  = psi1
         isgnsave = isgn
         xn       = x1
         yn       = y1
         costsq   = cost * cost
         sintsq   = sint * sint
         newti    = 0
         step     = serrt + 2._DP
         DO WHILE (ABS (step) .GT. serrt .AND. newti .LE. nitr8)
!***         call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,xn,yn,pds,
!*** .                                 ier,itty,ilsave,jlsave)
!             CALL dbcevl1 (x,nw,y,nh,cspln,nw,xn,yn,pds,ier,6)
            CALL my_dbcevl1 (x,nw,y,nh,cspln,nw,xn,yn,pds,ier,6)
             step  = -(pds(2)*cost+pds(3)*sint)/(costsq*pds(5) &
                        + 2.0*cost*sint*pds(4) + sintsq*pds(6))
             xn    = xn+step*cost
             yn    = yn+step*sint
             newti = newti+1
         END DO
         IF (newti .LE. nitr8) THEN          ! Newton's method converged
             isgn = -isgn
             IF (iflg .EQ. 1)  dyy = 0.5 * dyy
             IF (iflg .EQ. 0)  dxx = 0.5 * dxx
             x2   = xn
             y2   = yn
             psi2 = pds(1)
!
!            we have found the minimum in psi along the ray
!            this minimum must be less than or equal to the value for
!            which we are searching,psivl. If this is not the case then
!            we can't find psivl on the ray so abandon the search:
!
             IF (psi2 .GT. psivl)  go to 1000
         ELSE
           go to 1500            ! Newton's method failed, try to refine
         END IF
       END IF
       x1     = x2
       y1     = y2
       psi1   = psi2
       nsteps = nsteps+1
       go to 40
!
! ----------------------------------------------------------------------
!
! --- now have psivl between psi1 and psi2, converge using Newton-Raphson
!
   70 newti = 0
      IF (iflg .EQ. 1)  go to 75
      xn = x1+isgn*dxx * 0.5
      yn = a*xn+bincp
      go to 80
   75 yn = y1+isgn*dyy * 0.5
      xn = a*yn+bincp
!* 80 call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,xn,yn,pds,ier,3,
!*** .                          ilsave,jlsave)
!  80 CALL dbcevl1(x,nw,y,nh,cspln,nw,xn,yn,pds,ier,3)
  80 CALL my_dbcevl1(x,nw,y,nh,cspln,nw,xn,yn,pds,ier,3)
      IF (ier .NE. 0) THEN
        newti = nitr8
        go to 76
      END IF
      dpsids =  pds(2)*cost+pds(3)*sint
      dpsi   =  pds(1)-psivl
      serr   = -dpsi/dpsids
      IF (ABS (serr) .LT. serrt)  go to 90
      IF (psivl .EQ. zeroc.AND. ABS (dpsi) .LT. derrt)  go to 90
      IF (psivl .NE. 0.0) THEN
        IF (ABS (dpsi/psivl) .LT. derrt)  go to 90
      END IF
      delx  = serr*cost
      dely  = serr*sint
      xn    = xn+delx
      yn    = yn+dely
      newti = newti+1
   76 IF (newti .GE. nitr8) THEN
!***       if ( dxx .le. dxmin)  stop 'newti'
!***       if ( dyy .le. dymin)  stop 'newti'
           IF ( dxx .LE. dxmin)  go to 90
           IF ( dyy .LE. dymin)  go to 90
           IF (iflg .EQ. 0    )  dxx = 0.5 * dxx
           IF (iflg .EQ. 1    )  dyy = 0.5 * dyy
           IF (iqpt .EQ. 0    )  go to 40           ! 6/11/93
           thet = thet - dthet
           go to 10
      END IF
      go to 80
!
! ----------------------------------------------------------------------
! --- end of Newton iteration
! --- check for sufficient accuracy in point spacing as determined by thet
! --- accuracy test is based on a relative error in poloidal b field of bperr
! ----------------------------------------------------------------------
!
   90 bp2 = SQRT (pds(2)**2 + pds(3)**2) / xn
      IF (thet .EQ. tstart)  go to 100          ! first point of contour
!
! --- if relative change in bp is less than bperr store the point just found:
!
      IF (ABS (bp2-bp1) / MAX (bp2,bp1) .LT. bperr)  go to 100
!
! --- if dthet is currently at its smallest allowed value of dthetmin
! --- and relative change in bp is still greater than bperr continue anyway:
!
      IF (dthet .EQ. dthetmin)  go to 100
!
! --- spacing too large for grad psi. decrease theta and try again
!
!***  thet = thet-dthet             ! reset to previous value
!
      dthet = dthet * 0.5           ! take half the angle increment
      dthet = MAX (dthet,dthetmin)  ! dthet lower limit must be observed
!
! --- we found an acceptable point. collect it and set up for next point
!
  100 ipts = ipts+1
      IF (ipts .LE. iptsm-1)  go to 150
!
! --- we ran out of storage for the contour points. assuming iptsm is
! --- set to a reasonable value this should not happen unless bperr is
! --- set unrealistically small. we increase bperr by dbperr and try
! --- to generate the contour again with the relaxed error requirement.
!
      IF (bperr-bperrsav .GE. 2.0*dbperr)  go to 110
      bperr = bperr+dbperr
      go to 5


  110     CONTINUE

          IF(myid == master)                                        &
              WRITE  (itty, 1040)  iptsm,dtarcl,bperr,psivl
 1040     FORMAT &
                 ('  more points on plasma '                                     / &
                  '  boundary were found than dimension of xcontr,ycontr allows' / &
                  '  allowed number = ',i5                                       / &
                  '  adjustment in point spacing will be made', &
                  '  dtarcl,bperr,psivl =',5(2x,1pe12.4))
          lerrno = iomaxerr + 207_I4B
          CALL terminate(lerrno,nlog)

  150 xc(ipts)    = xn  ! found  new point; save and set up for next ray
      yc(ipts)    = yn
      bpmag(ipts) = bp2
      xemin       = MIN (xemin,xn)
      xemax       = MAX (xemax,xn)
      IF (xemax .EQ. xn)  yxmax = yn
      IF (xemin .EQ. xn)  yxmin = yn
      yemin       = MIN (yemin, yn)
      yemax       = MAX (yemax, yn)
      IF (yemax .EQ. yn)  xymax = xn
      IF (yemin .EQ. yn)  xymin = xn
!
      IF ( ipts .EQ. 1 ) THEN  ! get a radius estimate so we can set..
!                              ..the theta increment based on arc length
        estmrad  = SQRT ((xn-xaxd)**2+(yn-yaxd)**2)
        dthetmin = arclmin/estmrad
        dtarcl   = arcl/estmrad
      END IF
!
      IF (ABS (bp2-bp1) / MAX (bp2,bp1) .LE. 0.25 * bperr) &
      dthet = 2.0 * dthet
      bp1   = bp2                    ! needed for error test on next ray
      dthet = MIN (dthet, dtarcl  )  ! don't let dthet get too large
      dthet = MAX (dthet, dthetmin)  ! don't let dthet get too small
!
! ----------------------------------------------------------------------
! --- for contours sufficiently far removed from xaxd,yaxd,that the
! --- search along the ray is expensive we use the following alternative.
! --- an approximation to the new point,(xns,yns)is found by moving along
! --- the tangent line, for a distance arcl.  see subroutine NEW_POINT.
! --- if we are near an x point, as signaled by bpmag < bpmintol,
! --- then we switch back to the ray method, with small
! --- dx,dy increments to avoid crossing the x point.
! ----------------------------------------------------------------------
!
      bptol = bpmag(ipts) .GT. bpmintol
      IF (nsteps .GT. nstepc .AND. bptol) THEN
!
!          use tangent line method:
!
           CALL new_point (pds,arcl,xn,yn,thet,twopi,spiov4,xmin,xmax, &
                           ymin,ymax,yaxd,xaxd,ier,sint,cost, &
                           xns,yns,thetnew)
           IF (ier .GT. 0)  go to 10    ! failed, go back use ray method
           dthet = thetnew -  thet
           thet  = thet    + dthet      ! do not use thetnew here
           IF (thet .GT. twopi .AND. tstart .NE. 0.0) THEN
             thet = thet - twopi
             tend = tend - twopi
           END IF
           thet   = MIN (thet, tend)
           IF (thet .EQ. tend) THEN
             nothing = 1
             go to 200
           END IF
           sint  = SIN (thet)
           cost  = COS (thet)
           xn    = xns
           yn    = yns
           newti = 0
           iqpt  = 1
           go to 80                   ! skip directly to Newton's method
      ELSE                            ! use ray method
           go to 10
      END IF
!
! --- done with this contour. close contour
!
  200 ipts        = ipts+1
      xc(ipts)    = xc(1)
      yc(ipts)    = yc(1)
      bpmag(ipts) = bpmag(1)
!
      bperr = bperrsav
      IF (itry .NE. 0 .AND. iconvg .EQ. 1) THEN
        psiin  = psivl
        psivln = 0.5 * (psiout+psiin)
        IF (ABS (psivln-psivl) .LT. 1.0e-6) THEN
          IF (itty .NE. 0) THEN
            IF(myid == master) WRITE  (itty, 1045) psivln
 1045       FORMAT (' final converged value of psi = ', 1pe16.8)
          END IF
          RETURN
        ELSE
          psivl = psivln
          go to 5
        END IF
      END IF
      RETURN
!
! --- errors

! 1000 CALL dbcevl1 (x, nw, y, nh, cspln, nw, x2, y2, pds, ier, 1)
 1000 CALL my_dbcevl1 (x, nw, y, nh, cspln, nw, x2, y2, pds, ier, 1)
      IF (itty .NE. 0) THEN
        IF(myid == master)WRITE  (itty, 1010) xmin, xmax, x2, ymin, ymax, y2, &
                            psivl, pds(1), irestart
 1010   FORMAT (' open contour in CNTOUR1'        / &
                ' xmin, xmax, x = ', 3(2x, e16.8) / &
                ' ymin, ymax, y = ', 3(2x, e16.8) / &
                ' psivl, pds(1) = ', 2(2x, e16.8), '    irestart =', i4)
      END IF
      IF (iauto .EQ. 1)  go to 1030
      ierr = 1
      IF (itty .NE. 0) THEN
          IF(myid == master)WRITE  (itty, 2001)  psiaxd, psivl, ipts
 2001     FORMAT (' subroutine CNTOUR1 has detected an error:' / &
                '   psiaxd, psivl, ipts =', 2(2x, 1pe14.6), 2x, i5)
          IF (itry .GT. 100) THEN
            IF(myid == master)WRITE  (itty, 2002) itry
 2002       FORMAT (/ &
            ' The given input value of psi on the plasma boundary' / &
            ' is so inappropriate that after ', i5, ' tries'       / &
            ' we still could not find a closed flux surface'       / &
            ' consistent with the mhdgrid and limiter'             / &
            ' geometry.  WE ARE GIVING UP.  Please fix your value' / &
            ' of psilim in the input eqdsk file and try again.'    /)
          END IF
 2003     FORMAT ( ' ERROR: x1, x2, psivl =', 3(2x, 1pe14.6))
 2004     FORMAT ( ' ERROR: y1, y2, psivl =', 3(2x, 1pe14.6))
 2005     FORMAT ( ' ERROR: psi2, psivl =',   2(2x, 1pe14.6))
 2006     FORMAT ( ' ERROR: netwi, nitr8 =',  2(2x, i5))
 2007     FORMAT ( ' ERROR: cmult, xmult, ymult, xsave, ysave, psisave', &
                   ' nsteps =' / 6(2x, 1pe16.8), i5)
 2008     FORMAT ( ' ERROR: x1, y1, psi1, x2, y2, psi2, isgn = ' / &
                     6(2x,1pe16.8),i5)
          IF(myid == master)THEN
             IF (   x1 .EQ. x2   )  WRITE (itty, 2003) x1,x2,psivl
             IF (   y1 .EQ. y2   )  WRITE (itty, 2004) y1,y2,psivl
             IF ( psi2 .GT. psivl)  WRITE (itty, 2005) psi2 , psivl
             IF (newti .GT. nitr8)  WRITE (itty, 2006) newti, nitr8
             IF (cmult .LT. 0.99 )  WRITE (itty, 2007) cmult, xmult, ymult, &
                                        xsave, ysave, psisave, nsteps
             IF (cmult .LT. 0.99 )  WRITE (itty, 2008) x1, y1, psi1, x2, &
                                        y2, psi2, isgn
          ENDIF
      END IF
      RETURN
!
 1030 psivl0 = psivl
      psiout = psivl
      dapsi  = psiaxd - psivl0
      dmult  = 0.0005_DP
      IF (delta_psi .GT. 0.0) THEN
        maxtry=20
        dapsi=delta_psi
        dmult=1.0D0 / DFLOAT (maxtry)  ! increments psivl because of psi
      END IF
      psivl  = psivl0 + dapsi * dmult
      IF ((itry .GT. 0) .AND. (psiin .NE. 0.0) .AND. (iconvg .EQ. 1)) &
        psivl = 0.5 * (psiin + psiout)
!***  if (iconvg .eq. 0)  psivl = pds(1)  ! HSJ 8/28/96
      iautoc = 1
      IF (itty .NE. 0) THEN
        IF(myid == master)WRITE  (itty, 1020)  psivl0, psivl
 1020   FORMAT (' will try to correct improper setting of psilim' / &
                ' changing psilim from', e16.8,'  to  ',e16.8)
      END IF
      itry = itry + 1
      IF (itry .GT. 100) THEN
        iauto = 0           ! allows error exit return to be taken
        go to 1000          ! can't do it, so quit
      END IF
      go to 5               ! go back and try again
!
 1500 cmult = 0.5           ! set cmult < 1 so logic above kicks in
      x1    = xsave
      y1    = ysave
      psi1  = psisave
      isgn  = isgnsave      ! restore search direction
      go to 40
!
      END  SUBROUTINE cntour1

      SUBROUTINE cntour2 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin, &
                          yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin, &
                          xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y, &
                          nh,cspln,n2cspln,nh2,itty,iptsm,ierr,bpmag, &
                          iconvg,delta_psi)
! ----------------------------------------------------------------------
! this routine uses the bicubic spline representation of psi to get
! psi values at an arbitrary point on the MHD grid. the spline
! coefficients must exist (in cspln) before this subroutine is called.
! given (xaxd,yaxd) and a psi value, psivl, generate a contour of
! ordered points,(xc(i),yc(i)), i = 1,ipts, which has (xaxd,yaxd) as
! as an interior point. the contour must fully encircle (xaxd,yaxd).
! the search is limited to a rectangle specified by (xmin,xmax,ymin
! ymax).  dx and dy determine the basic increment in x and y for
! the coarse grid search.  they should be picked so that over the
! distance dx or dy, psi is single valued.  a sliding interval
! search rather than a binary search is used for the coarse grid so
! that non-monotonic psi can be handled.  dx and dy are in meters.
! dx = dy = 0.01 is suggested.  dang is angular step for the rays which
! emmanate from (xaxd,yaxd) in degrees.  dang should be set
! according to shape of equilibrium.  if dang is too small, many
! unnecessary points will be generated in (xc,yc).  if dang is too
! large for bperr(see below), the routine wastes computation time
! reducing the magnitude of dang.  dang = 10 deg. near psilim and
! dang = 30 deg. near psimax is suggested.  for highly elongated
! plasmas control with dang is difficult to set for all contours.
! therefore arcl may be used to set dang internally.  arcl is the
! arclength in meters taken from the current point to get to the
! next point.  set arcl to a large number (e.g. 10 meters) to
! overide this option and use dang only.  the angle increment will
! be the minimum of (arcl/rad,dang).  bperr is relative change in
! poloidal b field between (xc(i),yc(i)) and (xc(i+1),yc(i+1)).
! if the number of elements (xc,yc) exceeds the limit
! bperr is relaxed (up to twice) by a constant increment
! dbperr before an error exit is taken.
! bperr = 0.03 is suggested.  note: if dang yields bperr(computed)<
! bperr(input), then dang is used.  otherwise dang is successively
! reduced by 0.5 until this condition is meet.  (xemin,xemax,yemin,
! yemax) are min and max values of contour.  (yxmin,yxmax), are y
! values at x = xemin and x = xemax respectively.  (xymin,xymax) are x
! values at y = yemin and y = yemax respectively.  iauto is used as an
! "error recovery switch".  that is, if a flux surface is found to
! pass outside the search box an error exit is taken by way of label
! 1000, if iauto = 0.  if iauto = 1, the input value of psi, psivl, is
! increased toward the value of psi at (xaxd,yaxd) until a closed
! surface inside the search box is found.  obviously this option
! makes sense only for the limiter flux surface.  it could be used
! to find the plasma boundary if an appopriate search box were
! defined.  however its primary use here is to slightly modify
! psilim so that a closed limiter flux surface is found.  the option
! is necessary due to the fact some codes generate eqdsks
! using linear interpolation, which differs from the bicubic inter-
! polation used here.  the user is informed that psivl was changed
! by returning iautoc = 1.  if no change occured iautoc = 0.
! if no error ocurred ierr = 0 is returned. if an error ocurred then
! ierr = 1 is returned and the results of this calculation are not valid.
!
! --- input (see description above)
!  xaxd
!  yaxd
!  psivl   psi value for which contour is desired
!  dang
!  arcl   !!!!!!!!!!!!!!!!!!!no longer used
!  bperr   dang,arcl,bperr control point density on contour returned
!  dx
!  dy      defines step length taken
!          input dx and dy as 0.0 to set step length in this subroutine
!  xmin
!  xmax
!  ymin
!  ymax    define the search box
!             NOTE: the search box could be defined as the entire
!                   MHD grid (i.e., xmin = x(1),xmax = x(nw),ymin = y(1),
!                   ymax = y(nh)) However this will normally not work
!                   satisfactorily because there are local minima and
!                   maxima around the f coils which would confuse the
!                   contour tracing and require many more calculations.
!                   (these additional calculations are not done here!)
!                   The best way to define the search box is to set it
!                   equal to a rectangle,within which the plasma must
!                   reside under all circumstances. The only such
!                   rectangle to suggest itself is the one defined by
!                   the extremes of the limiter. Hence it is suggested
!                   that xmin be set to the min x value of the limiter
!                   points vector,xmax the max x value of limiter points
!                   vector and similarly for ymin,ymax.
!  iaouto  = 1 automatic error recovery, = 0 no recovery
!             (iauto may be set to 0 in this subroutine
!              if a failure occurs)
!  x(1..nw)
!  y(1..nh)   the MHD grid vectors
!  cspln(n2cspln,nw,nh2)  bicubic spline coefficient array
!  itty      iounit for diagnostic messages. set to 0 to suppress.
!  iptsm     max allowed number of points to be returned in
!            vectors xc,yc
!
! delta_psi        this is the psi grid spacing. Define delta_psi as
!                  psi(penultimate value)-psi(last(i.e., edge) value)
!                  If the value of psi to be traced leads to an error
!                  and iauto=1 then delta_psi is used to control the
!                  amount by which the value of psi is adjusted. Otherwise
!                  delta_psi is not used,so you set its value accordingly.
!
! --- output:
!  xemin
!  xemax
!  yemin 
!  yemax
!  yxmin
!  yxmax
!  xymax
!  xymin    limits of contour (see above)
!  xc(1,..ipts)
!  yc(1,..ipts)   the ipts contour points found
!  bpmag(1,..ipts)    magnitude of poloidal bfield at xc,yc
!  ipts
!  ierr   error flag, = 0 no error, = 1 error, output not usable.
!             if ierr .gt. 1 then ierr equals the number of points found
!             on the contour up to this point. Since this number will
!             be exceeded if the code continues,the subroutine must
!             stop. The error condition on bperr has allready been relaxed
!             before this condition causes the subroutine to terminate.
!             Hence in order to continue this subroutine could be called
!             again with a larger value of arcl. (Note: in practice we have
!             not found this to be a problem,provided that iptsm is large
!             enough to begin with toi give an reasonable representation
!             of the contour).
! ------------------------------------------------------------------ HSJ
!

!


      IMPLICIT  NONE
!
      REAL(DP)   pds(6), xc(*), yc(*), x(nw), y(nh),             &
                 cspln(n2cspln,nw,nh2), bpmag(*)
      REAL(DP)   pi, piov2, piov4, fpiov4, spiov4 ,twopi,        &
                 tpiov4, tpiov2

      REAL(DP)   xemin ,xemax,yemin ,yemax ,xymin,xymax,         &
                 yxmin , yxmax ,psiaxd,xmin,xmax,ymin,ymax,x2,   &
                 y2,ymult,sint,cost,xns,yns,thetnew,psival,dang, &
                 bperr,dx,dy,psivl,delta_psi,bpmintol,psiin,     &
                 bperrsav,dbperr,dthetmin,dthet,dxx,dyy,dxmin,   &
                 dymin,serrt,derrt,x1,y1,psi1,psi2,cmult,        &
                 dpsi,ysave,arclmin,bp1,bp2,sintsq,costsq,       &
                 dtarcl,xsave,psisave,xmult,dpsids,psivl0,dapsi, &
                 serr,dmult,tstart,psivln,psiout,delx,dely,step, &
                 estmrad,dthet0,thet1


      REAL(DP)   xaxd, yaxd, thet, a, bincp,tend,xn,yn,arcl

      INTEGER(I4B) nitr8,nstepc,itry,irestart,ierr,iautoc,ilsave,&
                   jlsave,ier,nothing,nw,nh,n2cspln,nh2, iflg,   &
                   isgn,iauto,ipts,iptsm,iconvg,nsteps,j,        &
                   maxtry,newti,iqpt,isgnsave,itty

      LOGICAL    bptol, reverse

      DATA       pi, piov2, piov4, fpiov4, spiov4                 &
       /3.141592654,1.570796327,0.7853981634,3.926990818,5.497787145/
      DATA       twopi,tpiov4, tpiov2, nitr8                      &
       /6.283185307,2.356194491,4.712388904,20/



      itry     = 0
      irestart = 0
      psiin    = zeroc
      ierr     = 0
      iautoc   = 0
      bperrsav = bperr
      dbperr   = bperr / 2.0_DP
!      iprintit = 77
!      iprintito = 0
!      iptsmon  = 64
!      miunit = 0 !avoids dec f90 complaint
!      IF (iprintit .EQ. 1) THEN
!          iprintito = 1
!LLL: had been:
!     open (unit = iprintit, file = 'cntourdebug.dat',
!    .    status = 'UNKNOWN')
!LLL: note write's to miunit below, but miunit is undefined. propose:
!          CALL getioun(miunit,77)
!          OPEN (unit = miunit, file = 'cntourdebug.dat', &
!          status = 'UNKNOWN')
!      END IF
!
! --- get psi at (xaxd,yaxd)
!
!      CALL dbcevl1 (x,nw,y,nh,cspln,nw,xaxd,yaxd,pds,ier,1)
      CALL my_dbcevl1 (x,nw,y,nh,cspln,nw,xaxd,yaxd,pds,ier,1)
      IF (ier .NE. 0) THEN
        IF(myid == master)     &
           WRITE (itty, 4) ier
    4      FORMAT (/ ' subroutine CNTOUR2 detects an error:'   / &
                  ' EVAL_BICUBIC_SPLINE returned ier =', i5 / &
                  ' while getting psi on axis'               )
        lerrno = iomaxerr + 208_I4B
        CALL terminate(lerrno,nlog)
      END IF

      psiaxd = pds(1)
      IF (psiaxd .LT. psivl) THEN
        IF(myid == master)                                        &
        WRITE (itty, 3)  psiaxd, psivl
    3   FORMAT (/ ' subroutine CNTOUR2 detects an ERROR:'       / &
                  ' psi on magnetic axis    ', 1pe12.6          / &
                  ' psi contour to be found ', 1pe12.6          / &
                  ' must have max psi on axis, non-recoverable' / &
                  ' execution must stop'                         )
        lerrno = iomaxerr + 209_I4B
        CALL terminate(lerrno,nlog)
      END IF
    5 xemin =  1.0e10
      xemax = -1.0e10
      yemin =  1.0e10
      yemax = -1.0e10
      xymin =  0.0
      xymax =  0.0
      yxmin =  0.0
      yxmax =  0.0
      ipts  =  0
      irestart = irestart + 1
      tstart   = fpiov4
      tend     = tstart + twopi
      thet     = tstart
      dthet0   = twopi * dang / 360.0
      dthetmin = 1.0e-3    ! min increment in theta from one ray to next
      IF (dthetmin .GT. dthet0)  dthetmin = 0.5 * dthet0
      dthet    = 0.0
      dxx      = dx
      dyy      = dy
      IF (dxx .EQ. 0.0)  dxx = x(3)-x(1)    ! 2 times grid spacing is ok
      IF (dyy .EQ. 0.0)  dyy = y(3)-y(1)
      dxmin = 0.05*(x(2)-x(1)) ! subroutine stops if dxx becomes < dxmin
      dx    = dxx
      dy    = dyy
      IF (dxmin .GE. dxx)  dxmin = 0.5 * dxx
      serrt = 3.5e-06
      derrt = 0.5e-07
!
! --- serrt is absolute error convergence criteria for Newton's method below
! --- loop over theta from tstart to tend( = tstart+twopi)
!
   10 thet   = thet + dthet
      IF (thet .GT. twopi .AND. tstart .NE. 0.0) THEN
        thet = thet - twopi
        tend = tend - twopi
      END IF
      thet = MIN (thet, tend)
      IF (thet .EQ. tend)  go to 200
      dxx = dx
      dyy = dy
!
! --- get equation of ray emanating from (xaxd, yaxd)
!
      IF (( piov4 .LE. thet) .AND. (thet .LE. tpiov4))  go to 20
      IF ((fpiov4 .LE. thet) .AND. (thet .LE. spiov4))  go to 20
!
! --- y as a function of x            y = a*x+bincp
!
      isgn = -1
      IF ((thet .LT. piov4) .OR. (thet .GT. spiov4))  isgn = 1
      a = TAN (thet)
      iflg = 0
      bincp = yaxd-a*xaxd
      go to 30
!
! --- x as a function of y            x = a*y+bincp
!
   20 isgn = 1
      IF (thet .GT. pi)  isgn = -1
      IF (isgn .EQ. -1)  go to 22
      thet1 = piov2 - thet
      IF (thet .GT. piov2)  thet1 = twopi - ABS (thet1)
      go to 25
   22 thet1 = tpiov2 - thet
      IF (thet .GT. tpiov2)  thet1 = pi - ABS (thet1)
   25 a     = TAN (thet1)
      iflg  = 1
      bincp = xaxd - a*yaxd
   30 CONTINUE
!
! --- now have y = a*x+bincp    (iflg=0)   or
! --- x = a*y+bincp             (iflg=1)
!
      x1   = xaxd
      y1   = yaxd
      cost = COS (thet)
      sint = SIN (thet)
      psi1 = psiaxd
!
! --- sliding interval search. max width of interval ~1.41*(dx or dy)
!
   40 IF (iflg .EQ. 1)  go to 50
!
! --- search in x
!
      xmult = 1.0
   41 x2 = x1+isgn*dxx*xmult         ! increment x
!      IF (iprintit .EQ. 1 .AND. ipts .EQ. iptsmon) &
!          WRITE (miunit, 3001) psivl,x1,y1,x2,y2,xmin,xmax,ymin,ymax
! 3001 FORMAT (2x,"psivlx,x1,y1,x2,y2,xmin,xmax,ymin,ymax  =" / &
!              5(2x,1pe16.8),4(2x,1pe16.8))
      IF (x2 .EQ. x1)  go to 1000    ! psivl is not within reach
      x2 = MAX (x2,xmin)
      x2 = MIN (x2,xmax)             ! limit the search
      IF (x2 .EQ. x1)  go to 1000    ! psivl is not within reach
      y2 = a*x2+bincp                ! get corresponding value of y
      IF (y2 .LT. ymin .OR. y2 .GT. ymax) THEN
        xmult = 0.5 * xmult
        go to 41
      END IF
      go to 60
!
! --- search in y
!
   50 ymult = 1.0
   51 y2 = y1+isgn*dyy*ymult          ! increment y
!      IF (iprintit .EQ. 1 .AND. ipts .EQ. iptsmon) &
!              WRITE(miunit,3000)psivl,x1,y1,x2,y2, &
!                                  xmin,xmax,ymin,ymax
!
 3000         FORMAT(2x,"psivly,x1,y1,x2,y2,xmin,xmax,ymin,ymax  =",/, &
                     5(2x,1pe16.8),4(2x,1pe16.8))
      IF (y2 .EQ. y1)  go to 1000     ! psivl is not within reach
      y2 = MAX (y2,ymin)
      y2 = MIN (y2,ymax)              ! limit the search
      IF (y2 .EQ. y1)  go to 1000     ! psivl is not within reach
      x2 = a*y2+bincp                 ! get corresponding value of x
      IF (x2 .LT. xmin .OR. x2 .GT. xmax) THEN
          ymult = ymult * 0.5
          go to 51
      END IF
!
!   60 CALL dbcevl1 (x,nw,y,nh,cspln,nw,x2,y2,pds,ier,1)
   60 CALL my_dbcevl1 (x,nw,y,nh,cspln,nw,x2,y2,pds,ier,1)
      IF (ier .NE. 0) THEN
        IF(myid == master)THEN
           WRITE (itty, *)  'CNTOUR2 reports, call to DBCEVL1 ,ier =', ier
           WRITE (itty, *)  'arguments x2,y2 =',x2,y2
           WRITE (itty, *)  'xgrid(1),xgrid(nw) =',x(1),x(nw)
           WRITE (itty, *)  'ygrid(1),ygrid(nh) =',y(1),y(nh)
        ENDIF
!        CALL STOP ('subroutine CNTOUR2: unspecified problem', 162)
        lerrno = 217 + iomaxerr
        call terminate(lerrno,nlog)
      END IF
      psi2 = pds(1)
      dpsi = (psivl-psi1)*(psivl-psi2)
!      IF (iprintit .EQ. 1 .AND. ipts .EQ. iptsmon) &
!              WRITE(miunit,3002)dpsi,psivl,psi1,psi2
! 3002 FORMAT(2x,"dpsi,psivl,psi1,psi2 =",/, &
!       4(2x,1pe16.8))
      IF (dpsi .LE. zeroc)  go to 70
      x1   = x2
      y1   = y2
      psi1 = psi2
      go to 40
!
! --- now have psivl between psi1 and psi2,converge using newton-raphson
!
   70 newti = 0
      IF (iflg .EQ. 1)  go to 75
      xn = x1+isgn*dxx * 0.5_DP
      yn = a*xn+bincp
      go to 80
   75 yn = y1+isgn*dyy * 0.5_DP
      xn = a*yn+bincp
!   80 CALL dbcevl1 (x,nw,y,nh,cspln,nw,xn,yn,pds,ier,3)
   80 CALL my_dbcevl1 (x,nw,y,nh,cspln,nw,xn,yn,pds,ier,3)
      IF (ier .NE. 0) THEN
        newti = nitr8
        go to 76
      END IF
      dpsids = pds(2)*cost+pds(3)*sint
      dpsi   = pds(1)-psivl
      serr   = -dpsi/dpsids
      IF (ABS (serr) .LT. serrt)  go to 90
      IF (psivl .EQ. zeroc.AND. ABS (dpsi) .LT. derrt)  go to 90
      IF (psivl .NE. zeroc) THEN
          IF (ABS (dpsi/psivl) .LT. derrt)  go to 90
      END IF
      delx  = serr*cost
      dely  = serr*sint
      xn    = xn+delx
      yn    = yn+dely
      newti = newti+1
   76 IF (newti .GE. nitr8) THEN
         IF (dxx .LE. dxmin)THEN         
           IF(myid == master)     &
              WRITE(ITTY,FMT='("subroutine CNTOUR2: newti problem,dxx .LE. dxmin")')dxx,dxmin
           lerrno = iomaxerr + 218_I4B
           CALL terminate(lerrno,nlog)
         ENDIF
        dxx = 0.5_DP * dxx
        dyy = 0.5_DP * dyy
!      IF (iprintit .EQ. 1 .AND. ipts .EQ. iptsmon) &
!              WRITE(miunit,3004)
! 3004 FORMAT(2x,"returning to label 40 with dxx,dyy = ",/, &
!             2(2x,1pe16.8))
        go to 40        ! 6/11/93
      END IF
!
      go to 80
!
! --- end of newton iteration
! --- check for sufficient accuracy in point spacing as determined by thet
! --- accuracy test is based on a relative error in poloidal b field of bperr
!
   90 bp2 = SQRT (pds(2)**2 + pds(3)**2) / xn
      IF (thet .EQ. tstart)  go to 100    ! first point of contour
!
! --- if relative change in bp is less than bperr store the point just found:
!
!      IF (iprintit .EQ. 1 .AND. ipts .EQ. iptsmon) &
!             WRITE(miunit,3005)bp2,bp1,bperr,dthetmin,dthet
! 3005 FORMAT(2x,"bp2,bp1,bperr,dthetmin,dthet =",/, &
!      5(2x,1pe16.6))
      IF (ABS (bp2-bp1) / MAX (bp2,bp1) .LT. bperr)  go to 100
!
! --- if dthet is currently at its smallest allowed value of dthetmin
! --- and relative change in bp is still greater than bperr continue anyway:
!
      IF (dthet .EQ. dthetmin)  go to 100
!
! --- spacing too large for grad psi. decrease theta and try again
!
      thet  = thet - dthet          ! reset to previous value
      dthet = dthet * 0.5             ! take half the angle increment
      dthet = MAX (dthet, dthetmin) ! dthet lower limit must be observed
      go to 10        ! go back and try with a new ray closer to old one
  100 bp1   = bp2
      ipts  = ipts+1
      IF (ipts .LE. iptsm-1)  go to 150
!
! --- we ran out of storage for the contour points. assuming iptsm is
! --- set to a reasonable value this should not happen unless bperr is
! --- set unrealistically small. we increase bperr by dbperr and try
! --- to generate the contour again with the relaxed error requirement.
!
      IF (bperr-bperrsav .GE. 2.0*dbperr)  go to 110
      bperr = bperr+dbperr
      go to 5
  110 IF (itty .NE. 0) THEN
        IF(myid == master)WRITE (itty, 1040)  iptsm, dthet0, bperr, psivl
 1040   FORMAT &
       ('  more points on plasma '                                     / &
        '  boundary were found than dimension of xcontr,ycontr allows' / &
        '  allowed number = ',i5                                       / &
        '  adjustment in point spacing will be made', &
        '  dthet0,bperr,psivl =',5(2x,1pe12.4))
      END IF
!
      ierr = ipts
      RETURN
!
  150 xc(ipts) = xn   ! found new point; save it and set up for next ray
      yc(ipts) = yn
      bpmag(ipts) = bp2
      xemin = MIN (xemin, xn)
      xemax = MAX (xemax, xn)
      IF (xemax .EQ. xn)  yxmax = yn
      IF (xemin .EQ. xn)  yxmin = yn
      yemin = MIN (yemin, yn)
      yemax = MAX (yemax, yn)
      IF (yemax .EQ. yn)  xymax = xn
      IF (yemin .EQ. yn)  xymin = xn
      IF (ABS (bp2-bp1) / MAX (bp2,bp1) .LE. 0.25 * bperr) &
      dthet = 2.0 * dthet
      dthet = MIN (dthet, dthet0)        ! don't let dthet get too large
      dthet = MAX (dthet, dthetmin)      ! don't let dthet get too small
      go to 10                           ! continue with the next ray
!
! --- done with this contour. close contour
!
  200 ipts        = ipts+1
      xc(ipts)    = xc(1)
      yc(ipts)    = yc(1)
      bpmag(ipts) = bpmag(1)
      bperr       = bperrsav
!
      IF (itry .NE. 0 .AND. iconvg .EQ. 1) THEN
        psiin  = psivl
        psivln = 0.5 * (psiout+psiin)
        IF (ABS (psivln-psivl) .LT. 1.0e-6) THEN
          IF (itty .NE. 0) THEN
            IF(myid == master)WRITE  (itty, 1045) psivln
 1045       FORMAT (' final converged value of psi = ', 1pe16.8)
          END IF
          RETURN
        ELSE
          psivl = psivln
          go to 5
        END IF
      END IF
      RETURN
!
! --- errors
!
! 1000 CALL dbcevl1 (x, nw, y, nh, cspln, nw, x2, y2, pds, ier, 1)
 1000 CALL my_dbcevl1 (x, nw, y, nh, cspln, nw, x2, y2, pds, ier, 1)
!      IF (iprintito .EQ. 1) THEN
!        iprintit = miunit
!        CALL giveupus(iprintit)
!        CLOSE (unit = iprintit)
!      END IF
      IF (itty .NE. 0) THEN
        IF(myid == master)WRITE  (itty, 1010) xmin, xmax, x2, ymin, ymax, y2, &
                            psivl, pds(1), irestart
 1010   FORMAT (' open contour in CNTOUR2'       / &
                ' xmin, xmax, x = ', 3(2x, e16.8) / &
                ' ymin, ymax, y = ', 3(2x, e16.8) / &
                ' psivl, pds(1) = ', 2(2x, e16.8), '    irestart =', i4)
      END IF
      IF (iauto .EQ. 1)  go to 1030
      ierr   = 1
      RETURN
 1030 psivl0 = psivl
      psiout = psivl
      dapsi  = psiaxd - psivl0
      dmult  = 0.0005
      IF (delta_psi .GT. 0.0) THEN
        maxtry=99
        dapsi=delta_psi
        dmult = 1.0D0 / DFLOAT (maxtry)! increments psivl because of psi
      END IF
      psivl = psivl0 + dapsi * dmult
      IF ((itry .GT. 0) .AND. (psiin .NE. 0.0) .AND. (iconvg .EQ. 1)) &
        psivl = 0.5 * (psiin+psiout)
!***  if (iconvg .eq. 0)  psivl = pds(1)  ! HSJ 8/28/96
      iautoc = 1
      IF (itty .NE. 0 .AND. myid == master) &
               WRITE  (itty, 1020)  psivl0, psivl
 1020 FORMAT (' boundary search, will change psilim from' / &
                e16.8, '  to  ', e16.8, '  and try again')
      itry = itry + 1
      IF (itry .GT. 50) THEN
        iauto = 0
        go to 1000        ! can't do it, quit
      END IF
      go to 5             ! go back and try again
!
      END SUBROUTINE cntour2

      SUBROUTINE findeqlin (xaxd, yaxd, thet, a, bincp, iflg, isgn)
!

      USE nrtype,                                ONLY : Dp,I4B

      IMPLICIT  NONE

      REAL(DP) pi,piov2,piov4,fpiov4,spiov4,twopi, tpiov4, tpiov2,   &
               a,bincp,thet,thet1,xaxd,yaxd

      INTEGER(I4B) isgn,iflg
!
! ----------------------------------------------------------------------
! --- find eq(uation)lin(e) finds the equation of the straight line
! --- passing through (xaxd,yaxd) and with slope whose angle is thet radians.
! --- input
!
!   xaxd
!   yaxd
!   thet
!
! --- output
!
!   a
!   bincp
!   iflg  =0 or 1 depending on how the line is defined, see below
!   isgn  +1 or -1,gives direction in which to move along the line to go
!           from (xaxd,yaxd) toward the boundary in direction thet
!
! ----------------------------------------------------------------------
!
      DATA pi,piov2,piov4,fpiov4,spiov4 &
          /3.141592654,1.570796327,0.7853981634,3.926990818,5.497787145/
      DATA twopi, tpiov4, tpiov2 &
          /6.283185307,2.356194491,4.712388904/
!
! --- get equation of ray emanating from (xaxd,yaxd)
!
      IF (( piov4 .LE. thet) .AND. (thet .LE. tpiov4))  go to 20
      IF ((fpiov4 .LE. thet) .AND. (thet .LE. spiov4))  go to 20
!
! --- y as a function of x            y = a*x+bincp
!
      isgn = -1
      IF ((thet .LT. piov4) .OR. (thet .GT. spiov4))  isgn = 1
      a = TAN (thet)
      iflg = 0
      bincp = yaxd-a*xaxd
      RETURN
!
! --- x as a function of y            x = a*y+bincp
!
   20 isgn = 1
      IF (thet .GT. pi)  isgn = -1
      IF (isgn .EQ. -1)  go to 22
      thet1 = piov2 - thet
      IF (thet .GT. piov2)  thet1 = twopi - ABS (thet1)
      go to 25
   22 thet1 = tpiov2-thet
      IF (thet .GT. tpiov2)  thet1 = pi - ABS (thet1)
   25 a     = TAN (thet1)
      iflg  = 1
      bincp = xaxd-a*yaxd
      RETURN
!
! --- now have y = a*x+bincp    (iflg=0)   or
! --- x = a*y+bincp             (iflg=1)
!
      END SUBROUTINE findeqlin

      SUBROUTINE fixedcntour (rplasbdry,zplasbdry,nplasbdry, &
                              rcontr,zcontr,ncontr, &
                              rcmin,rcmax,zcmin,zcmax, &
                              rzcmin,rzcmax,zrcmin,zrcmax, &
                              rmhdgrid,zmhdgrid,nw,nh, &
                              bpcontr,cspln,n2cspln,nh2,pds)
! ----------------------------------------------------------------------
! --- SUBROUTINE handles the fixed boundary plasma contour CASE
! ------------------------------------------------------------------ HSJ

      USE nrtype,                       ONLY : DP,I4B

      IMPLICIT  NONE

      INTEGER(I4B) nw,nh,ncontr,nplasbdry,nh2,n2cspln,j,ier

      REAL(DP) rplasbdry(*),zplasbdry(*),rcontr(*),zcontr(*),     &
               pds(*),cspln(n2cspln,nw,nh2),rmhdgrid(nw),         &
               zmhdgrid(nh),bpcontr(*),rcmin,rcmax,zcmin,zcmax,   &
               rzcmin,rzcmax,zrcmin,zrcmax



          rcmin = rplasbdry(1)
          rcmax = rplasbdry(1)
          zcmin = zplasbdry(1)
          zcmax = zplasbdry(1)
          DO j=1,nplasbdry
              rcontr(j) = rplasbdry(j)
              zcontr(j) = zplasbdry(j)
              rcmin = MIN (rcmin,rcontr(j))
              IF (rcmin .EQ. rcontr(j))  zrcmin = zcontr(j)
              rcmax = MAX (rcmax,rcontr(j))
              IF (rcmax .EQ. rcontr(j))  zrcmax = zcontr(j)
              zcmin = MIN (zcmin,zcontr(j))
              IF (zcmin .EQ. zcontr(j))  rzcmin = rcontr(j)
              zcmax = MAX (zcmax,zcontr(j))
              IF (zcmax .EQ. zcontr(j))  rzcmax = rcontr(j)
!              CALL dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,nw, &
!                            rcontr(j),zcontr(j),pds,ier,3)
              CALL my_dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,nw, &
                            rcontr(j),zcontr(j),pds,ier,3)
              bpcontr(j) = (1.0_DP/rcontr(j)) * SQRT (pds(2)**2+pds(3)**2)
          END DO
          ncontr = nplasbdry

      RETURN

      END SUBROUTINE fixedcntour


      SUBROUTINE getrmaj (zelev, nf, ierr, rmajor,    &
                          psirmaj, isigncur, psiedge,cspln,n2cspln,nw,nh,nh2,pds)
! ------------------------------------------------------------------ HSJ
! --- at elevation zelev,get a uniformly-spaced major radius array, rmajor,
! --- starting at the inside edge of the plasma and working to the
! --- outside edge. for each value of rmajor(j),get a corresponding
! --- value of psi,stored in psirmaj(j).
! --- ROUTINE ASSUMES BICUBIC SPLINE ARRAY CSPLN IS SET UP AND READY FOR USE
! --- BY DBCEVL1 !
! ------------------------------------------------------------------ HSJ
      USE nrtype,                                   ONLY : DP,I4B

      USE tport_mhd_grid_data,                      ONLY : rmhdgrid,zmhdgrid,xmagn1

      IMPLICIT NONE

      INTEGER(I4B) j,nf,isigncur,ierr,n2cspln,nw,nh,nh2,ind

      REAL(DP) rmajor(nf), psirmaj(nf),psiedge,zelev,         &
               toler,ra,rb,rg,pds(*),psirzp,dr

      REAL(DP) cspln(n2cspln,nw,nh2)

      ierr   = izero

      toler  = 1.0e-06_DP

! --- first get inside radius

      ra = rmhdgrid(1)
      rb = xmagn1

  100 rg = 0.5_DP * (ra+rb)
      CALL chkinout (rg, zelev, psirzp, ind, isigncur, psiedge,    &
                     cspln,n2cspln,nw,nh,nh2,pds)
      IF (ABS (ra-rb) .LT. toler)  go to 120
      IF (ind .EQ. 0) THEN
        ra = rg
      ELSE
        rb = rg
      END IF
      go to 100

  120 rmajor(1) = rg

! --- next get outside radius

      ra = xmagn1
      rb = rmhdgrid(nw)

  140 rg = 0.5 * (ra+rb)
      CALL chkinout (rg, zelev, psirzp, ind, isigncur, psiedge,  &
           cspln,n2cspln,nw,nh,nh2,pds)
      IF (ABS (ra-rg) .LT. toler)  go to 130
      IF (ind .EQ. 0) THEN
        rb = rg
      ELSE
        ra = rg
      END IF
      go to 140

  130 rmajor(nf) = rg
      dr = (rmajor(nf)-rmajor(1))/(nf-1)
      DO 150 j=2,nf-1
        rmajor(j) = rmajor(j-1)+dr
!        CALL dbcevl1 (rmhdgrid, nw, zmhdgrid, nh, cspln, nw, rmajor(j), &
!                     zelev, pds, ierr, 1)
        CALL my_dbcevl1 (rmhdgrid, nw, zmhdgrid, nh, cspln, nw, rmajor(j), &
                     zelev, pds, ierr, 1)
  150 psirmaj(j)  = pds(1)
      psirmaj(1)  = psiedge
      psirmaj(nf) = psiedge

      RETURN

      END SUBROUTINE getrmaj


      SUBROUTINE limiter_check(rcminm,rcmaxm,zcminm,zcmaxm,xlimiter,    &
                         ylimiter,nlimiter)
!-------------------------------------------------------------------
!     check to make sure that limiter remains consistent with
!     fixed boundary calculations when the boundary is allowed
!     to move slightly ( limiter may have been set tightly up
!     against the plasma,not allowing for any variation)
!     Return the limiter points modified to accomodate the new
!     extremes in the boundary. Dont allow limiter to ,move off
!     the computational grid however. (such an adjustment will require
!     a new starting eqdsk).
!     All quantities are in meters  
!     
!---------------------------------------------------------HSJ-1-20-01
      USE nrtype,                                   ONLY : DP,I4B
      USE Plasma_properties ,                       ONLY : dischg
      USE tport_mhd_grid_data,                      ONLY : rmhdgrid,zmhdgrid

      IMPLICIT NONE
      REAL(DP) rcminm,rcmaxm,zcminm,zcmaxm
      REAL(DP) ylimiter(*),xlimiter(*)
      REAL(DP) xdim, ydim, ymid, redge
      INTEGER(I4B) nlimiter,nw,nh

      nw = dischg%nr_mhd    ; nh = dischg%nz_mhd


      xlimiter(nlimiter+1) = MIN(xlimiter(nlimiter+1),rcminm)
!      xlimiter(nlimiter+1) = MAX(xlimiter(nlimiter+1),rmhdgrid(1))

      xlimiter(nlimiter+2) = MAX(xlimiter(nlimiter+2),rcmaxm)
!      xlimiter(nlimiter+2) = MIN(xlimiter(nlimiter+2),rmhdgrid(nw))

      ylimiter(nlimiter+1) = MIN(ylimiter(nlimiter+1),zcminm)
!      ylimiter(nlimiter+1) = MAX(ylimiter(nlimiter+1),zmhdgrid(1))

      ylimiter(nlimiter+2) = MAX(ylimiter(nlimiter+2),zcmaxm)
!      ylimiter(nlimiter+2) = MIN(ylimiter(nlimiter+2),zmhdgrid(nh))

      RETURN 
      END  SUBROUTINE limiter_check


      SUBROUTINE new_point (pds, arcl, xn, yn, thet, twopi, spiov4, &
                            xmin, xmax, ymin, ymax, yaxd, xaxd, ier, &
                            sint, cost, xns, yns, thetnew)
! ----------------------------------------------------------------------
! --- subroutine used by CNTOUR to get approximation for next point on
! --- contour,given the current point and the gradient of psi .
! --- the tangent line at the point (xn,yn) is
!             df = 0=df/dx*dx+df/dy*dy
!   a circle of radius arcl,centered at (xn,yn) is
!             arcl**2 = dr**2+dz**2
!   this gives us two eq. in two unknowns, dr and dz.
!   the new point is xn +/- dx and yn +/- dy where the signs have to
!   be picked so that thet increases. Note the special treatment
!   required as thetnew crosses zero (outside this routine).
!
! ------------------------------------------------------------------ HSJ
!
      USE nrtype,                                  ONLY : Dp,I4B


      IMPLICIT  NONE

      INTEGER(I4b) ier

      REAL(DP)  pds(*),dpx,dpy, alpha,dy,dx, arcl, xn, yn, thet, twopi, &
                spiov4,xmin, xmax, ymin, ymax, yaxd, xaxd, sint,cost, xns,  &
                yns, thetnew,proj
!
      ier = 0
      dpx = pds(2)
      dpy = pds(3)
      IF (ABS (dpx) .GT. ABS (dpy)) THEN
        IF (dpx .EQ. zeroc) THEN
          ier = 1
          RETURN
        END IF
        alpha = dpy/dpx
        dy    = arcl / SQRT (1.0+alpha*alpha)
        dx    = -alpha * dy
      ELSE
        IF (dpy .EQ. zeroc) THEN
          ier = 1
          RETURN
        END IF
        alpha = dpx/dpy
        dx    = arcl / SQRT (1.0_DP + alpha*alpha)
        dy    = -alpha * dx
      END IF
!
! ----------------------------------------------------------------------
! --- the sign on dx,dy must be taken so that thet increases.
! --- a unit vector in the direction of increasing thet
! --- (i.e., thet counterclockwise) is (-SIN (thet),COS (thet))
! --- the displacement vector is (dx,dy). its projection on the above
! --- vector must be positive and equals -dx*SIN (thet)+dy*COS (thet)
! ----------------------------------------------------------------------
!
      proj = -dx*sint+dy*cost
      IF (proj .LT. zeroc) THEN
        dx = -dx
        dy = -dy
      END IF
      xns = xn+dx
      yns = yn+dy
      IF (xns .LT. xmin .OR. xns .GT. xmax .OR. &
          yns .LT. ymin .OR. yns .GT. ymax) THEN
        ier = 1
        RETURN
      END IF
      thetnew = ATAN2 (yns-yaxd, xns-xaxd)
      IF (thetnew .LT. 0.0) THEN
        thetnew = thetnew+twopi
      ELSE IF (thet .GT. spiov4) THEN
        thetnew = thetnew+twopi
      END IF
      RETURN
!
      END SUBROUTINE new_point


      SUBROUTINE volcalc (rcontr, zcontr, ncontr, rma,zma, volume, area)
! ----------------------------------------------------------------------
! get the plasma volume, using the contour given by
! rcontr(j),zcontr(j), j = 1,2..ncontr
! use  vol = 2 * pi * contour integral(r*z*dr)
! use area = contour integral(z*dr)
! where the contour is given in rcontr, zcontr
! -----------------------------------------------------------HSJ--------

      USE nrtype,                                     ONLY : DP,I4B

      IMPLICIT  NONE 

      REAL(DP)   rcontr(*), zcontr(*)
      REAL(DP) xym,area,xyma,zzm,volsum,dx,zzp,slope,             &
               rminzm,rmaxzm,volume,rma,zma,xyp,xypa

      INTEGER(I4B)  ncontr,j

 

      xym    = rcontr(1)*zcontr(1)
      area   = zeroc
      xyma   = zcontr(1)
      zzm    = zma-zcontr(1)
      volsum = zeroc
      DO j=2,ncontr
        xyp    = rcontr(j)*zcontr(j)
        xypa   = zcontr(j)
        dx     = rcontr(j)-rcontr(j-1)
        volsum = volsum+(xyp+xym)/2.0*dx
        area   = area+(xyma+xypa)/2.0*dx
        xym    = xyp
        xyma   = xypa
        zzp    = zma-zcontr(j)
        IF (zzp*zzm .LE. 0.0) THEN
          slope = (rcontr(j)-rcontr(j-1))/(zcontr(j)-zcontr(j-1))
          IF (rcontr(j) .LT. rma)  rminzm = rcontr(j) + zzp*slope
          IF (rcontr(j) .GT. rma)  rmaxzm = rcontr(j) + zzp*slope
        END IF
        zzm = zzp
      END DO
      volume = ABS (volsum)*twopi
      area   = ABS (area)
      RETURN

      END SUBROUTINE volcalc



     SUBROUTINE chkinout (rpt, zpt, psirzp, ind, isigncur, psiedge, &
                          cspln,n2cspln,nw,nh,nh2,pds)
! ----------------------------------------------------------------------
! set ind = 1 if (rpt,zpt) is inside plasma ...
! set ind = 0 if (rpt,zpt) is outside plasma ...
! return psirzp,value of psi at (rpt,zpt) if ind = 1
! ----------------------------------------------------------------------- HSJ
!
      USE nrtype,                                    ONLY : DP,I4B

      USE tport_mhd_grid_data,                                ONLY : rmhdgrid,zmhdgrid, &
                                                                     rplasmin,rplasmax, &
                                                                     zplasmin,zplasmax

      IMPLICIT  NONE
      REAL(DP)  psichk,psiedge,rpt,zpt,pds(*),psirzp
      REAL(DP) cspln(n2cspln,nw,nh2)
      INTEGER(I4B) ind,isigncur,n2cspln,nw,nh2,ier,nh

      psichk = psiedge
      ind    = 1
      IF (rpt .LT. rplasmin)  ind = 0
      IF (rpt .GT. rplasmax)  ind = 0
      IF (zpt .LT. zplasmin)  ind = 0
      IF (zpt .GT. zplasmax)  ind = 0
      IF (ind .EQ. 0       )  RETURN
!      CALL dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,nw,rpt,zpt,pds,ier,1)
      CALL my_dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,nw,rpt,zpt,pds,ier,1)
      psirzp = pds(1)
      IF (isigncur .GT. 0 .AND. psirzp .LT. psichk)  ind = 0
      IF (isigncur .LT. 0 .AND. psirzp .GT. psichk)  ind = 0

      RETURN

      END SUBROUTINE chkinout

  END MODULE contour
