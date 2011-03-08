      subroutine read_eqdsk(cnin,iflag,message)
!***********************************************************************
!READ_EQDSK reads an EQDSK file from EFIT
!  J. Kinsey 2/25/03
!Input:
!  cnin-input file name
!Output:
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!Available output thru common block, efitdat:
!  nx_xy-number of x points on psi(x,y) grid
!  ny_xy-number of y points on psi(x,y) grid
!  nr_r-number of radial points
!  a0-minor radius, half diameter of boundary flux surface (m)
!  r0-major radius = center of boundary flux suface (m)
!  bt0-toroidal field at r0 (T)
!  rmag-major radius of magnetic axis (m)
!  zmag-vertical position of magnetic axis (m)
!  bmag-axial magnetic field (T)
!  x_xy(i)-horizontal grid for 2-D poloidal flux (m)
!  y_xy(i)-vertical grid for 2-D poloidal flux (m)
!  psi_xy(i,j)-poloidal flux/(2*pi) on 2-D grid (Wb/rad)
!  psi_r(i)-poloidal flux/2*pi (Wb/rad)
!  rhop_r-normalized poloidal flux grid proportional to psi
!  rhot_r(i)-normalized toroidal flux grid proportional to (Phi)**0.5
!  rin_r(i)-major radius grid on inside of torus in axis plane (m)
!  rout_r(i)-major radius grid on outside of torus in axis plane (m)
!  frb_r(i)-R*B_t (m*T)
!  ffp_r(i)-frb_r*dfrb_r/dpsi_r (rad*T)
!  p_r(i)-kinetic pressure (N/m**2)
!  pp_r(i)-dp_r/dpsi_r (rad*N/m**2/Wb)
!  q_r(i)-safety factor
!  grho1_r(i)-<|grad(rhot_r)|>
!  grho2_r(i)-<|grad(rhot_r)|**2>
!  bpout_r(i)-poloidal field at rout_r(i) (T)
!  btout_r(i)-toroidal field at rout_r(i) (T)
!  n_bdry-number of points on plasma boundary
!  x_bdry(i)-major radius of boundary point (m)
!  y_bdry(i)-vertical position of boundary point (m)
!  n_lim-number of points on limiter
!  x_lim(i)-horizontal positions of limiter points (m)
!  y_lim(i)-vertical positions of limiter points (m)
!***********************************************************************
      implicit none
c
      include '../inc/efitdat.m'
c
c Declaration of input variables
      character*(*)  cnin
c
c... Radial grid-psi
      integer       kpsi
      parameter     (kpsi=65) 
      real          psival(kpsi), torflux(kpsi), rho(kpsi)
      real           aspline(kpsi), bspline(kpsi), cspline(kpsi),
     &               cs2spline(kpsi,3), dspline(kpsi), espline(kpsi),
     &               fspline(kpsi)
c... Two dimensional poloidal flux-xy
      integer        mxnx_xy, mxny_xy
      parameter     (mxnx_xy=130, mxny_xy=130)
      real           x_xy(mxnx_xy), y_xy(mxny_xy), 
     &               psi_xy(mxnx_xy,mxny_xy)
c
c Declaration of output variables
      character*(*)  message
      integer        iflag, iflag1, ier
c
c Declaration of internal variables
      integer        i, j, nin, ncrt
      integer        nr_e, nlimiter, iunfrm, nw
      real           zpi, btor, rmajor, bmag, dum
      real           eqdskvol, eqdskarea, dpsi, dx, dy
      real           tension, tmax, bpar(4)
c ... all these in efitdat.m:
c     real           gth_r(mxnr_r), gph_r(mxnr_r),
c    &               vp_r(mxnr_r), r2_r(mxnr_r), rm2_r(mxnr_r),
c    &               rm1_r(mxnr_r), fm_r(mxnr_r), grth_r(mxnr_r),
c    &               b2_r(mxnr_r), bm2_r(mxnr_r), gr2bm2_r(mxnr_r),
c    &               ftrap_r(mxnr_r), fhat_r(mxnr_r)
      real           r2_r(mxnr_r)
      real           work(nr_r)
c
      zpi = atan2(0.0,-1.0)
      ncrt = 6
      nin = 20
c
c Open the EQDSK file
      open(unit=nin,status='old',file=cnin,form='formatted')
c
c Read the EQDSK file
      read(nin,'(a52,2i4)') etitle, nx_xy, ny_xy
c Same number of points in the radial and horizontal grids
      nr_e=nx_xy
c Check whether data exceeds maximum dimensions
      if(nx_xy.gt.mxnx_xy) then
c       Horizontal grid points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:x grid dimension exceeded'
        goto 1000
      endif
      if(ny_xy.gt.mxny_xy) then
c       Vertical grid points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:y grid dimension exceeded'
        goto 1000
      endif
      if(nr_e.gt.mxnr_r) then
c       Radial grid points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:r grid dimension exceeded'
        goto 1000
      endif
c
c Read in EQDSK data
c
c  rdim-width of computational region (m)
c  zdim-height of computational region (m)
c  rmin-major radius at inside edge of computational region (m)
c  zmid-vertical center of computational region (m)
c  psimag-poloidal flux at magnetic axis/(2*pi) (Wb/rad)
c  psilim-poloidal flux at limiter/(2*pi) (Wb/rad)
c  torcur-toroidal current (A)
      read(nin,'(5e16.9)') rdim,zdim,r0,rmin,zmid
      read(nin,'(5e16.9)') rmag,zmag,psimag,psilim,bt0
      read(nin,'(5e16.9)') current
      read(nin,'(5e16.9)') dum1
c Profile data
c  The grid is assumed to be equally spaced in poloidal flux
c  Node 1 is the axis and node nr_r is the edge
      read(nin,'(5e16.9)') (f_e(i),i=1,nr_e)     ! fpsi
      read(nin,'(5e16.9)') (p_e(i),i=1,nr_e)     ! ppsi
      read(nin,'(5e16.9)') (ffp_e(i),i=1,nr_e)   ! ffprime
      read(nin,'(5e16.9)') (pp_e(i),i=1,nr_e)    ! pprime
      read(nin,'(5e16.9)') ((psi_xy(i,j),i=1,nx_xy),j=1,ny_xy)  ! psi
      read(nin,'(5e16.9)') (q_e(i),i=1,nr_e)     ! qpsi
c Boundary and limiter data
      read(nin,'(2i5)') n_bdry, n_lim
      if(n_bdry.gt.mxn_bdry) then
!       Boundary points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:bdry grid dim exceeded'
        goto 1000
      endif
      if(n_lim.gt.mxn_lim) then
!       Limiter points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:lim grid dim exceeded'
        goto 1000
      endif
      read(nin,'(5e16.9)') (x_bdry(i),y_bdry(i),i=1,n_bdry)
      read(nin,'(5e16.9)') (x_lim(i),y_lim(i),i=1,n_lim)
c
c from subroutine reqdsk in cray209 of ONETWO
c
      nlimiter=99
      if (nlimiter .gt. 0) then
        call volcalc (x_bdry,y_bdry,n_bdry,
     .                rmag,zmag,eqdskvol,eqdskarea)
        write (ncrt, '(2x, a, f12.6)')
     .      'volume determined from eqdsk boundary values = ', eqdskvol
        write (ncrt, '(2x, a, f12.6)')
     .      'area   determined from eqdsk boundary values = ', eqdskarea
      endif
c
c set up psival,with psival(1) = edge, psival(nxeqd) = axis
c
      iunfrm=1
      call psiset (nx_xy,psimag,psilim,psival,iunfrm)
c
c calculate the toroidal flux, rho grid
c
      nw = nx_xy
      btor = abs(bt0)
      rmajor = r0
      call trap2 (psival, q_e, torflux, nw)
      do i=1,nw
        torflux(i) = 2.0*zpi*(torflux(i)-torflux(nw))
        rho(i) = sqrt(torflux(i)/(zpi*btor))
        bp_e(i) = rho(i)*btor/(rmajor*q_e(i))
c       write(ncrt,'i3,3f12.6') i, rho(i), q_e(i), bp_e(i)
      enddo
c
c Derived quantities
c  0-D :
      bmag=f_e(1)/rmag
c  1-D grids - map from uniform in Psi to uniform in Psi^0.5
      dpsi=(psilim-psimag)/(nr_e-1)
      psi_e(1)=psimag
      rhop_e(1)=0.0
      psi_e(nr_e)=psilim
      rhop_e(nr_e)=1.0
      do i=2,nr_e-1
        psi_e(i)=psi_e(i-1)+dpsi 
        rhop_e(i)=(psi_e(i)-psi_e(1))/(psi_e(nr_e)-psi_e(1))
      enddo
      psi_r(1)=psimag
      rhop_r(1)=0.0
      psi_r(nr_r)=psilim
      rhop_r(nr_r)=1.0
      dpsi=psilim-psimag
      do i=2,nr_r-1
        psi_r(i)=psi_r(1)+dpsi*(float(i-1)/float(nr_r-1))**2.0
        rhop_r(i)=(psi_r(i)-psi_r(1))/(psi_r(nr_r)-psi_r(1))
      enddo
      call w_lin_interp(nr_e,rhop_e,f_e,nr_r,rhop_r,frb_r,iflag,message)
      if(iflag.ne.0) then
        i=LEN(message)
        message='READ_EFIT_EQDSK(1)/'//message(1:i-19)
        if(iflag.eq.1) goto 1000
      endif
      call w_lin_interp(nr_e,rhop_e,p_e,nr_r,rhop_r,p_r,iflag,message)
      if(iflag.ne.0) then
        i=LEN(message)
        message='READ_EFIT_EQDSK(2)/'//message(1:i-19)
        if(iflag.eq.1) goto 1000
      endif
      call w_lin_interp(nr_e,rhop_e,ffp_e,nr_r,rhop_r,ffp_r,iflag,
     &                  message)
      if(iflag.ne.0) then
        i=LEN(message)
        message='READ_EFIT_EQDSK(3)/'//message(1:i-19)
        if(iflag.eq.1) goto 1000
      endif
      call w_lin_interp(nr_e,rhop_e,pp_e,nr_r,rhop_r,pp_r,iflag,message)
      if(iflag.ne.0) then
        i=LEN(message)
        message='READ_EFIT_EQDSK(4)/'//message(1:i-19)
        if(iflag.eq.1) goto 1000
      endif
      call w_lin_interp(nr_e,rhop_e,q_e,nr_r,rhop_r,q_r,iflag,message)
      if(iflag.ne.0) then
        i=LEN(message)
        message='READ_EFIT_EQDSK(5)/'//message(1:i-19)
        if(iflag.eq.1) goto 1000
      endif
c  2-D grid
      x_xy(1)=rmin
      dx=rdim/(nx_xy-1)
      do i=2,nx_xy
        x_xy(i)=x_xy(i-1)+dx
      enddo
      y_xy(1)=zmid-0.5*zdim
      dy=zdim/(ny_xy-1)
      do i=2,ny_xy
        y_xy(i)=y_xy(i-1)+dy
      enddo
c
c call FLUXAV to get metrics
c
      call fluxav_efit(bmag,rmag,zmag,nx_xy,ny_xy,x_xy,y_xy,psi_xy,
     &                 n_lim,x_lim,y_lim,nr_r,psi_r,frb_r,ffp_r,q_r,
     &                 rin_r,rout_r,rhot_r,elong_r,gth_r,gph_r,vol_r,
     &                 vp_r,r2_r,rm2_r,phit_r,fm_r,grth_r,b2_r,bm2_r,
     &                 grho1_r,grho2_r,gr2bm2_r,ftrap_r,fhat_r,bpout_r,
     &                 btout_r,a0,rm1_r,iflag1,message)
c
      write(ncrt,'(a,f8.4)') 'R0  = ',r0
      write(ncrt,'(a,f8.4)') 'a   = ',a0
      write(ncrt,'(a,f8.4)') 'Bt0 = ',bt0
      write(ncrt,'(a,f8.4)') 'rho(a) = ',rho(1)
      write(ncrt,100)
      do j=1,nr_r
        rhor_r(j)=rhot_r(j)*SQRT(phit_r(nr_r)/(zpi*bmag))
        grho1_r(j)=grho1_r(j)/a0
        grho2_r(j)=grho2_r(j)/a0**2.0
        write(ncrt,'(i3,8f8.4,2f10.5)') j,rhot_r(j),rin_r(j),rout_r(j), 
     &                          grho1_r(j), grho2_r(j), vol_r(j), 
     &                          elong_r(j), q_r(j), psi_r(j), phit_r(j)
      enddo
c
      tension = 0.0    ! don't use tension option of tspline
      tmax    = 0.0    ! max allowed tension
      bpar(1) = 0.0    ! set boundary conditions on spline
      bpar(2) = 0.0
      bpar(3) = 0.0
      bpar(4) = 0.0
c     call tspline (psival,grho1npsi,nr_e,bpar,cs2spline,kpsi,
c    &              ier,tension,aspline,bspline,cspline,dspline,
c    &              espline,fspline,tmax,psi_e,work,nr_r)
c
c     do j=1,nr_e
c       write(*,200) j, psi_e(j)
c     enddo
c
  100 format(6x,'rho',5x,'Rin',5x,'Rout',3x,'grho1',3x,'grho2',
     &       4x,'vol',4x,'kappa',5x,'q',8x,'psi',6x,'phi')
  200 format(i4,1pe14.5,0pf14.5)
 1000 close(unit=nin)
      return
      end
!***********************************************************************
      subroutine volcalc (rcontr,zcontr,ncontr,rma,zma,volume,area)
c
      implicit  integer (i-n), real (a-h, o-z)
c
c ----------------------------------------------------------------------
c get the plasma volume, using the contour given by
c rcontr(j),zcontr(j), j = 1,2..ncontr
c use  vol = 2 * pi * contour integral(r*z*dr)
c use area = contour integral(z*dr)
c where the contour is given in rcontr, zcontr
c ----------------------------------------------------------------------
c
      dimension  rcontr(*), zcontr(*)
      data       twopi /6.28318530718 /
c
      xym    = rcontr(1)*zcontr(1)
      area   = 0.0
      xyma   = zcontr(1)
      zzm    = zma-zcontr(1)
      volsum = 0.0
      do j=2,ncontr
        xyp    = rcontr(j)*zcontr(j)
        xypa   = zcontr(j)
        dx     = rcontr(j)-rcontr(j-1)
        volsum = volsum+(xyp+xym)/2.0*dx
        area   = area+(xyma+xypa)/2.0*dx
        xym    = xyp
        xyma   = xypa
        zzp    = zma-zcontr(j)
        if (zzp*zzm .le. 0.0) then
          slope = (rcontr(j)-rcontr(j-1))/(zcontr(j)-zcontr(j-1))
          if (rcontr(j) .lt. rma)  rminzm = rcontr(j) + zzp*slope
          if (rcontr(j) .gt. rma)  rmaxzm = rcontr(j) + zzp*slope
        end if
        zzm = zzp
      end do
      volume = ABS (volsum)*twopi
      area   = ABS (area)
      return
c
      end
!***********************************************************************
      subroutine psiset (npsi, psiaxis, psibdry, psival, iunfrm)
c
      implicit  integer (i-n), real (a-h, o-z)
c
c ----------------------------------------------------------------------
c calculate a vector of psi values, psival. note that
c psival(1) = psi on plasma edge,psival(npsi) = psi at magnetic axis
c ----------------------------------------------------------------------
c
      dimension psival(*)
c
      if (iunfrm .eq. 0) then
        do j=2,npsi-1
          psival(j) = psiaxis-(psiaxis-psibdry)*((npsi-j)/(npsi-1.0))**2
        end do
      else
        dpsi = (psibdry-psiaxis)/(npsi-1)
        do j=2,npsi-1
          psival(j) = psibdry-(j-1)*dpsi
        end do
      end if
      psival(1)    = psibdry
      psival(npsi) = psiaxis
      return
c
      end
!***********************************************************************
      subroutine trap2 (x, y, yint, npts)
c
      implicit  integer (i-n), real (a-h, o-z)
c
c ----------------------------------------------------------------------
c calculates integral y*dx using trapezoidal rule
c ----------------------------------------------------------------------
c
      dimension  x(*), y(*), yint(*)
c
      yint(1) = 0.0
      do 10 i=2,npts
   10 yint(i) = (y(i)+y(i-1))*(x(i)-x(i-1))*0.5+yint(i-1)
      return
c
      end
!***********************************************************************
      subroutine tspline (x, y, nx, bpar, cs, ic, ier, t, a, b, c,
     .                    fpp, r, dx, tmax, rgrid, tspl, npts)
c
      implicit  integer (i-n), real (a-h, o-z)
c
c --- tspline calculates the second derivatives of the spline at the knots.
c --- these coefficients define the spline and are stored
c --- in array cs(i,j),i = 1,#knots,j=1,3
c --- subroutine EVTSPLN is used to evaluate the spline,
c
c --- input
c
c      x      vector,length nx,of knot locations (must be in ascending order)
c      y          function values at the knots
c      n          #points in x and y
c    rgrid(i)     i = 1,2...npts values of x at which evaluation is desired
c      bpar       array used to specifiy boundary conditions
c      ic         exact row dimension of matrix c
c      t          tension parameter
c
c      a,b,c,     temporary work vectors of length nx
c      fppr,dx
c
c --- output
c
c      cs         array of spline coefficients (see above)
c --- the value of the spline for x in the interval (x(i),x(i+1)) is
c ---           (the range on i is 1 to nx-1 )
c          f(x) = cs(i,1)*cs(i,2) * SINH (t*(x(i+1)-x))+
c               (y(i)*t*t-cs(i,1))*cs(i,3)*(x(i+1)-x)
c               +cs(i+1,1)*cs(i,2) * SINH (t*(x-x(i)))+
c               (y(i+1)*t*t-cs(i+1,1))*cs(j,3)*(x-x(i))
c       tmax   max allowed value of tension parameter
c       tspl(i)      i = 1,2..npts values of y at points  rgrid(i)
c
c ------------------------------------------------------------------ HSJ
c
      logical    reverseset
      dimension  a(*),b(*),c(*),fpp(*),r(*),dx(*),rgrid(*),tspl(*)
      dimension  x(*),y(*),bpar(*),cs(ic,3)
c
c --- check input
c
      reverseset = .false.
      dxmin = ABS (x(nx)-x(1))
      nm1 = nx-1
      ier = 0
      if (ic .lt. nx)  ier = 129
      if (nx .lt. 2)  ier = 130
      if (ier .ne. 0)  go to 1000
      do 15 j=1,nm1
        dxmin = MIN (dxmin, ABS (x(j+1)-x(j)))
        if (x(j) .lt. x(j+1))  go to 15
        ier = 131
        go to 1000
   15 continue
c
c --- some initialization
c
   16 bp1   = bpar(1)
      bp2   = bpar(2)
      bp3   = bpar(3)
      bp4   = bpar(4)
      tmax  = 30.0 / dxmin    ! max tension, avoids overflow
      t     = MIN (t, tmax)
      tsq   = t*t
      nm1   = nx-1
      dx(1) = x(2) - x(1)
c
c --- in case smtdx is incorrectly set
c
      if (t .eq. 0.0)  smtdx = 0.0
c
c --- set up vectors a,b,c
c --- vector a is lower diagonal
c --- vector b is diagonal
c --- vector c is upper diagonal
c
      if (t .eq. 0.0)  go to 50     ! zero tension case must be separate
      do 20 j=2,nm1
        dx(j) = x(j+1)-x(j)
        a(j) = (1.0/dx(j-1)-t / SINH (t*dx(j-1)))/tsq
        b(j) = t * COSH (t*dx(j-1)) / SINH (t*dx(j-1))
        b(j) = b(j)-1.0/dx(j-1)-1.0/dx(j)
        b(j) = (b(j)+t * COSH (t*dx(j)) / SINH (t*dx(j)))/tsq
        c(j) = (1.0/dx(j)-t / SINH (t*dx(j)))/tsq
   20   r(j) = (y(j+1)-y(j))/dx(j)-(y(j)-y(j-1))/dx(j-1)
c
c --- boundary conditions on spline
c
        th2 = t*dx(1)
        thn = t*dx(nm1)
        b1 = (1.0/t)*(1.0 / SINH (th2)-1.0/th2)
        a1 = -1.0/tanh(th2)+1.0/th2
        a1 = a1/(1.0 / SINH (th2)-1.0/th2)
        an = (1.0/t)*(-1.0 / SINH (thn)+1.0/thn)
        bn = 1.0 / tanh(thn)-1.0/thn
        bn = bn/(-1.0 / SINH (thn)+1.0/thn)
        if (bpar(1) .ne. -1.0e30)  go to 30
        bp1 = 1.0
****    bp2 = (-1.0/(dx(1)*b1))*(dx(1)*bpar(2)+y(1)-y(2))
        bp2 = (+1.0/(dx(1)*b1))*(dx(1)*bpar(2)+y(1)-y(2))
   30   if (bpar(3) .ne. -1.0e30)  go to 40
        bp3 = 1.0
        bp4 = (1.0/(dx(nm1)*an))*(dx(nm1)*bpar(4)+y(nm1)-y(nx))
   40   b(1) = a1
        c(1) = bp1
        a(nx) = bp3
        b(nx) = bn
        r(1) = bp2
        r(nx) = bp4
        go to 120
c
c --- for tension sufficiently small go over into cubic spline
c
   50 do 110 j=2,nm1
      dx(j) = x(j+1)-x(j)
      a(j) = dx(j-1)/6.0
      b(j) = (dx(j)+dx(j-1))/3.0
      c(j) = dx(j)/6.0
      if (t .eq. 0.0)  go to 110
      a(j) = a(j)-1.94444444e-02*tsq*dx(j-1)**3
      b(j) = b(j)-tsq*(dx(j-1)**3+dx(j)**3)/45.0
      c(j) = c(j)-1.94444444e-02*tsq*dx(j)**3
  110 r(j) = (y(j+1)-y(j))/dx(j)-(y(j)-y(j-1))/dx(j-1)
c
c --- boundary conditions
c
      if (bpar(1) .ne. -1.0e30)  go to 60
      bp1 = 1.0
      bc = dx(1)*(-0.166666667+1.94444444e-02*dx(1)**2*tsq)
      bp2 = (1.0/(dx(1)*bc))*(y(1)-y(2)+dx(1)*bpar(2))
   60 if (bpar(3) .ne. -1.0e30)  go to 70
      bp3 = 1.0
      ac = dx(nm1)*(0.166666667-1.94444444*dx(nm1)**2*tsq)
      bp4 = (1.0/(dx(nm1)*ac))*(dx(nm1)*bpar(4)+y(nm1)-y(nx))
   70 a(nx) = bp3
      b( 1) = 2.0 + dx(1  )**2*tsq / 10.0
      c( 1) = bp1
      b(nx) = 2.0 + dx(nm1)**2*tsq / 10.0
      r( 1) = bp2
      r(nx) = bp4
c
c --- note a(1) and c(nx) are not used
c --- r is vector of rhs
c
c --- solve the tridiagonal system
c --- forward elimination
c
  120 do j=2,nx
        pv = a(j)/b(j-1)
        b(j) = b(j)-c(j-1)*pv
        r(j) = r(j)-r(j-1)*pv
      end do
c
c --- back substitution
c
      fpp(nx) = r(nx)/b(nx)
      do j=1,nm1
        i = nx-j
        fpp(i) = (r(i)-c(i)*fpp(i+1))/b(i)
      end do
c
c --- now have vector of second derivatives fpp. set up the vector cs(i,j)
c --- to be used in evaluating the tension spline.
c --- convention used for natural spline will not work here.
c --- however we still pack all the necessary information into array cs.
c
      do 150 j=1,nx
        cs(j,1) = fpp(j)
        if (j .eq. nx )  go to 150
        if (t .eq. 0.0)  go to 150
        cs(j,2) = 1.0 / (SINH (t*dx(j))*t*t)
        cs(j,3) = 1.0 / (dx(j)*t*t)
  150 continue
c
      call evtspl (x, y, nx, cs, ic, rgrid, tspl, npts, t)
      if (.not. reverseset)  return
c
 1000 do j=1,nx/2
        xhold     = x(j)
        x(j)      = x(nx-j+1)
        x(nx-j+1) = xhold
        yhold     = y(j)
        y(j)      = y(nx-j+1)
        y(nx-j+1) = yhold
      end do
c
      reverseset  = .not. reverseset
      if (reverseset)  go to 16
      return
c
      end
!***********************************************************************
      subroutine evtspl (x, y, nx, cs, ic, r, tspl, npts, t)
c
      implicit  integer (i-n), real (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- evaluate the tension spline at npts points of r, return results in tspl
c ----------------------------------------------------------------------
c
      dimension  x(*), y(*), cs(ic,3), r(*), tspl(*)
c
      do j=1,npts
        z = r(j)
        do i=1,nx-1
          if (z .le. x(i+1))  go to 30
        end do
        if (ABS ((z-x(nx))/x(nx)) .lt. 1.0e-8)  go to 30
c
c       call STOP ('subroutine EVTSPL: unspecified problem', 124)
c
c       evaluate the spline
c
   30   dxl = z-x(i)
        dxr = x(i+1)-z
        if (t .ne. 0.0) then
             tspl(j) = cs(i,1)*cs(i,2) * SINH (t*dxr)+
     .                 (y(i)*t*t-cs(i,1))*cs(i,3)*dxr
     .                +cs(i+1,1)*cs(i,2) * SINH (t*dxl)+
     .                 (y(i+1)*t*t-cs(i+1,1))*cs(i,3)*dxl
        else
             delx    = x(i+1)-x(i)
             tspl(j) = (dxr*(-cs(i  ,1) * dxl*(dxr+delx)+6.0 * y(i  ))
     .                + dxl*(-cs(i+1,1) * dxr*(dxl+delx)+6.0 * y(i+1)))
     .                 /(6.0 * delx)
        end if
      end do
      return
c
      end
