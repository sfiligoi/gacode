c---------------------------
      subroutine eqi_frho1(ifnum,ict,ivec,ivecd,zans,indx,zparam,
     >   zhrho,zhrhoi)
c
c  interpolation on one function, any order
c
      use eq_xyz_cmn
      implicit NONE
c
      integer ifnum                     ! function id number
      integer ivec                      ! input vector dimension
      integer ivecd                     ! output vector dimension
      integer ict(*)                    ! value/derivative selector array
      real*8 zans(ivecd,*)              ! answer(s) returned
c
      integer indx(ivec)                ! rho-zone index
      real*8 zparam(ivec)               ! normalized displacement w/in zone
      real*8 zhrho(ivec)                ! zone width
      real*8 zhrhoi(ivec)               ! inverse zone width
c
c  the last 4 arguments set by a prior "eqi_lkup" call...
c
c--------------------------------
c
      integer iadr,istype,ix,isize,i
      real*8 ztmp(ivec)
c
c--------------------------------
c
      iadr=xp_items(ifnum)%mag_addr
      istype=xp_items(ifnum)%mag_order
      ix=xp_items(ifnum)%mag_rho
      isize=axes(ix)%size
c
      if(istype.eq.3) then
c  explicit spline
         ztmp=zhrho*zparam
         call r8cspevfn(ict,ivec,ivecd,zans,indx,ztmp,eqbuf(iadr),isize)
      else if(istype.eq.2) then
c  spline
         call r8fvspline(ict,ivec,ivecd,zans,indx,zparam,zhrho,zhrhoi,
     >      eqbuf(iadr),isize)
      else if(istype.eq.1) then
c  hermite
         call r8herm1fcn(ict,ivec,ivecd,zans,indx,zparam,zhrho,zhrhoi,
     >      eqbuf(iadr),isize)
c  piecewise linear
      else if(istype.eq.0) then
         call r8pc1fcn(ict,ivec,ivecd,zans,indx,zparam,zhrho,zhrhoi,
     >      eqbuf(iadr),isize)
c  step fcn
      else
         do i=1,ivec
            zans(i,1)=czero
            if(ict(1).eq.1) zans(i,1)=eqbuf(iadr+indx(i)-1)
         enddo
         if(ict(2).eq.1) write(lunerr,*)
     >      ' %f(rho): no 1st derivative, ',xp_names(ifnum),
     >      ' is a step function.'
         if(ict(3).eq.1) write(lunerr,*)
     >      ' %f(rho): no 2nd derivative, ',xp_names(ifnum),
     >      ' is a step function.'
      endif
c
      return
      end
      subroutine eqi_frho(ivec,zrho,ifnum,ict,ivecd,zans,ier)
c
      use eq_xyz_cmn
c
c  **vectorized**
c  evaluate 1d spline (and derivatives) vs. rho
c
c  input:
c
      IMPLICIT NONE
c
      integer ivec                      ! vector dimension
      REAL*8 zrho(ivec)                 ! argument -- where to evaluate
      integer ifnum                     ! function to evaluate
      integer ict(*)                    ! fcn/derivative info
c
c  output:
c
      integer ivecd                     ! output vector dimension
      REAL*8 zans(ivecd,*)              ! result of evaluation(s)
      integer ier                       ! exit code, 0= normal
c
c----------------------------------
c
      REAL*8 zparam(ivec),zhrho(ivec),zhrhoi(ivec)
      integer indx(ivec),iadr,istype,iftype,ix1
c
c--------------------------------------------------------------------
c
      iftype=xp_items(ifnum)%mag_dim
      if(iftype.ne.1) then
         write(lunerr,*)
     >      ' ?f(rho) interpolation:  not a function of just rho:  ',
     >      xp_names(ifnum)
         ier=1
         return
      endif
c
      ier=0
      ix1=xp_items(ifnum)%mag_rho
      if(ix1.ne.iaxis_rho) then
         if(ix1.eq.0) then
            ier=1
         else
            if(axes(ix1)%kin.ne.iaxis_rho) then
               ier=1
            endif
         endif
      endif
      if(ier.eq.1) then
         write(lunerr,*)
     >      ' ?f(rho) interpolation:  invalid rho axis w/fcn:  ',
     >      xp_names(ifnum)
         return
      endif
c
      call eqi_binget(ix1,ivec,zrho,
     >   indx,zparam,zhrho,zhrhoi,ier)
      if(ier.ne.0) return
c
c  interpolation...
c
      call eqi_frho1(ifnum,ict,ivec,ivecd,zans,indx,zparam,zhrho,zhrhoi)
c
      return
      end
