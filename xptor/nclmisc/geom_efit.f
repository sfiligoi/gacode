      SUBROUTINE GEOM_EFIT(device,idshot,time,nr_r,rhop_r,f_r,q_r,
     &                     rhot_r,rin_r,rout_r,elong_r,vol_r,vp_r,
     &                     phit_r,fm_r,grth_r,gph_r,gth_r,b2_r,bm2_r,
     &                     grho1_r,grho2_r,gr2bm2_r,ftrap_r,fhat_r,
     &                     bpout_r,btout_r,psi_r,a0,r0,bt0,rm1_r,
     &                     rm2_r,rhor_r,iflag,message)
!***********************************************************************
!GEOM_EFIT generates plasma geometry information for FORCEBAL by calling
!  getting EQDSK information about the MHD equilibrium and then calling
!  FLUXAV_EFIT to generate additional flux surface 
!  quantities 
!  W.A.Houlberg 3/2000
!Input:
!  device-experimental device
!  idshot-shot number/id
!  time-analysis time
!Output:
!  nr_r-number of radial points
!  rhop_r-normalized poloidal flux grid proportional to psi
!  f_r(i)-2*pi*R*B_t/mu0 (A)
!  q_r(i)-safety factor
!  rhot_r(i)-normalized toroidal flux grid proportional to (Phi)**0.5
!  rin_r(i)-major radius grid on inside of torus in axis plane (m)
!  rout_r(i)-major radius grid on outside of torus in axis plane (m)
!  elong_r(i)-elongation
!  vol_r(i)-volume enclosed (m**3)
!  vp_r(i)-d vol_r/d rhot_r/a0 (m**2)
!  phit_r(i)-toroidal flux (Wb)
!  fm_r(3,i)-geometric factor
!  grth_r(i)-<n.grad(theta)> (1/m)
!  gph_r(i)-poloidal flux metric (-)
!  gth_r(i)-toroidal flux metric (-)
!  b2_r(i)-<B**2> (T**2)
!  bm2_r(i)-<1/B**2> (/T**2)
!  grho1_r(i)-a0*<|grad(rhot_r)|>
!  grho2_r(i)-a0**2*<|grad(rhot_r)|**2>
!  gr2bm2_r(i)-a0**2*<|grad(rhot_r)|**2/B**2> (1/T**2)
!  ftrap_r(i)-trapped fraction
!  bpout_r(i)-poloidal field at rout_r(i) (T)
!  btout_r(i)-toroidal field at rout_r(i) (T)
!  psi_r(i)-poloidal flux (Wb/rad)
!  a0-minor radius, half diameter of boundary flux surface (m)
!  r0-major radius, center of boundary flux suface (m)
!  bt0-toroidal field at r0 (T)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  device
      INTEGER        idshot
      REAL           time
!Declaration of output variables
!  0-D
      CHARACTER*(*)  message
      INTEGER        nr_r,                    iflag
      REAL           a0,                      r0,
     &               bt0
!  1-D
      REAL           b2_r(*),                 bm2_r(*),
     &               bpout_r(*),              btout_r(*),
     &               elong_r(*),              f_r(*),
     &               fhat_r(*),               fm_r(3,*),
     &               ftrap_r(*),              gph_r(*),
     &               grho1_r(*),              grho2_r(*),
     &               gr2bm2_r(*),             grth_r(*),
     &               gth_r(*),                phit_r(*),
     &               psi_r(*),                q_r(*),
     &               rhop_r(*),               rhot_r(*),
     &               rin_r(*),                rout_r(*),
     &               vol_r(*),                vp_r(*),
     &               rhor_r(*)
!Declaration of internal variables
      CHARACTER      cidshot*9,               ctimems*5,
     &               cnin*25
      INTEGER        i,                       itimems
      REAL           bmag,                    const,
     &               rmag,                    zmag
      REAL           z_mu0,                   z_pi
!  Two dimensional poloidal flux-xy
      INTEGER        mxnx_xy,                 mxny_xy,
     &               nx_xy,                   ny_xy
      PARAMETER     (mxnx_xy=130,             mxny_xy=130)
!  Radial grid-psi
      INTEGER        mxnr_r
      PARAMETER     (mxnr_r=300)
!  Limiter surface-lim
      INTEGER        mxn_lim,                 n_lim
      PARAMETER     (mxn_lim=200)
!  Plasma boundary-bdry
      INTEGER        mxn_bdry,                n_bdry
      PARAMETER     (mxn_bdry=1500)
!  1-D psi
      REAL           ffp_r(mxnr_r),           frb_r(mxnr_r),
     &               p_r(mxnr_r),             pp_r(mxnr_r),
     &               r2_r(mxnr_r),            rm2_r(mxnr_r),
     &               rm1_r(mxnr_r)
!  2-D x,y                                                    
      REAL           x_xy(mxnx_xy),           y_xy(mxny_xy),
     &               psi_xy(mxnx_xy,mxny_xy)
!  Limiter
      REAL           x_lim(mxn_lim),          y_lim(mxn_lim)
!  Plasma boundary
      REAL           x_bdry(mxn_bdry),        y_bdry(mxn_bdry)
!Initialization
!  Physical and conversion constants
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
!Get EFIT data
      IF(device.eq.'cmod') THEN
!       Read MDS files for CMOD EFIT MHD equilibrium information
c--        CALL READ_EFIT_MDS(idshot,time,mxnx_xy,mxny_xy,mxnr_r,mxn_bdry,
c--     &                     mxn_lim,nx_xy,ny_xy,nr_r,r0,bt0,rmag,zmag,
c--     &                     bmag,x_xy,y_xy,psi_xy,psi_r,rhop_r,frb_r,
c--     &                     ffp_r,p_r,pp_r,q_r,n_bdry,x_bdry,y_bdry,
c--     &                     n_lim,x_lim,y_lim,iflag,message)
      ELSEIF(device.eq.'jet') THEN
!       Read PPF files for JET EFIT MHD equilibrium information
c--        CALL READ_EFIT_JET(idshot,time,mxnx_xy,mxny_xy,mxnr_r,mxn_bdry,
c--     &                     mxn_lim,nx_xy,ny_xy,nr_r,r0,bt0,rmag,zmag,
c--     &                     bmag,x_xy,y_xy,psi_xy,psi_r,rhop_r,frb_r,ffp_r,
c--     &                     p_r,pp_r,q_r,n_bdry,x_bdry,y_bdry,n_lim,x_lim,
c--     &                     y_lim,iflag,message)
      ELSE
!       Read EQDSK file to get EFIT MHD equilibrium information
        itimems=1.0e3*(time+.00001)
        WRITE(ctimems,'(i5)') itimems
        DO i=1,5
          IF(ctimems(i:i).eq.' ') ctimems(i:i)='0'
        ENDDO
        WRITE(cidshot,'(i9)') idshot
        DO i=1,9
          IF(cidshot(i:i).eq.' ') cidshot(i:i)='0'
        ENDDO
        cnin='g'//cidshot(4:9)//'.'//ctimems
        CALL READ_EFIT_EQDSK(cnin,mxnx_xy,mxny_xy,mxnr_r,mxn_bdry,
     &                       mxn_lim,nx_xy,ny_xy,nr_r,r0,bt0,rmag,zmag,
     &                       bmag,x_xy,y_xy,psi_xy,psi_r,rhop_r,frb_r,
     &                       ffp_r,p_r,pp_r,q_r,n_bdry,x_bdry,y_bdry,
     &                       n_lim,x_lim,y_lim,iflag,message)
      ENDIF
!Check warning and error flag
      IF(iflag.gt.1) THEN
        i=LEN(message)
        message='GEOM_EFIT(1)/'//message(1:i-13)
        GOTO 1000
      ENDIF
!Call FLUXAV_EFIT to generate metrics from EFIT MHD equilibrium
      CALL FLUXAV_EFIT(bmag,rmag,zmag,nx_xy,ny_xy,x_xy,y_xy,psi_xy,
     &                 n_lim,x_lim,y_lim,nr_r,psi_r,frb_r,ffp_r,q_r,
     &                 rin_r,rout_r,rhot_r,elong_r,gth_r,gph_r,vol_r,
     &                 vp_r,r2_r,rm2_r,phit_r,fm_r,grth_r,b2_r,bm2_r,
     &                 grho1_r,grho2_r,gr2bm2_r,ftrap_r,fhat_r,bpout_r,
     &                 btout_r,a0,rm1_r,rhor_r,iflag,message)
!Check warning and error flag
      do i=1,nr_r
c       write(*,222) i, rhot_r(i), rhop_r(i)
c       write(*,222) i, rhop_r(i), psi_r(i), p_r(i), q_r(i)
c       write(*,222) i, rhot_r(i), rm1_r(i), rm2_r(i), frb_r(i)
      enddo
 222  format(2x,i2,1p4e12.4,' geom_efit')
      IF(iflag.gt.1) THEN
        i=LEN(message)
        message='GEOM_EFIT(2)/'//message(1:i-13)
        GOTO 1000
      ENDIF
!At this point q_r is always positive
!Correct sign on q_r for coordinate consistency
      const=btout_r(nr_r)/bpout_r(nr_r)*ABS(bpout_r(nr_r)/btout_r(nr_r))
      CALL RARRAY_SCALE(nr_r,const,q_r,1)
!Get poloidal current by scaling frb=r*Bt
      const=2.0*z_pi/z_mu0
      CALL RARRAY_COPY(nr_r,frb_r,1,f_r,1)
      CALL RARRAY_SCALE(nr_r,const,f_r,1)
 1000 RETURN
      END
