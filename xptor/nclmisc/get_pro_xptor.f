      SUBROUTINE GET_PRO_XPTOR(nout,idshot,igrid_ncl,time,
     #                         bt0,izi0,izim1,izim2,nr_r,
     #                         mxnr_r,rho_r,rhot_r,te_r,ti_r,den_r,
     #                         deni_r,zeff_ex_r,denim_r,denim2_r,
     #                         denfst_r,
     #                         omega_ex_r,vp_im_ex_r,vt_im_ex_r,
     #                         xk_im_ex_r,zeff_r,e_par_ex_r,xj_ex_r,
     #                         xj_nb_ex_r,iflag)
!***********************************************************************
!GET_PRO_XPTOR gets profile data from XPTOR variables
!  J. Kinsey 6/29/00
!***********************************************************************
      IMPLICIT NONE
      include '../inc/input.m'
      include '../inc/tport.m'
!Declaration of input variables
      INTEGER        idshot,                  izi0(*),
     #               izim1,                   izim2,
     #               mxnr_r,                  nout,
     #               nr_r,                    igrid_ncl
      REAL           bt0,                     rhot_r(*),
     #               time
!Declaration of output variables
      INTEGER        iflag
      REAL           den_r(*),                deni_r(mxnr_r,*),
     #               denim_r(*),              denim2_r(mxnr_r,*),
     #               e_par_ex_r(*),           omega_ex_r(*),
     #               te_r(*),                 ti_r(*),
     #               vp_im_ex_r(*),           vt_im_ex_r(*),
     #               zeff_ex_r(*),            zeff_r(*),
     #               xj_ex_r(*),              xj_nb_ex_r(*),
     #               xk_im_ex_r(*),           denfst_r(*),
     #               rho_r(*)
!Declaration of local variables
!  Profiles
      INTEGER        mxnr_4d,                 nr_4d
      PARAMETER     (mxnr_4d=101)
      REAL           r_4d(mxnr_4d),           y_4d(mxnr_4d)
      INTEGER        mxn
      PARAMETER     (mxn=130)
      REAL           denb_ex_r(mxn)
      real           r_xpt(jmaxmt+1), te_xpt(jmaxmt+1),
     #               ti_xpt(jmaxmt+1), den_xpt(jmaxmt+1),
     #               deni_xpt(jmaxmt+1), zeff_xpt(jmaxmt+1),
     #               nimp_xpt(jmaxmt+1), vphi_xpt(jmaxmt+1),
     #               nfst_xpt(jmaxmt+1)
!  Other
      CHARACTER*120 msg
      CHARACTER      cidshot*9,               cnin*25,
     #               ctimems*5
      INTEGER        i,                       ierr,
     #               itimems,                 jpro,
     #               k
      REAL           fscale,                  yim,
     #               yim2
      iflag=0
!Initialize arrays
      CALL RARRAY_ZERO(nr_r,rho_r)
      CALL RARRAY_ZERO(nr_r,te_r)
      CALL RARRAY_ZERO(nr_r,ti_r)
      CALL RARRAY_ZERO(nr_r,den_r)
      CALL RARRAY_ZERO(nr_r,denfst_r)
      CALL RARRAY_ZERO(nr_r,zeff_ex_r)
      CALL RARRAY_ZERO(nr_r,omega_ex_r)
      CALL RARRAY_ZERO(nr_r,denim_r)
      CALL RARRAY_ZERO(mxnr_r*izim2,denim2_r)
      CALL RARRAY_ZERO(nr_r,vt_im_ex_r)
      CALL RARRAY_ZERO(nr_r,vp_im_ex_r)
      CALL RARRAY_ZERO(nr_r,xk_im_ex_r)
      CALL RARRAY_ZERO(nr_r,e_par_ex_r)
      CALL RARRAY_ZERO(nr_r,xj_ex_r)
      CALL RARRAY_ZERO(nr_r,xj_nb_ex_r)
!
      if(igrid_ncl.eq.0) then  ! interpolate to MHD grid
        do i=0,jmaxm
          r_xpt(i+1)=rho(i)
          te_xpt(i+1)=te_m(i)
          ti_xpt(i+1)=ti_m(i)
          den_xpt(i+1)=ne_m(i)*1.e19
          deni_xpt(i+1)=ni_m(i)*1.e19  ! total ion density
          zeff_xpt(i+1)=zeff_exp(i)
          nimp_xpt(i+1)=nz_exp(i)*1.e19
          nfst_xpt(i+1)=nfast_exp(i)*1.e19
          vphi_xpt(i+1)=rmajor_exp*angrot_exp(i)
        enddo
        call w_lin_interp(jmaxm+1,r_xpt,te_xpt,nr_r,rhot_r,te_r,
     &     iflag,msg)
        call w_lin_interp(jmaxm+1,r_xpt,ti_xpt,nr_r,rhot_r,ti_r,
     &     iflag,msg)
        call w_lin_interp(jmaxm+1,r_xpt,den_xpt,nr_r,rhot_r,den_r,
     &     iflag,msg)
c       call wlterp(jmaxmt+1,r_xpt,deni_xpt,nr_r,rhot_r,deni_r,iflag)
        call w_lin_interp(jmaxm+1,r_xpt,zeff_xpt,nr_r,rhot_r,zeff_ex_r,
     &     iflag,msg)
        call w_lin_interp(jmaxm+1,r_xpt,nimp_xpt,nr_r,rhot_r,denim_r,
     &     iflag,msg)
        call w_lin_interp(jmaxm+1,r_xpt,nfst_xpt,nr_r,rhot_r,denfst_r,
     &     iflag,msg)
        call w_lin_interp(jmaxm+1,r_xpt,vphi_xpt,nr_r,rhot_r,vt_im_ex_r,
     &     iflag,msg)
c
c... use XPTOR profile grid
c    Note: rhot_r still on MHD grid, not reset until after
c    interpolation in forcebal.f
c
      else
        do i=1,nr_r
          rho_r(i)=rho(i-1)
          te_r(i)=te_m(i-1)
          ti_r(i)=ti_m(i-1)
          den_r(i)=ne_m(i-1)*1.e19
          deni_r(i,1)=ni_m(i-1)*1.e19  ! total ion density
          zeff_ex_r(i)=zeff_exp(i-1)
          denim_r(i)=nz_exp(i-1)*1.e19
          denfst_r(i)=nfast_exp(i-1)*1.e19
          vt_im_ex_r(i)=rmajor_exp*angrot_exp(i-1)
        enddo
      endif
c
!Add an ad hoc test species at small percentage of electron density
      do k=2,5
        if(izi0(k).ne.0) then
          do i=1,nr_r
            deni_r(i,k)=1.0e-5*den_r(i)
          enddo
        endif
      enddo
!Break electron density into ion and impurity components
!Here, ion density is thermal only with fast ions removed in deni_r(i,1)
      do i=1,nr_r
!       Diagnostic impurity
        yim=izim1*denim_r(i)
        yim2=izim1**2*denim_r(i)
!       Other than main ion
        do k=2,5
          yim=yim+izi0(k)*deni_r(i,k)
          yim2=yim2+izi0(k)**2*deni_r(i,k)
        enddo
!       Multiple charge state impurity
        do k=1,izim2
          yim=yim+k*denim2_r(i,k)
          yim2=yim2+k**2*denim2_r(i,k)
        enddo
c        deni_r(i,1)=(den_r(i)-yim-denb_ex_r(i)-denfst_r(i))/izi0(1)
        deni_r(i,1)=(den_r(i)-yim-denfst_r(i))/izi0(1)
        zeff_r(i)=(izi0(1)*deni_r(i,1)+yim2+denfst_r(i))/den_r(i)
      enddo
c
c     write(*,*) 'jmaxm = ',jmaxm
c     write(*,*) 'nr_r = ',nr_r
c     do i=0,jmaxm
c        write(*,100) i, rho(i),ne_m(i),ni_m(i),nz_exp(i)
c     enddo
c     do i=0,jmaxm
c       write(*,105) i, rho(i), rmajor_exp*angrot_exp(i)
c     enddo
c     do i=1,nr_r
c       write(*,105) i, rhot_r(i), den_r(i), deni_r(i,1), denim_r(i)
c   .               denfst_r(i), denim_r(i), zeff_ex_r(i), zeff_r(i)
c     enddo
c     do i=1,nr_r
c       write(*,105) i, rhot_r(i), vt_im_ex_r(i)
c       write(*,105) i, rhot_r(i), te_r(i), ti_r(i)
c     enddo
c
  50  format(i2,2x,0p1f4.2,0p6f10.5,2x,'get_pro_xptor')
 100  format(i2,2x,0p6f10.5)
 105  format(i2,2x,0p1f10.5,1p6e13.4)
 1000 RETURN
      END
