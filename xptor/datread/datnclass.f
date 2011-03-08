c@datclass.f
c jek 18-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... calls NCLASS neoclassical using exp data in _d variables
c    use_xneo_m=2(exp profiles) or =3(model profiles)
c    iexb=3 to use NCLASS for ExB shear (Hahm-Burrell)
c    v_nt = nondiffusive pinch velocity, v_eb=Ware pinch velocity
c    read in NCLASS geometry variables in readxp.f if ncl_flag=1
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine datnclass
c
      implicit none
      include 'mpif.h'
      include '../inc/tport.m'
      include '../inc/data.m'
      include '../inc/model.m'
      include '../inc/input.m'
      include '../inc/glf.m'
c
      integer j, xneo, jgrid, iflag_ncl, mx_ni, mxnr_r, iflag
      parameter(mx_ni=5,mxnr_r=300)
      include '../inc/pamx_ms.inc'
      include '../inc/pamx_mz.inc'
      include '../inc/comfbl.inc'
      character*3 ci0(mx_ni), cim1, cim2
      character*15 device
      character*120 msg
      real time_ncl
      real rhot_r(mxnr_r), rhop_r(mxnr_r), q_r(mxnr_r), frb_r(mxnr_r), 
     &     rm1_r(mxnr_r), rm2_r(mxnr_r), vol_r(mxnr_r)
      real*8 chie_ncl_r(mxnr_r), chii_ncl_r(mxnr_r),
     &       vphi_ncl_r(mxnr_r), vpol_ncl_r(mxnr_r),
     &       di_ncl_r(mxnr_r), de_ncl_r(mxnr_r),
     &       ero_ncl_r(mxnr_r), omexb_ncl_r(mxnr_r),
     &       vne_ncl_r(mxnr_r), vni_ncl_r(mxnr_r),
     &       frb_ncl_r(mxnr_r)
      real*8 game_ncl_r(mxnr_r), gami_ncl_r(mxnr_r),
     &       qe_ncl_r(mxnr_r), qi_ncl_r(mxnr_r)
      real*8 chie_ncl_d(nj), chii_ncl_d(nj),
     &       vphi_ncl_d(nj), vpol_ncl_d(nj),
     &       di_ncl_d(nj), de_ncl_d(nj),
     &       ero_ncl_d(nj), omexb_ncl_d(nj),
     &       vne_ncl_d(nj), vni_ncl_d(nj)
      real*8 game_ncl_d(nj), gami_ncl_d(nj),
     &       qe_ncl_d(nj), qi_ncl_d(nj)
      real*8 r_mh(mxnr_r,2), vol_mh(mxnr_r,2)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      if(i_proc.eq.0) then
        if(use_xneo_m.ge.2.and.iexb.eq.3) then
          write(*,40) use_xneo_m
        elseif(use_xneo_m.ge.2.and.iexb.ne.3) then
          write(*,42) use_xneo_m
        else
          write(*,44) use_xneo_m
        endif
      endif
c
c get densities and temperatures needed by get_pro_xptor
c Note: no ascale and bscale factors included here
c
      xneo=0
      jgrid=1     ! interpolate MHD info onto transport grid
      device='xptor'
      ishot_d=shotno
      time_ncl=xp_time
      rmajor_exp = rmajor_d
      do j=1,mx_ni
        ci0(j)='   '
      enddo
      ci0(1)='D'
      ci0(2)='N'
      cim1='C'
      cim2='   '
c
      do j=1,nj_d
        zeff_exp(j-1)=zeff_d(j)
        nz_exp(j-1)=1.D-19*en_d(j,nprim_d+1)
        nfst_exp(j-1)=1.D-19*enbeam_d(j)
        angrot_exp(j-1)=angrot_d(j)
      enddo
c
      do j=0,jmaxm
        te_m(j)=te_d(j+1)
        ti_m(j)=ti_d(j+1)
        ne_m(j)=1.D-19*ene_d(j+1)
        ni_m(j)=1.D-19*en_d(j+1,1)
      enddo
c
      call forcebal(jgrid,jmaxm,xneo,device,ishot_d,time_ncl,
     >              ci0(1),cim1,cim2,nr_xp,rhot_r,rhop_r,
     >              q_r,frb_r,rm1_r,rm2_r,vol_r,iflag_ncl)
c
c...mhd variables
c   1=zone center, 2=zone boundary
c
      do j=1,nr_xp-1
         r_mh(j,1)=(rhot_r(j)+rhot_r(j+1))/2.D0
         r_mh(j,2)=rhot_r(j)
         vol_mh(j,1)=(vol_r(j)+vol_r(j+1))/2.D0
         vol_mh(j,2)=vol_r(j)
      enddo
       r_mh(nr_xp,2)=rhot_r(nr_xp)
       vol_mh(nr_xp,2)=vol_r(nr_xp)
c
      do j=1,nr_xp
         rho_ncl_r(j)=rhot_r(j)
         game_ncl_r(j)=gam_r(j,1)/1.D19
         gami_ncl_r(j)=gam_r(j,2)/1.D19
         qe_ncl_r(j)=q_con_r(j,1)
         qi_ncl_r(j)=q_con_r(j,2)
c        de_ncl_r(j)=d_eff_g_r(j,1)
c        di_ncl_r(j)=d_eff_g_r(j,2)
         de_ncl_r(j)=d_n_g_r(j,1) ! elec diagonal ptcle diffusivity (m**2/s)
         di_ncl_r(j)=d_n_g_r(j,2) ! D ion diagonal ptcle diffusivity (m**2/s)
         chie_ncl_r(j)=chi_eff_g_r(j,1) ! effective elec heat diffusivity (m**2/s)
         chii_ncl_r(j)=chi_eff_g_r(j,2) ! effective D ion heat diffusivity (m**2/s)
         vne_ncl_r(j)=v_nt_g_r(j,1)+v_eb_g_r(j,1) ! elec ptcle pinch velocity (m/s)
         vni_ncl_r(j)=v_nt_g_r(j,2)+v_eb_g_r(j,2) ! D ptcle pinch velocity (m/s)
         vphi_ncl_r(j)=v_t_o_r(j,2) ! v_t D rotation (m/s)
         vpol_ncl_r(j)=v_p_o_r(j,2) ! v_p D rotation (m/s)
         ero_ncl_r(j)=e_rad_o_r(j,1)
         omexb_ncl_r(j)=omexb_o_r(j) ! (rad/s)
         frb_ncl_r(j)=abs(frb_r(j))
c          write(*,60) j, rhot_r(j), q_con_r(j,1), q_con_r(j,2),
c    >                 q_con_r(j,3), q_con_r(j,4)
c           write(*,60) j, rho_ncl_r(j), gami_ncl_r(j)
c          write(*,60) j, rho_ncl_r(j), omexb_ncl_r(j)
c          write(*,60) j, rho_ncl_r(j), chie_ncl_r(j), chii_ncl_r(j)
c          write(*,60) j, rho_ncl_r(j), de_ncl_r(j), di_ncl_r(j)
c          write(*,60) j, rho_ncl_r(j), vne_ncl_r(j), vni_ncl_r(j)
c          write(*,60) j, rho_ncl_r(j), q_con_r(j,1), q_con_r(j,2)
c          write(*,60) j, rho_ncl_r(j), v_t_o_r(j,2), v_p_o_r(j,2)
c          write(*,60) j, rho_ncl_r(j), v_p_o_r(j,2), v_t_o_r(j,2),
c    >                 v_t_o_r(j,3), v_t_o_r(j,4),
c    >                 rmajor_d*angrot_d(j)
      enddo
c
c...interpolate to transport grid
      if(jgrid.eq.0) then
c
        call w_lin_interp_r8(nr_xp,rho_ncl_r,game_ncl_r,
     >                    nj_d,rho_d,game_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,gami_ncl_r,
     >                    nj_d,rho_d,gami_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,qe_ncl_r,
     >                    nj_d,rho_d,qe_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,qi_ncl_r,
     >                    nj_d,rho_d,qi_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,chie_ncl_r,
     >                    nj_d,rho_d,chie_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,chii_ncl_r,
     >                    nj_d,rho,chii_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,de_ncl_r,
     >                    nj_d,rho_d,de_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,di_ncl_r,
     >                    nj_d,rho_d,di_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,vne_ncl_r,
     >                    nj_d,rho,vne_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,vni_ncl_r,
     >                    nj_d,rho,vni_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,vphi_ncl_r,
     >                    nj_d,rho,vphi_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,vpol_ncl_r,
     >                    nj_d,rho,vpol_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,ero_ncl_r,
     >                    nj_d,rho,ero_ncl_d,iflag,msg)
        call w_lin_interp_r8(nr_xp,rho_ncl_r,omexb_ncl_r,
     >                    nj_d,rho,omexb_ncl_d,iflag,msg)
        do j=1,nj_d
           chieneo_exp(j-1)=chie_ncl_d(j)
           chiineo_exp(j-1)=chii_ncl_d(j)
           deneo_exp(j-1)=de_ncl_d(j)
           dineo_exp(j-1)=di_ncl_d(j)
           veneo_exp(j-1)=vne_ncl_d(j)
           vineo_exp(j-1)=vni_ncl_d(j)
           vphi_ncl_exp(j-1)=vphi_ncl_d(j)
           vpol_ncl_exp(j-1)=vpol_ncl_d(j)
           omexb_ncl_exp(j-1)=omexb_ncl_d(j)
c          write(*,60) j, rho_d(j), de_ncl_d(j), di_ncl_d(j)
c          write(*,60) j, rho_d(j), chie_ncl_d(j), chii_ncl_d(j)
        enddo
      else
        do j=1,nj_d
           chieneo_exp(j-1)=chie_ncl_r(j)
           chiineo_exp(j-1)=chii_ncl_r(j)
           deneo_exp(j-1)=de_ncl_r(j)
           dineo_exp(j-1)=di_ncl_r(j)
           veneo_exp(j-1)=vne_ncl_r(j)
           vineo_exp(j-1)=vni_ncl_r(j)
           vphi_ncl_exp(j-1)=vphi_ncl_r(j)
           vpol_ncl_exp(j-1)=vpol_ncl_r(j)
           omexb_ncl_exp(j-1)=omexb_ncl_r(j)
c          write(*,60) j, rho_d(j), de_ncl_d(j), di_ncl_r(j)
c          write(*,60) j, rho_d(j), chie_ncl_r(j), chii_ncl_r(j)
        enddo
      endif
c
c     do j=0,jmaxm
c       write(*,60) j, rho(j), chieneo_exp(j), chiineo_exp(j)
c       write(*,60) j, rho(j), dineo_exp(j), vineo_exp(j)
c     enddo
      if (itest_ncl.gt.0 .and. i_proc.eq.0) then
        if (i_proc.eq.0) write(*,55)
        do j=1,jmaxm
           write(*,60) j, rho(j), chieneo_exp(j), chiineo_exp(j)
        enddo
        if (i_proc.eq.0) write(*,56)
        do j=1,jmaxm
           write(*,60) j, rho(j), deneo_exp(j), dineo_exp(j)
        enddo
c
        write(*,65)
        do j=1,nr_xp
           write(*,60) j, rho_ncl_r(j), v_t_o_r(j,2), v_p_o_r(j,2),
     >     e_rad_o_r(j,1)
        enddo
        write(*,70)
        do j=1,nj_d
           write(*,60) j, rho_d(j), vphi_ncl_d(j), vpol_ncl_d(j),
     >     ero_ncl_d(j)
        enddo
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 40   format(' Using NCLASS for neoclassical transport, ExB shear, ',
     >       'use_xneo_m =',i2)
 42   format(' Using NCLASS for neoclassical transport only, ',
     >       'use_xneo_m =',i2)
 44   format(' Using NCLASS for ExB shear only, ',
     >       'use_xneo_m =',i2)
 55   format(3x,'j',6x,'rho',8x,'chie',11x,'chii',8x,'# NCLASS')
 56   format(3x,'j',6x,'rho',8x,'de',13x,'di',10x,'# NCLASS')
 57   format(3x,'j',6x,'rho',8x,'chie',11x,'chii',8x,'# KAPISN')
 60   format(2x,i2,2x,0p1f10.6,1p6e15.7)
 65   format(/,3x,'j',6x,'rho',8x,'vt-D',11x,'vp-D',
     >        10x,'Er-tot-D',8x,'# NCLASS')
 70   format(/,3x,'j',6x,'rho',8x,'vt-D',11x,'vp-D',
     >        10x,'Er-tot-D',8x,'# NCLASS-interpolated')
c
      end
