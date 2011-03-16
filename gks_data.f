c@gridsetup.f
c jek 18-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... setup transport grid
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine gridsetup(igrid,mxgrid)
c
      implicit none
c      include 'mpif.h'
c      include '../inc/tport.m'
      include 'data_d.m'
      include 'data_exp.m'
c      include '../inc/model.m'
      include 'input.m'
      include 'glf.m'
c      include '../inc/ptor.m'
c
      integer j, aj, igrid, nex, ngg, mxgrid
      real*8 rhox(0:nj-1)
      real, allocatable, dimension(:) :: r_new,
     &   fcap_new, gcap_new, hcap_new, rbp_new, bp0_new,
     &   te_new, ti_new, q_new, ene_new,
     &   en1_new, en2_new, sion_new, 
     &   srecom_new, scx_new, sbcx_new,
     &   s_new, dudtsv_new, enbeam_new,
     &   enn1_new, ennw1_new, ennv1_new,
     &   sbion_new, sbeam_new, 
     &   curden_new, zeff_new, 
     &   angrot_new, chieinv_new, chiinv_new, 
     &   dpedtc_new, dpidtc_new, 
     &   qconde_new, qcondi_new,
     &   qconve_new, qconvi_new, qbeame_new,
     &   qbeami_new, qdelt_new, qrad_new,
     &   qohm_new, qrfe_new, qrfi_new, 
     &   qione_new, qioni_new, qcx_new, 
     &   qfuse_new, qfusi_new, 
     &   rmaj_new, rmin_new, 
     &   psivolp_new, elongx_new, deltax_new, 
     &   sfareanpsi_new, grho1_new,
     &   grho2_new, xkineo_new,
     &   xb2_new, xbm2_new, xngrth_new,
     &   xgrbm2_new, fm1_new, fm2_new,
     &   fm3_new, fhat_new, torque_new,
     &   ptot_new, pfast_new, er_new
      allocate (r_new(nj),
     &   fcap_new(nj), gcap_new(nj), hcap_new(nj), rbp_new(nj),
     &   bp0_new(nj), te_new(nj), ti_new(nj), q_new(nj), ene_new(nj),
     &   en1_new(nj), en2_new(nj), sion_new(nj), 
     &   srecom_new(nj), scx_new(nj), sbcx_new(nj),
     &   s_new(nj), dudtsv_new(nj), enbeam_new(nj),
     &   enn1_new(nj), ennw1_new(nj), ennv1_new(nj),
     &   sbion_new(nj), sbeam_new(nj), 
     &   curden_new(nj), zeff_new(nj), 
     &   angrot_new(nj), chieinv_new(nj), chiinv_new(nj), 
     &   dpedtc_new(nj), dpidtc_new(nj), 
     &   qconde_new(nj), qcondi_new(nj),
     &   qconve_new(nj), qconvi_new(nj), qbeame_new(nj),
     &   qbeami_new(nj), qdelt_new(nj), qrad_new(nj),
     &   qohm_new(nj), qrfe_new(nj), qrfi_new(nj), 
     &   qione_new(nj), qioni_new(nj), qcx_new(nj), 
     &   qfuse_new(nj), qfusi_new(nj), 
     &   rmaj_new(nj), rmin_new(nj), 
     &   psivolp_new(nj), elongx_new(nj), deltax_new(nj), 
     &   sfareanpsi_new(nj), grho1_new(nj),
     &   grho2_new(nj), xkineo_new(nj),
     &   xb2_new(nj), xbm2_new(nj), xngrth_new(nj),
     &   xgrbm2_new(nj), fm1_new(nj), fm2_new(nj),
     &   fm3_new(nj), fhat_new(nj), torque_new(nj),
     &   ptot_new(nj), pfast_new(nj), er_new(nj) )
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c      if(i_proc.eq.0) then
         write(*,10) ' nj_d = ',nj_d
        if(nj_d.gt.nj) write(*,*) 'Warning: nj_d > nj'
c      endif
      jmaxm=mxgrid
c
c... experimental rho_hat grid
c
      do j=0,nj
        rhox(j)=0.D0
      enddo
c
      do j=0,nj_d-1
        aj=j
        rhox(j)=REAL(aj)/REAL(nj_d-1)
      enddo
      rhox(0)=1.D-6
      rhox(nj_d-1)=rhox(nj_d-1)-1.D-6
c
c... setup transport grid
c     igrid = 0 for uniform grid
c     igrid > 0 for nonuniform grid
      if(igrid.eq.0) then
        do j=0,nj_d-1
c          aj=j
          rho(j)=REAL(j)/REAL(nj_d-1)
        enddo
        rho(0)=1.D-6
        rho(jmaxm)=rho(jmaxm)-1.D-6
      else
      write(*,*) ' Using non-unform grid'
        do j=0,40
          aj=j
          rho(j)=0.9D0*aj/40
        enddo
        do j=41,50
          aj=j-40
          rho(j)=0.1D0*aj/(50-40)+0.9D0
        enddo
      endif
c
      do j=1,nj_d
       rho_d(j)=rho(j-1)
      enddo
c
c... Interpolate profiles to new grid if 
c    jmaxm > nj_d-1 OR nonuniform grid
c
      if (jmaxm.gt.(nj_d-1) .or. igrid.gt.0) then
        write(*,25)
        nex=jmaxm+1
        if(igrid.eq.0) then
          do j=1,nex
            aj=REAL(j-1)
            rho_d(j)=aj/REAL(jmaxm)
          enddo
        else
          do j=1,40
            aj=j-1
            rho_d(j)=0.9D0*aj/40
          enddo
          do j=41,51
            aj=j-41
            rho_d(j)=0.1D0*aj/10+0.9D0
          enddo
        endif
        rho_d(1)=1.0D-12
c        do j=0,51
c          write(*,'(i3,3f8.4)') j, rho_d(j), rhox(j), r_d(j)
c        enddo
c
        ngg=nj_d
c       call w_lin_interp_r8(ngg,rhox,r_d,nex,rho_d,r_new,iflag,msg)
        call inter_cspl(ngg,rhox,r_d,nex,rho_d,r_new)
        call inter_cspl(ngg,rho,fcap_d,nex,rho_d,fcap_new)
        call inter_cspl(ngg,rho,gcap_d,nex,rho_d,gcap_new)
        call inter_cspl(ngg,rho,hcap_d,nex,rho_d,hcap_new)
        call inter_cspl(ngg,rho,rbp_d,nex,rho_d,rbp_new)
        call inter_cspl(ngg,rho,bp0_d,nex,rho_d,bp0_new)
        call inter_cspl(ngg,rhox,te_d,nex,rho_d,te_new)
        call inter_cspl(ngg,rhox,ti_d,nex,rho_d,ti_new)
        call inter_cspl(ngg,rhox,q_d,nex,rho_d,q_new)
        call inter_cspl(ngg,rhox,ene_d,nex,rho_d,ene_new)
        call inter_cspl(ngg,rhox,en_d(1:j,1),nex,rho_d,en1_new)
        call inter_cspl(ngg,rhox,en_d(1:j,2),nex,rho_d,en2_new)
        call inter_cspl(ngg,rhox,sion_d(1:j,1),nex,rho_d,sion_new)
        call inter_cspl(ngg,rhox,srecom_d(1:j,1),nex,rho_d,srecom_new)
        call inter_cspl(ngg,rhox,sbcx_d(1:j,1),nex,rho_d,sbcx_new)
        call inter_cspl(ngg,rhox,scx_d(1:j,1),nex,rho_d,scx_new)
        call inter_cspl(ngg,rhox,s_d(1:j,1),nex,rho_d,s_new)
        call inter_cspl(ngg,rhox,dudtsv_d(1:j,1),nex,rho_d,dudtsv_new)
        call inter_cspl(ngg,rhox,enbeam_d,nex,rho_d,enbeam_new)
        call inter_cspl(ngg,rhox,enn_d(1:j,1),nex,rho_d,enn1_new)
        call inter_cspl(ngg,rhox,ennw_d(1:j,1),nex,rho_d,ennw1_new)
        call inter_cspl(ngg,rhox,ennv_d(1:j,1),nex,rho_d,ennv1_new)
        call inter_cspl(ngg,rhox,sbion_d,nex,rho_d,sbion_new)
        call inter_cspl(ngg,rhox,sbeam_d,nex,rho_d,sbeam_new)
        call inter_cspl(ngg,rhox,curden_d,nex,rho_d,curden_new)
        call inter_cspl(ngg,rhox,zeff_d,nex,rho_d,zeff_new)
        call inter_cspl(ngg,rhox,angrot_d,nex,rho_d,angrot_new)
        call inter_cspl(ngg,rhox,chieinv_d,nex,rho_d,chieinv_new)
        call inter_cspl(ngg,rhox,chiinv_d,nex,rho_d,chiinv_new)
        call inter_cspl(ngg,rhox,dpedtc_d,nex,rho_d,dpedtc_new)
        call inter_cspl(ngg,rhox,dpidtc_d,nex,rho_d,dpidtc_new)
        call inter_cspl(ngg,rhox,qconde_d,nex,rho_d,qconde_new)
        call inter_cspl(ngg,rhox,qcondi_d,nex,rho_d,qcondi_new)
        call inter_cspl(ngg,rhox,qconve_d,nex,rho_d,qconve_new)
        call inter_cspl(ngg,rhox,qconvi_d,nex,rho_d,qconvi_new)
        call inter_cspl(ngg,rhox,qbeame_d,nex,rho_d,qbeame_new)
        call inter_cspl(ngg,rhox,qbeami_d,nex,rho_d,qbeami_new)
        call inter_cspl(ngg,rhox,qdelt_d,nex,rho_d,qdelt_new)
        call inter_cspl(ngg,rhox,qrfe_d,nex,rho_d,qrfe_new)
        call inter_cspl(ngg,rhox,qrfi_d,nex,rho_d,qrfi_new)
        call inter_cspl(ngg,rhox,qione_d,nex,rho_d,qione_new)
        call inter_cspl(ngg,rhox,qioni_d,nex,rho_d,qioni_new)
        call inter_cspl(ngg,rhox,qcx_d,nex,rho_d,qcx_new)
        call inter_cspl(ngg,rhox,qfuse_d,nex,rho_d,qfuse_new)
        call inter_cspl(ngg,rhox,qfusi_d,nex,rho_d,qfusi_new)
        call inter_cspl(ngg,rhox,qrad_d,nex,rho_d,qrad_new)
        call inter_cspl(ngg,rhox,qohm_d,nex,rho_d,qohm_new)
        call inter_cspl(ngg,rhox,rmajavnpsi_d,nex,rho_d,rmaj_new)
        call inter_cspl(ngg,rhox,rminavnpsi_d,nex,rho_d,rmin_new)
        call inter_cspl(ngg,rhox,psivolp_d,nex,rho_d,psivolp_new)
        call inter_cspl(ngg,rhox,elongx_d,nex,rho_d,elongx_new)
        call inter_cspl(ngg,rhox,deltax_d,nex,rho_d,deltax_new)
        call inter_cspl(ngg,rhox,sfareanpsi_d,nex,rho_d,sfareanpsi_new)
        call inter_cspl(ngg,rhox,grho1npsi_d,nex,rho_d,grho1_new)
        call inter_cspl(ngg,rhox,grho2npsi_d,nex,rho_d,grho2_new)
        call inter_cspl(ngg,rhox,ptot_d,nex,rho_d,ptot_new)
        call inter_cspl(ngg,rhox,pfast_d,nex,rho_d,pfast_new)
        if (ncl_flag .eq. 1) then
          call inter_cspl(ngg,rhox,xb2_d,nex,rho_d,xb2_new)
          call inter_cspl(ngg,rhox,xbm2_d,nex,rho_d,xbm2_new)
          call inter_cspl(ngg,rhox,xngrth_d,nex,rho_d,xngrth_new)
          call inter_cspl(ngg,rhox,xgrbm2_d,nex,rho_d,xgrbm2_new)
          call inter_cspl(ngg,rhox,fm1_d,nex,rho_d,fm1_new)
          call inter_cspl(ngg,rhox,fm2_d,nex,rho_d,fm2_new)
          call inter_cspl(ngg,rhox,fm3_d,nex,rho_d,fm3_new)
          call inter_cspl(ngg,rhox,fhat_d,nex,rho_d,fhat_new)
        endif
        if(itorque.ne.0) then
          call inter_cspl(ngg,rhox,torque_d,nex,rho_d,torque_new)
        endif
c
c... map interpolated variables back to data (_d) variables
c
c        do j=1,jmaxm+1
c          write(*,50) j, rhox(j-1), r_d(j), ti_d(j), 
c     &                rho_d(j), r_new(j), ti_new(j)
c        enddo
        do j=1,nex
          r_d(j)=r_new(j)
          fcap_d(j)=fcap_new(j)
          gcap_d(j)=gcap_new(j)
          hcap_d(j)=hcap_new(j)
          rbp_d(j)=rbp_new(j)
          bp0_d(j)=bp0_new(j)
          te_d(j)=te_new(j)
          ti_d(j)=ti_new(j)
          q_d(j)=q_new(j)
          ene_d(j)=ene_new(j)
          en_d(j,1)=en1_new(j)
          en_d(j,2)=en2_new(j)
          sion_d(j,1)=sion_new(j)
          srecom_d(j,1)=srecom_new(j)
          sbcx_d(j,1)=sbcx_new(j)
          scx_d(j,1)=scx_new(j)
          s_d(j,1)=s_new(j)
          dudtsv_d(j,1)=dudtsv_new(j)
          enbeam_d(j)=enbeam_new(j)
          enn_d(j,1)=enn1_new(j)
          ennw_d(j,1)=ennw1_new(j)
          ennv_d(j,1)=ennv1_new(j)
          sbion_d(j)=sbion_new(j)
          sbeam_d(j)=sbeam_new(j)
          curden_d(j)=curden_new(j)
          zeff_d(j)=zeff_new(j)
          angrot_d(j)=angrot_new(j)
          chieinv_d(j)=chieinv_new(j)
          chiinv_d(j)=chiinv_new(j)
          dpedtc_d(j)=dpedtc_new(j)
          dpidtc_d(j)=dpidtc_new(j)
          qconde_d(j)=qconde_new(j)
          qcondi_d(j)=qcondi_new(j)
          qconve_d(j)=qconve_new(j)
          qconvi_d(j)=qconvi_new(j)
          qbeame_d(j)=qbeame_new(j)
          qbeami_d(j)=qbeami_new(j)
          qdelt_d(j)=qdelt_new(j)
          qrfe_d(j)=qrfe_new(j)
          qrfi_d(j)=qrfi_new(j)
          qione_d(j)=qione_new(j)
          qioni_d(j)=qioni_new(j)
          qcx_d(j)=qcx_new(j)
          qfuse_d(j)=qfuse_new(j)
          qfusi_d(j)=qfusi_new(j)
          qrad_d(j)=qrad_new(j)
          qohm_d(j)=qohm_new(j)
          rmajavnpsi_d(j)=rmaj_new(j)
          rminavnpsi_d(j)=rmin_new(j)
          psivolp_d(j)=psivolp_new(j)
          elongx_d(j)=elongx_new(j)
          deltax_d(j)=deltax_new(j)
          sfareanpsi_d(j)=sfareanpsi_new(j)
          grho1npsi_d(j)=grho1_new(j)
          grho2npsi_d(j)=grho2_new(j)
          torque_d(j)=torque_new(j)
          ptot_d(j)=ptot_new(j)
          pfast_d(j)=pfast_new(j)
c          er_d(j)=er_new(j)
          if(itorque.eq.0) torque_d(j)=0.D0
          if(iptotr.eq.0) pfast_d(j)=0.D0
        enddo
        if (ncl_flag .eq. 1) then
          do j=1,nex
            xb2_d(j)=xb2_new(j)
            xbm2_d(j)=xbm2_new(j)
            xngrth_d(j)=xngrth_new(j)
            xgrbm2_d(j)=xgrbm2_new(j)
            fm1_d(j)=fm1_new(j)
            fm2_d(j)=fm2_new(j)
            fm3_d(j)=fm3_new(j)
            fhat_d(j)=fhat_new(j)
          enddo
        endif
c
c deallocate temp arrays
c
      deallocate( r_new, 
     &   fcap_new, gcap_new, hcap_new, rbp_new, bp0_new,
     &   te_new, ti_new, q_new, ene_new,
     &   en1_new, en2_new, sion_new, 
     &   srecom_new, scx_new, sbcx_new,
     &   s_new, dudtsv_new, enbeam_new,
     &   enn1_new, ennw1_new, ennv1_new,
     &   sbion_new, sbeam_new, 
     &   curden_new, zeff_new, 
     &   angrot_new, chieinv_new, chiinv_new, 
     &   dpedtc_new, dpidtc_new, 
     &   qconde_new, qcondi_new,
     &   qconve_new, qconvi_new, qbeame_new,
     &   qbeami_new, qdelt_new, qrad_new,
     &   qohm_new, qrfe_new, qrfi_new, 
     &   qione_new, qioni_new, qcx_new, 
     &   qfuse_new, qfusi_new, 
     &   rmaj_new, rmin_new, 
     &   psivolp_new, elongx_new, deltax_new, 
     &   sfareanpsi_new, grho1_new,
     &   grho2_new, xkineo_new,
     &   xb2_new, xbm2_new, xngrth_new,
     &   xgrbm2_new, fm1_new, fm2_new,
     &   fm3_new, fhat_new, torque_new,
     &   ptot_new, pfast_new, er_new)
c
c... reset rho grid after interpolation
c
      if(igrid.eq.0) then
        do j=0,jmaxm
          aj=DFLOAT(j)
          rho(j)=aj/DFLOAT(jmaxm)
        enddo
      else
        do j=0,40
          aj=j
          rho(j)=0.9D0*aj/40
        enddo
        do j=41,50
          aj=j-40
          rho(j)=0.1D0*aj/(50-40)+0.9D0
        enddo
      endif
      rho(0)=1.D-6
      rho(jmaxm)=rho(jmaxm)-1.D-6
c
cjek nj_d -> mxgrid+1 for arbitrary grid size
c
      nj_d=mxgrid+1
c
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 10   format(a9,i3)             ! common write/read  format
 25   format(' Interpolating ONETWO iterdb grid to new grid ...')
 50   format(i2,2x,0p1f4.2,0p6f10.5)
c
      end
c@u2tread.f   derived from u2read.f by G. M. Staebler 6/28/2000
c reads full 2d time series from u-file.c@readiterdb.f
c jek 19-Jan-11 version 2.0
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Reads in data from Onetwo ascii iterdb file
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine readiterdb(mxgrid)
c
      implicit none
      include 'input.m'
      include 'data_d.m'
      include 'data_exp.m'
      include 'glf.m'
c
      integer iflag
      character*15 extension
      character*50 iterdbfile
      character cdfile*130
c
      integer niterdb, i, j, jj, jn, mxgrid
      integer idchar, ifchar, ichar, ichmax
      character*132 headerline
      character*2 stflg
c
      real*8 pi_m, pgasa, pgasz, bgasa, bgasz
     & , rhostar, drm, dvoldr_p, dvoldr_m, dummy1, aomega
c
      parameter(niterdb=11)
      namelist/datafiles/nstk,xp_time,cudir,tok,shot,phase,ismooth
c      real, allocatable, dimension(:) :: blank_d, bblank_d
c      allocate (blank_d(nj), bblank_d(1000))
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      pi_m=3.1415926D0
c      xi=(0.D0,1.D0)
      kevdsecpmw=1.6022D-19*1.D3*1.D-6
      jmaxm=mxgrid
c
c     Assume the worki and beam ions are same and the impurity
c     is carbon for the moment.
c
      pgasa=amassgas_exp
      pgasz=1
      bgasa=pgasa
      bgasz=1
c     pimpa=12
c     pimpz=6
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... read ONETWO iterdb file from external directory
c    JEK 3/19/07
c
      extension = shot
      iterdbfile = 'iterdb.'//extension
c      write(*,*) 'iterdbfile = ',iterdbfile
c
c... count characters in directory name -> idchar
c
      ichmax = 60
      idchar = 0
      do j=1,ichmax
        if ( cudir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 5
        endif
      enddo
 5    continue
c
c... count characters in iterdb filename -> ifchar
c
      ifchar = 0
      do j=1,ichmax
        if ( iterdbfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 6
        endif
      enddo
  6   continue
c
c       write (*,*) 'idchar = ',idchar
c       write (*,*) 'ifchar = ',ifchar
       if (    cudir(idchar:idchar) .ne. '/') then
          ichar = idchar + 1 + ifchar
          cdfile = cudir(1:idchar) // '/' // iterdbfile(1:ifchar)
       else
          ichar  = idchar + ifchar
          cdfile = cudir(1:idchar) // iterdbfile(1:ifchar)
       endif
       write (*,*) ' iterdb = ',cdfile(1:ichar)
c
       open(unit=niterdb,status='old',access='sequential',
     &      file=cdfile(1:ichar))
c old coding to read in same dir as xptor
c     &      file='iterdb.'//extension)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
         read(niterdb,'(a)')headerline  

c ishot : shot number
         read(niterdb,'(a)') stflg
         read(niterdb,9) ishot_d

c nj : the size of the vectors printed in this file
         read(niterdb,'(a)') stflg
         read(niterdb,9) nj_d
         if(nj_d.gt.nj+1)then
           write(*,900) nj_d,nj+1
           stop
         else
           jmaxm=nj_d-1
           mxgrid=jmaxm
         endif

c nion : the number of ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nion_d

c nprim : the number of primary ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nprim_d

c nimp : the number of impurity ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nimp_d

c nneu : the number of neutral ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nneu_d

c ibion : index of beam species
         read(niterdb,'(a)')stflg
         read(niterdb,9)ibion_d

c namep : name(s) of primary ion species
         read(niterdb,'(a)')stflg
         read(niterdb,8)(namep_d(i),i=1,nprim_d)

c namei : name(s) of impurity ion species
         read(niterdb,'(a)')stflg
         read(niterdb,8)(namei_d(i),i=1,nimp_d)

c namen : name(s) of neutral ion species
         read(niterdb,'(a)')stflg
         read(niterdb,8)(namen_d(i),i=1,nneu_d)

c time :  time at which data is printed
         read(niterdb,'(a)')stflg
         read(niterdb,10)time_d

c Rgeom : major radius of geometric
         read(niterdb,'(a)')stflg
         read(niterdb,10)rgeom_d                   !meters

c Rmag :  major radius of mag axis, meters
         read(niterdb,'(a)')stflg
         read(niterdb,10)rmag_d                    !meters

c R0 : major radius of vaccuum btor ref
         read(niterdb,'(a)')stflg
         read(niterdb,10)rmajor_d                  !meters

c kappa : plasma elongation
         read(niterdb,'(a)')stflg
         read(niterdb,10)kappa_d

c delta : plasma triangularity
         read(niterdb,'(a)')stflg
         read(niterdb,10)deltao_d

c pindent : plasma indentation
         read(niterdb,'(a)')stflg
         read(niterdb,10)pindento_d

c volo : plasma volume,meters**3
         read(niterdb,'(a)')stflg
         read(niterdb,10)volo_d                    !meters**3

c cxareao :plasma cross sectional area, meters**2
         read(niterdb,'(a)')stflg
         read(niterdb,10)areao_d                   !meters**2

c Btor : vaccuum toroidal field at rmajor, tesla
         read(niterdb,'(a)')stflg
         read(niterdb,10)btor_d                    !tesla

c total,ohmic,bootstrap,beam,and rf currents, amps
         read(niterdb,'(a)')stflg
         read(niterdb,10)tocur_d,totohm_d,totboot_d,totbeam_d,totrf_d !amps

c betap : poloidal beta
         read(niterdb,'(a)')stflg
         read(niterdb,10)betap_d

c beta : toroidal beta
         read(niterdb,'(a)')stflg
         read(niterdb,10)beta_d

c ali : plasma inductance
         read(niterdb,'(a)')stflg
         read(niterdb,10)ali_d

c te0 : central electron temperature
         read(niterdb,'(a)')stflg
         read(niterdb,10)te0_d                        !kev

c ti0 : central ion temperature
         read(niterdb,'(a)')stflg
         read(niterdb,10)ti0_d                        !kev

c---psi on rho grid,volt*sec/rad

          read(niterdb,'(a)')stflg
          read(niterdb,10)(psir_d(j),j=1,nj_d)      !volt*se/rad
cx          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---rho grid, meters

          read(niterdb,'(a)')stflg
          read(niterdb,10)(r_d(j),j=1,nj_d)          !meters

c---fcap, (ie f(psilim)/f(psi) )

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fcap_d(j),j=1,nj_d)
cx          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---gcap, (ie <(grad rho)**2*(R0/R)**2> )

          read(niterdb,'(a)')stflg
          read(niterdb,10)(gcap_d(j),j=1,nj_d)
cx          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---hcap, (ie (dvolume/drho)/(4*pi*pi*R0*rho))

          read(niterdb,'(a)')stflg
          read(niterdb,10)(hcap_d(j),j=1,nj_d)

c*** read in NCLASS quantities if ncl_flag=1 ***
          if(ncl_flag.eq.1) then
            write(*,25)

c--xb2, <B^2>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xb2_d(j),j=1,nj_d)

c--xbm2, <1/B^2>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xbm2_d(j),j=1,nj_d)

c--xngrth, <n.grad theta>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xngrth_d(j),j=1,nj_d)

c--xgrbm2, <grad-rho^2/B^2>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xgrbm2_d(j),j=1,nj_d)

c--fm1, 1st poloidal moment

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fm1_d(j),j=1,nj_d)

c--fm2, 2nd poloidal moment

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fm2_d(j),j=1,nj_d)

c--fm3, 3rd poloidal moment

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fm3_d(j),j=1,nj_d)

c--fhat, muo*F/dpsidr

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fhat_d(j),j=1,nj_d)
c
          endif
c
c---electron temperature, kev

          read(niterdb,'(a)')stflg
          read(niterdb,10)(te_d(j),j=1,nj_d)           !kev

c---ion temperature, kev

          read(niterdb,'(a)')stflg
          read(niterdb,10)(ti_d(j),j=1,nj_d)           !kev
c
c---q (ie safety factor) profile

          read(niterdb,'(a)')stflg
          read(niterdb,10)(q_d(j),j=1,nj_d)
c
c---electron density,#/m**3 

          read(niterdb,'(a)')stflg
          read(niterdb,10)(ene_d(j),j=1,nj_d)           !#/meter**3

c---primary ion density,#/m**3,species
c---impurity ion density,#/m**3,species

       do jj=1,nion_d
         do j=1,nj_d
           en_d(j,jj)=0.D0
         enddo
       enddo
c      jp=0
c      ji=0
       do jj=1,nion_d
c          if(jj .le. nprim_d)jp=jp+1
c          if(jj .gt. nprim_d)ji=ji+1
               read(niterdb,'(a)')stflg
               read(niterdb,10)(en_d(j,jj),j=1,nj_d)    !#/meter**3
       enddo

c---sion   : source due to ionization,        #/(m**3*sec),species
c---srecom : source due to recombination,     #/(m**3*sec),species
c---scx    : source due to cx thermal neut.,  #/(m**3*sec),species
c---sbcx   : sink due to cx with beam neut.   #/(m**3*sec),species
c---s      : total source rate,               #/(m**3*sec),species
c---dudt   : s dot,                           #/(m**3*sec),species
 
       do jj=1,nprim_d
               read(niterdb,'(a)')stflg
               read(niterdb,10)(sion_d(j,jj),j=1,nj_d)    !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(srecom_d(j,jj),j=1,nj_d)  !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(scx_d(j,jj),j=1,nj_d)     !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(sbcx_d(j,jj),j=1,nj_d)    !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(s_d(j,jj),j=1,nj_d)       !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(dudtsv_d(j,jj),j=1,nj_d)  !#/meter**3/sec
          enddo

c---fast ion density, #/m**3, species
       
          read(niterdb,'(a)')stflg  
          read(niterdb,10)(enbeam_d(j),j=1,nj_d)          !#/meter**3
       
c---neutral density,  #/m**3,species

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(enn_d(j,jn),j=1,nj_d)     !#/meter**3
       enddo

c---neutral density from wall source, #/m**3, species

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(ennw_d(j,jn),j=1,nj_d)     !#/meter**3
       enddo

c---neutral density from volume source, #/m**3,species (recomb and beam cx)

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(ennv_d(j,jn),j=1,nj_d)     !#/meter**3
       enddo

c---volume source of neutrals, #/(m**3*sec),species

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(volsn_d(j,jn),j=1,nj_d)!#/(m**3*sec)  
       enddo

c---sbion : beam electron source,  #/(m**3*sec
           
          read(niterdb,'(a)')stflg
          read(niterdb,10)(sbion_d(j),j=1,nj_d)    !#/(m**3*sec)

c---sbeam : beam thermal ion source,  #/(m**3*sec)
           
          read(niterdb,'(a)')stflg   
          read(niterdb,10)(sbeam_d(j),j=1,nj_d)      !#/(m**3*sec)

c---total current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curden_d(j),j=1,nj_d)        !amps/meter**2
      
c---ohmic current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curohm_d(j),j=1,nj_d)        !amps/meter**2

c---bootstrap current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curboot_d(j),j=1,nj_d)       !amps/meter**2
       
c--- beam driven current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curbeam_d(j),j=1,nj_d)      !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
                
c---rf current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(currf_d(j),j=1,nj_d)         !amps/meter**2

c---rho*bp0*fcap*gcap*hcap, tesla*meters
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(rbp_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)         !tesla*meters      

c---zeff profile:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(zeff_d(j),j=1,nj_d)

c---angular rotation speed profile, rad/sec
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(angrot_d(j),j=1,nj_d)       !rad/sec

c---electron thermal diffusivity, meters*2*/sec on half grid
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(chieinv_d(j),j=1,nj_d)       !meters**2/sec
          
c---ion thermal diffusivity,meters*2*/sec  on half grid 
      
          read(niterdb,'(a)')stflg 
          read(niterdb,10)(chiinv_d(j),j=1,nj_d)        !meters**2/sec

c---ion neocl.  thermal conductivity, 1/(m*sec) on half grid
          
          read(niterdb,'(a)')stflg
          read(niterdb,10)(xkineo_d(j),j=1,nj_d)          !1/(m*sec)

c---wdot,electrons, watts/m**3 d(electron energy)/dt profile:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(dpedtc_d(j),j=1,nj_d)        !watts/meter**3
     
c---wdot,ions,watts/m**3 d(ion energy)/dt profile:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(dpidtc_d(j),j=1,nj_d)        !watts/meter**3

c---electron conduction, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qconde_d(j),j=1,nj_d)        !watts/meter**3
      
c---ion conduction, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qcondi_d(j),j=1,nj_d)        !watts/meter**3

c---electron convection,watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qconve_d(j),j=1,nj_d)        !watts/meter**3

c---ion convection, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qconvi_d(j),j=1,nj_d)        !watts/meter**3

c---power to elec.from beam, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qbeame_d(j),j=1,nj_d)        !watts/meter**3

c---qdelt : electron-ion equilibration, watts/m**3
      
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qdelt_d(j),j=1,nj_d)         !watts/meter**3
c
c---power to ions from beam, watts/m**3:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qbeami_d(j),j=1,nj_d)        !watts/meter**3
     
c---qrfe, rf electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qrfe_d(j),j=1,nj_d)          !watts/meter**3
      
c---qrfi, rf ion heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qrfi_d(j),j=1,nj_d)          !watts/meter**3

c---qione, recombination and impact ionization, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qione_d(j),j=1,nj_d)         !watts/meter**3
       
c---qioni, recombination and impact ionization, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qioni_d(j),j=1,nj_d)         !watts/meter**3
       
c---qcx,  neutral-ion charge exchange, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qcx_d(j),j=1,nj_d)           !watts/meter**3

c---2d MHD equil. electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qe2d_d(j),j=1,nj_d)        !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
          
c---2d MHD equil.ion heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qi2d_d(j),j=1,nj_d)         !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
          
c---fusion electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qfuse_d(j),j=1,nj_d)        !watts/meter**3
                
c---fusion ion heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qfusi_d(j),j=1,nj_d)        !watts/meter**3

c---beam fusion electron heating,watts/m**3
    
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfuse_d(j),j=1,nj_d)       !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                 
c---beam fusion ion heating profile
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfusi_d(j),j=1,nj_d)       !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                 
c---qmag electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qmag_d(j),j=1,nj_d)         !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                
c---sawtooth electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawe_d(j),j=1,nj_d)        !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                
c---sawtooth ion  heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawi_d(j),j=1,nj_d)        !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---radiated power density, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qrad_d(j),j=1,nj_d)           !watts/meter**3

c---(electron) ohmic power density,watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qohm_d(j),j=1,nj_d)           !watts/meter**3

c---avg major radius of each flux surface, meters, at elevation of mag. axis
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(rmajavnpsi_d(j),j=1,nj_d)   !meters

c---avg minor radius of each flux surface,meters, at elevation of mag. axis

          read(niterdb,'(a)')stflg
          read(niterdb,10)(rminavnpsi_d(j),j=1,nj_d)   !meters
c make gradient smooth near axis so that drhodr will be smooth
          amin_d=rminavnpsi_d(nj_d)   !meters
          
c---volume of each flux surface,  meters**3
                 
          read(niterdb,'(a)')stflg
          read(niterdb,10)(psivolp_d(j),j=1,nj_d)      !meters**3
c
c--- elongation of each flux surface
c
          read(niterdb,'(a)')stflg
          read(niterdb,10)(elongx_d(j),j=1,nj_d)
c
c---triangularity of each flux surface
c
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(triangnpsi_d(j),j=1,nj_d)
          read(niterdb,10)(deltax_d(j),j=1,nj_d)
c
c---indentation of each flux surface
c 
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(pindentnpsi(j),j=1,nj_d)
          read(niterdb,10)(blank_d(j),j=1,nj_d)
c
c--- flux surface area,  meters**2 4.*pi*pi*R0*hcap*rho*<abs(grad rho)>
c
          read(niterdb,'(a)')stflg
          read(niterdb,'(a)')stflg
          read(niterdb,10)(sfareanpsi_d(j),j=1,nj_d)  !meters**2
c
c---cross sectional area of each flux surface, meters**2
c    
          read(niterdb,'(a)')stflg
          read(niterdb,10)(cxareanpsi_d(j),j=1,nj_d)  !meters**2
c      
c---flux surface average absolute grad rho
c     
         read(niterdb,'(a)')stflg
         read(niterdb,10)(grho1npsi_d(j),j=1,nj_d)
c
c---flux surface average ( grad rho)**2
c 
         read(niterdb,'(a)')stflg
         read(niterdb,10)(grho2npsi_d(j),j=1,nj_d)
c
c---plasma boundary:
c            nplasbdry : no. pts on plasma bdry
c            r pts for plasma boundary, meters
c            z pts for plasma boundary, meters
c
             read(niterdb,'(a)')stflg
             read(niterdb,9)nplasbdry_d
c
             read(niterdb,'(a)')stflg
cx           read(niterdb,10)(rplasbdry_d(j), j=1,nplasbdry_d)
             read(niterdb,10)(bblank_d(j),j=1,nplasbdry_d)
c
             read(niterdb,'(a)')stflg
cx           read(niterdb,10)(zplasbdry_d(j), j=1,nplasbdry_d)
             read(niterdb,10)(bblank_d(j),j=1,nplasbdry_d)
c
c---torque density  Nt-m/m^3 (old iterdb files were in dyne-cm/cm^3 = 10 Nt-m/m^3)
c
             if(itorque.ne.0)then
               read(niterdb,'(a)')stflg
               read(niterdb,10)(torque_d(j), j=1,nj_d)
c               write(*,10) (torque_d(j),j=1,nj_d)
c this rescaling is not needed anymore. 
c For old iterdb files in cgs units set s0(4)=0.1, s0(5)=0.1, s1(4)=0.1, s1(5)=0.1 
c to rescale the torque. 
c                if(idata.eq.1) then
c                 write(*,*) 'Rescaling torque density ...'
c                 do j=1,nj_d
c                  torque_d(j) = torque_d(j) / 10.D0  ! convert to N-M/M**3
c                 enddo
c                endif
             endif
c
c---total and fast ion pressure (keV/m^3)
c   Pa = N/m**2, keV/m**3=Pa/1.602e-16
c   If iptotr=0, then total thermal pressure computed in pressure.f
c   Note: pfast_exp has 1.e19 factored out since alpha_exp,m
c   has 1.e19 factored out of densities
c
       if(iptotr.eq.0)then
         write(*,'(a32)')
     >   'computing total pressure from thermal pressure'
         do j=1,nj_d
           ptot_d(j) = 1.6022D-16*ene_d(j)*(te_d(j)+ti_d(j))   ! Pascals
           pfast_d(j) = 0.0
         enddo
       elseif(iptotr.eq.1) then
         write(*,'(a32)') 'Reading ptot only' 
         read(niterdb,'(a)')stflg
         read(niterdb,10)(ptot_d(j), j=1,nj_d)
         do j=1,nj_d
           pfast_d(j)=ptot_d(j)
     >     -1.6022D-16*ene_d(j)*(te_d(j)+ti_d(j))  !Pascals
         enddo
       elseif(iptotr.eq.2) then
         write(*,'(a32)') 'Reading ptot and pfast'
         read(niterdb,'(a)')stflg
         read(niterdb,10)(pfast_d(j), j=1,nj_d)
         read(niterdb,'(a)')stflg
         read(niterdb,10)(ptot_d(j), j=1,nj_d)
       endif
c
       close(unit=niterdb,status='keep')
c      close(unit=niterdb)
c
      do j=1,nj_d
       bp0_d(j)=rbp_d(j)/( max(r_d(j),1.d-6)*fcap_d(j)*
     >          gcap_d(j)*hcap_d(j) )
c       write(*,50) j,rho_d(j), r_d(j), fcap_d(j),gcap_d(j),
c     >            hcap_d(j), bp0_d(j)
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
  8   format(5(2x,a))      !common character variable write/read format
  9   format(5(2x,i6))     ! common integer write/read format
 10   format(5(2x,1pe14.4))! common write/read  format
 12   format(a9,i3)! common write/read  format
 25   format(' Reading geometric variables for NCLASS ...')
 50   format(i2,2x,0p1f4.2,0p6f10.5)
 51   format(i3,2x,0p1f6.4,1p6e12.4)
 52   format(i2,2x,i2,2x,0p1f4.2,0p6f10.5)
 53   format(i2,2x,0p1f4.2,0p6e13.5)
 54   format(i2,2x,0p6e13.5)
 85   format(' Smoothing electron density profile ...',0p1f4.2)
 86   format(' Smoothing ion density profile ...',0p1f4.2)
 87   format(' Smoothing Zeff profile ...',0p1f4.2)
 88   format(' Smoothing q-profile ...',0p1f4.2)
 89   format(' Smoothing angrot profile ...',0p1f4.2)
 90   format(' Smoothing fast ion density profile ...',0p1f4.2)
 91   format(' Smoothing triangularity profile ...',0p1f4.2)
 95   format(' Smoothing Qnb profile ...',0p1f4.2)
 100  format(' Error opening iterdb file')
 900  format(' Warning: nj_d larger than maximum pts, nj_d = ',i3,
     > 'max=',i3)
c
       end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine trapl (r, y, nj, xint)
c
      implicit integer (i-n), real*8 (a-h, o-z)
c
c this subroutine integrates y(r) with respect to r from zero to rminor
c the trapezoidal rule is used
c
      dimension  r(*), y(*)
c
      xint = 0.0
      do 10 j=2,nj
 10   xint = xint + 0.5 * (y(j)+y(j-1)) * (r(j)-r(j-1))
      return
c
      end
      subroutine readufiles(tok,shot,phase,cudir,mxgrid,ismooth,
     &                    btscale,xp_time,endtime_pt,time_series,
     &                    iexp_exch,iexp_q,
     &                    p_glob_exp,gtauth_exp,gtautot_exp,iproc)
c************************************************************************
c Reads experimental profiles and splines from ITER PDB 
c 1D and 2D Ufiles.   
c
c The following variables denote the mass and charge of the
c main, beam, and impurity ion species:
c     pgasa      plasma working gas atomic number
c     pgasz      plasma working gas charge
c     bgasa      beam species atomic number
c     bgasz      beam species charge
c     pimpa      impurity species atomic number
c     pimpz      impurity species charge
c 
c Global quantities vs. time in 1D ufiles:
c
c     amin       horizontal minor radius (m)
c     bt         vacuum toroidal field (T) at RGEO
c     delta      triangularity 
c     indent     indentation 
c     ip         plasma current (A)
c     kappa      elongation
c     li         induction 
c     q95        plasma safety factor evaluated at the flux surface that 
c                encloses 95% of the total poloidal flux
c     rgeo       plasma geometrical major radius (m)
c     wtot       total stored energy
c     zeff       ave zeff
c
c Profile quantities vs. rho, time in 2D ufiles:
c
c     chie         experimentally inferred chi_e
c     chii         experimentally inferred chi_i
c     curnbi       NBI currents??
c     curtot       Total current??
c     grho1        < abs(grad rho)>
c     grho2        <(grad rho)**2>
c     ne           electron density
c     nfast        fast ion density
c     nimp         impurity density
c     nm1          main ion density
c     qnbie        power from NBI to electrons
c     qnbii        power from NBI to ions
c     qohm         ohmic heating power
c     qrad         radiated power
c     rmajor       avg major radius of flux surface
c     rminor       avg minor radius of flux surface
c     snbie        electron particle source from NBI
c     snbii        ion particle source from NBI
c     swall        thermal ion particle source due to ionisation
c     surf         flux surface surface area
c     te           T_e
c     ti           T_i
c     volume       flux surface volume
c     zeffr        Z_eff profile
c     ptot         total pressure profile
c
c The data is assumed to be in directory cudir with a name of the form:
c machineNdshotnum.FIELD
c
c The ufiles are assumed to have names beginning with the tokamak,
c then 1d or 2d then the shot number. For example,
c
c /u/usernumber/iterpdb/d3d/69627/d3d1d69627.AMIN
c /u/usernumber/iterpdb/d3d/69627/d3d2d69627.TI
c
c A number of global quantites (e.g. Wtot) are often provided in a
c single 0D file. For example,
c
c /u/usernumber/iterpdb/d3d/69627/d3d_69627_0d.dat
c
c The 0D file is not in ufile format. It is in either in column format
c (idatzero=0) or using a column delimited format (idatzero=1, 2, or 3).
c
c If a file is not present, the relevant variables are either not
c set or some default value (e.g. zeff=1.0)
c
c nrmax = max number of radial grid points in experimental data
c         changed array declarations using nrmax -> nj
c nr   = actual number of radial grid points in experimental data
c rho_xp  = radial variable read from experimental data

c u3dx = radial positions of measured quantities like TEXP, TEXBR, etc.
c u3dy = values of measured quantities like TEXP, TEXBR, etc.
c nx3  = number of points in u3dx and u3dy.
c
c************************************************************************
      implicit none
c
      include 'input.m'
      include 'data_d.m'
c      include 'data_exp.m'
c
      integer iproc, mxgrid, idat, nrmax, ntmax, n1d, n2d, n3d, nur
c
      parameter(nrmax=301) !max number of radial points in experimental grid
      parameter(ntmax=3800) !max number of time points in experimental grid
      parameter(n1d=13)   !number of 1d fields read;
      parameter(n2d=48) 
      parameter(n3d=6)   !number of fields read for experimental data
c      include '../inc/ut.m'  ! need ntmax defined first
c
      character*6 u0phase
      character*7 fields1d(n1d), fields2d(n2d), fields3d(n3d)
      character*50 ufile, erfile
      character*100 udfile
      character cdfile*130
      character*67 u0names(200)
      character*11 u0values(200)
      character*(40) shot
      character*(10) tok
      character*(6) phase
      character*50 cudir
c
      integer i, ig, is, i2d, ierr, time_channel
     &  , j, i1d, id, ifchar, ichar, ifail, k, jstart, jend
     &  , nnexp_d, ntexp_d, ntixp_d, jj, jn, niterdb, linblnk, ier
     &  , ichmax, nr_er, iread, nsmooth, time_flag, itt, itime_shift
      integer mxnr_r
      integer nx3(n3d)
      parameter(mxnr_r=300)
c
      real*8 rfix, sfix, q01_exp, qa1_exp, zeff_t, zbrac, alamda
     &  , tau_e, ds, z2, alfj1_exp, xxq, alpha, sum, cur0, cursum
     &  , rhostar, shift_shaf, dummy1, aomega
     &  , zhalf, zthird, z2thrd, z4thrd, z3qtr
      real*8 xp_time, endtime_pt, btscale, gtauth_exp,gtautot_exp
c
      real*8 rgeo_0d, rmin_0d, bt_0d, ip_0d, nebar_0d, pnbi_0d
     &  , ti0_0d, te0_0d, pich_0d, pech_0d, q0_0d, q95_0d
     &  , kap_0d, delt_0d, wth_0d, wtot_0d
     &  , rho_norm_loc, aj
      real*8 pi_m, kevdsecpmw, p_glob_exp
      real*8 lnlam,eta,pohsum,pradsum,drm,dvoldr_p,dvoldr_m
     &  , prfesum, prfisum, conv_torq
      real*8 apgasa,apgasz,abgasa,abgasz,apimpa,apimpz
c
      real*8 xnexp_d(nrmax), ynexp_d(nrmax)
     &  , xnexpeb_d(nrmax), ynexpeb_d(nrmax)
     &  , xtexp_d(nrmax), ytexp_d(nrmax)
     &  , xtexpeb_d(nrmax), ytexpeb_d(nrmax)
     &  , xtixp_d(nrmax), ytixp_d(nrmax)
     &  , xtixpeb_d(nrmax), ytixpeb_d(nrmax)
c
      integer nr, nut, time_series, ismooth, iexp_exch, iexp_q
      real*8 rho_xp(nrmax+1), ut(ntmax)
      real*8 u1d_t(ntmax,n1d), u2d_t(nrmax,ntmax,n2d)
      real*8 rhox(0:nrmax-1)
c
      real*8 xtime_d(ntmax), ytime_d(nrmax,ntmax)
      real*8 rho_er_d(nrmax), er_raw_d(nrmax)
      real*8 ur(nrmax),rho_u(nrmax)
      real*8 pxp(nrmax),par(nrmax,501)
      real*8 u1d(n1d),u2d(nrmax,n2d),u1d_int(ntmax)
      real*8 u3dx(nrmax,n3d),u3dy(nrmax,n3d)
c
      data fields1d/'AMIN   ','BT     ','DELTA  ','INDENT ',
     .	'IP     ','KAPPA  ','LI     ','Q95    ','RGEO   ',
     .  'WTOT   ','ZEFF   ','POHM   ','VSURF  '/
c
      data fields2d/'CHIE   ','CHII   ','CURNBI ','CURTOT ',
     .	'GRHO1  ','GRHO2  ','NE     ','NFAST  ','NIMP   ','NM1    ',
     .	'QNBIE  ','QNBII  ','QOHM   ','QRAD   ','RMAJOR ','RMINOR ',
     .	'SNBIE  ','SNBII  ','SURF   ','ZEFFR  ','VOLUME ','TI     ',
     .	'TE     ','DELTAR ','INDENTR','KAPPAR ','Q      ','QEI    ',
     .	'VROT   ','SWALL  ','DWER   ','DWIR   ','QECHE ','QICRHE ',
     .	'QECHI  ','QICRHI ','QWALLE ','QWALLI ','QFUSE  ','QFUSI  ',
     .	'DNER','NM2 ','NM3 ','PTOTR','QLHE ','TORQ',
     .  'CHIENEO','CHIINEO'/
c
      data fields3d/'NEXP   ','NEXPEB ',
     >  'TEXP   ','TEXPEB ','TIXP   ','TIXPEB '/
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
       do i=1,mxgrid+1
         do j=1,n2d
            u2d(i,j)=0.D0
         enddo
       enddo
c
c...Misc switches
c
      idat=0
      rho_norm_loc=1.D0
      conv_torq=1.0
      nsmooth = ismooth
      jmaxm=mxgrid
      pi_m=3.1415926D0
      kevdsecpmw=1.6022D-19*1.D3*1.D-6
      zhalf = 1.D0/2.D0
      zthird = 1.D0/3.D0
      z2thrd = 2.D0/3.D0
      z4thrd = 4.D0/3.D0
      z3qtr = 3.D0/4.D0
c
c... Setup grids and constants
c
      nr=mxgrid+1
      do j=0,jmaxm
        aj=dfloat(j)
        rhox(j)=aj/dfloat(jmaxm)
        rho_d(j+1)=rhox(j)
      enddo
      rhox(0)=1.D-6
      rhox(jmaxm)=rhox(jmaxm)-1.D-6
      rho_d(1)=rhox(0)
      rho_d(jmaxm+1)=rhox(jmaxm)
c
c... Special shots needing different format reads
c
      is=linblnk(shot,10)
      if(shot(1:is).eq.'35156') idat=0
      if(shot(1:is).eq.'35171') idat=0
      if(shot(1:is).eq.'37379') idat=0
      if(shot(1:is).eq.'37944') idat=0
      if(shot(1:is).eq.'38285') idat=0
      if(shot(1:is).eq.'38287') idat=0
      if(shot(1:is).eq.'46664') idat=0
      if(shot(1:is).eq.'52009') idat=0
      if(shot(1:is).eq.'52014') idat=1
      if(shot(1:is).eq.'52022') idat=0
      if(shot(1:is).eq.'52025') idat=1
      if(shot(1:is).eq.'52096') idat=0
      if(shot(1:is).eq.'50844') idat=1
      if(shot(1:is).eq.'53299') idat=1
      if(shot(1:is).eq.'52014') idat=0
      if(shot(1:is).eq.'52015') idat=0
      if(shot(1:is).eq.'52979' .and. idatzero.eq.0) idat=1
      if(shot(1:is).eq.'55935') idat=1
      if(shot(1:is).eq.'57987') idat=0
      if(shot(1:is).eq.'58159') idat=0
      if(shot(1:is).eq.'58323') idat=0
      if(shot(1:is).eq.'60933') idat=0
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Loop over the 2d files first:
c
      do i=1,n2d
c
c       compose the file name:
         ig=len_trim(tok)
         is=len_trim(shot)         
c         i2d=linblnk(fields2d(i),7)
         i2d=len_trim(fields2d(i))
         ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields2d(i)(1:i2d)
c
         ierr=0
         time_flag=0
         if (i.eq.time_channel) then
           time_flag=1
           if(iproc.eq.0) write(6,*) 'Reading timedata for channel=',i
         endif
         if(time_series.eq.0)then
c          write(*,*) 'ufile = ',ufile
c          write(*,*) 'u2d = ',u2d(1,i)
          call u2read(idat,cudir,ufile,88,rho_xp,nr,xp_time,u2d(1,i),
     &               ierr,ntmax,xtime_d, ytime_d, ntime_d, time_flag,
     &               itime_shift,nsmooth,islice_d,i,iproc)
         else
          if(i.eq.1 .and. iproc.eq.0) write(6,60)
          do k=1,mxgrid+1
            pxp(k) = 0.D0
            do j=1,ntmax
              u2d_t(k,j,i) = 0.D0
            enddo
          enddo
          call u2tread(cudir,ufile,88,nur,nut,ur,ut,par,nsmooth,i,ierr)
          if(ierr.eq.0)then
           if(nut.lt.2)then
             if(iproc.eq.0)write(*,*)"error nut < 2 ",nut
             return
           endif
           if(ut(1).gt.endtime_pt.or.ut(nut).lt.xp_time
     >      .or.ut(nut).lt.endtime_pt.or.ut(1).gt.xp_time)then
             if(iproc.eq.0)
     >       write(*,*)"error: data time range is ",ut(1),ut(nut)
             stop
           endif
           do k=1,nur
             rho_u(k) = ur(k)
           enddo
c          rho_u(1) = 0.D-10
c          rho_u(nur) = 1.D0-1.D-10
           jstart = nut+1
           jend = 0
c
           do j=1,nut
c select times in the range xp_time to endtime_pt
             if(j.lt.nut.and.ut(j+1).lt.xp_time)go to 9
             if(j.gt.1.and.ut(j-1).gt.endtime_pt)go to 9
               if(jstart.gt.j)jstart=j
               if(jend.lt.j)jend=j
c fix the boundaries to 0,1
c              if(ur(1).gt.0.025)then
c                par(1,j) = par(1,j) - 
c     >        (ur(1)-rho_u(1))*(par(2,j)-par(1,j))/(ur(2)-ur(1))
                par(nur,j) = par(nur,j) + 
     >             (ur(nur)-rho_u(nur))*(par(nur,j)-par(nur-1,j))/
     >             (ur(nur)-ur(nur-1))
c              endif
c interpolate onto ptor data grid
              call INTER_CSPL(nur,rho_u,par(1,j),nr,rho_d,pxp)
c save timeseries
              do k=1,nr
               u2d_t(k,j,i) = pxp(k)
              enddo
 9            continue
           enddo
c
           do j=jstart,jend
             jj=j
             if(j.eq.1)jj=j+1    
             if(ut(jj).ge.xp_time.and.ut(jj-1).lt.xp_time)then
c interpolate to the time of interest
               do k=1,nr
                 u2d(k,i) = u2d_t(k,jj-1,i)+
     >           (xp_time-ut(jj-1))*(u2d_t(k,jj,i)-u2d_t(k,jj-1,i))/
     >           (ut(jj)-ut(jj-1))
c                u2d(k,i) = u2d_t(k,jj,i)
               enddo
             endif
           enddo
         endif
        endif
        if(ierr.eq.1) then
           do k=1,mxgrid+1
              u2d(k,i)=0.D0
           enddo
        else if(ierr.eq.2) then
           if(iproc.eq.0) write(*,*) 'array dimensions 
     &                     incompatible, i = ',i
           stop
        endif
      enddo
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Smooth data
c
      if(ismooth_data.ne.0)call average7_1d(u2d(1,i),nr)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Read expdata files
c
c     loop over the 2d files that contain experimental data
c     do not interpolate to the 51 pts grid 
c
      do i=1,n3d
c
c     compose the file name:
         ig=linblnk(tok,10)
         is=linblnk(shot,10)
         i2d=linblnk(fields3d(i),7)
         ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields3d(i)(1:i2d)
c
c Read the 3D ufile:
c
         ierr=0
         nx3(i)=0
         call u3read(cudir,ufile,88,u3dx(1,i),nx3(i),
     >               xp_time,u3dy(1,i),ierr)
c
c...Reset nx3(i) to zero if an error
c
         if(ierr.eq.1) then
            if (iproc.eq.0) write(*,57) ierr,ufile
            nx3(i)=0
            do j=1,mxgrid+1
                  u3dx(j,i)=0.D0
                  u3dy(j,i)=0.D0
            enddo
         else if(ierr.eq.2) then
            if(iproc.eq.0) write(*,*) 'array dimensions 
     &                      incompatible, i = ',i
            stop
         endif
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Loop over the 1d files next:
c
      do i=1,n1d
c
c     compose the file name:
         ig=linblnk(tok,10)
         is=linblnk(shot,10)
         i1d=linblnk(fields1d(i),7)
         ufile=tok(1:ig)//'1d'//shot(1:is)//'.'//fields1d(i)(1:i1d)
c
c Read the 1D ufile:
c
         ierr=0
         if(time_series.eq.0)then
c          call u1read(cudir,ufile,88,xp_time,u1d(i),ierr)
           call u1tread(cudir,ufile,88,nut,ut,u1d_t(1,i),nsmooth,ierr)
           if(ierr.eq.0)then
            do k=1,nut
             j=k
             if(xp_time.le.ut(1)) then
               xp_time=ut(1)+1.D-6
               write(*,35) ut(1)
             endif
c            if(i.eq.1) write(*,*) 'ut,xp_time = ',ut(j),xp_time
             if(k.eq.1)j=k+1
             if(ut(j).ge.xp_time.and.ut(j-1).lt.xp_time)then
c interpolate to the time of interest
                u1d(i) = u1d_t(j-1,i) + (u1d_t(j,i)-u1d_t(j-1,i))*
     >            (xp_time-ut(j-1))/(ut(j)-ut(j-1))
             endif
            enddo 
          endif
         else
           if(i.eq.1 .and. iproc.eq.0) write(6,65)
           do j=1,ntmax
             u1d_t(j,i) = 0.D0
           enddo
           call u1tread(cudir,ufile,88,nut,ut,u1d_t(1,i),nsmooth,ierr)
c
           if(ierr.eq.0)then
            do k=1,nut
             j=k
             if(k.eq.1)j=k+1
             if(ut(j).ge.xp_time.and.ut(j-1).lt.xp_time)then
c interpolate to the time of interest
c                 u1d(i) = u1d_t(j,i)
                u1d(i) = u1d_t(j-1,i) + (u1d_t(j,i)-u1d_t(j-1,i))*
     >            (xp_time-ut(j-1))/(ut(j)-ut(j-1))
             endif
            enddo 
          endif 
         endif
         if(ierr.eq.1) then
	    continue
            if (iproc.eq.0) write(6,57) ierr,ufile
            u1d(i)=0.D0
         else if(ierr.eq.2) then
            if(iproc.eq.0) write(*,*) 'array dimensions 
     &                      incompatible, i = ',i
            stop
         endif
      enddo
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Finally, read the 0D data:
c
         ig=linblnk(tok,10)
         is=linblnk(shot,10)
c
         ufile=tok(1:ig)//'_'//shot(1:is)//'_0d.dat'

         id=linblnk(cudir,50)
c        ifchar=ig+is+2
         ifchar=ig+is+8
         if ( id .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
            ichar  = ifchar            
            udfile = ufile
         else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cudir(id:id) .ne. '/'
     &    .and. cudir(id:id) .ne. ']' ) then
          ichar = id + 1 + ifchar
          udfile = cudir(1:id) // '/' // ufile(1:ifchar)
        else
          ichar  = id + ifchar
          udfile = cudir(1:id) // ufile(1:ifchar)
        endif
c
      endif
c
c..protect against names too long
c
      if ( ichar .gt. 130 ) then
       if (iproc .eq. 0 ) then
         write (*,*)
         write (*,*) 'directory and file names too long'
         write (*,*) ichar
     &        ,' = id + ifchar .gt. 130 in sbrtn read'
         write (*,*) ifchar,' = ifchar'
         write (*,*) id,' = idchar'
         write (*,*) cudir,' = dir'
         write (*,*) ufile,' = ufile'
       endif
      endif
c
c...  Read in 0d ufile data
c     u0read reads each element of a 0d file into a string
c     idatzero=0 old ITER PDB 0d file format
c     idatzero=1 for new 0D file formats
c     idatzero=3 for new 0D file format without DWDIA,DWMHD
c     Note that phase = OHM, L, LHLHL, H, HSELM,
c                       HGELM, HGELMH, VH, PEP
      ifail=0
      call u0read(idatzero,udfile,u0names,u0values,ifail)
c     open(unit=51,file=ufile(1:ichar), err=90)
      if (ifail.ne.0 .and. iproc .eq. 0) then
         write(*,*) 'Failure opening 0D ufile'
      else
         if(iproc.eq.0) then
           write(*,*) ' 0D file read successful, ifail=',ifail
         endif
      endif
      if(iproc.eq.0) write(*,*) ' idatzero = ',idatzero
c
c      do j=1,165
c        write(*,*) j, u0values(j)
c      enddo
c
      p_glob_exp=0.D0
      phase=' '
      if (ifail.eq.0) then
        if(idatzero.eq.0) then
          read(u0values(14),1010) apgasa
          read(u0values(15),1010) apgasz
          read(u0values(16),1010) abgasa
          read(u0values(17),1010) abgasz
          read(u0values(25),1010) apimpa
          read(u0values(26),1010) apimpz
          read(u0values(61),*) bt_0d
          read(u0values(62),1010) ip_0d
          read(u0values(64),1010) q95_0d
c Estimated thermal power loss (not corrected for cx and orbit losses)
          read(u0values(133),1010,err=1011) p_glob_exp
          read(u0values(135),1010,err=1011) gtauth_exp
          read(u0values(134),1010,err=1011) gtautot_exp
          p_glob_exp = p_glob_exp*1.D-6  ! [MW]
          u0phase = u0values(7)
c        write(*,*) 'q95_0d = ',q95_0d
        endif
        if(idatzero.eq.1) then
          read(u0values(94),*) apgasa
          read(u0values(95),*) apgasz
          read(u0values(99),1010) apimpa
          read(u0values(100),1010) apimpz
          read(u0values(8),*) abgasa
          read(u0values(10),*) abgasz
          u0phase = u0values(96)
          read(u0values(115),*) rgeo_0d
          read(u0values(1),*) rmin_0d
          read(u0values(15),*) bt_0d
          read(u0values(51),*) ip_0d
          read(u0values(56),*) kap_0d
          read(u0values(19),*) delt_0d
          read(u0values(64),*) nebar_0d
          read(u0values(106),*) pnbi_0d
          read(u0values(98),*) pich_0d
          read(u0values(91),*) pech_0d
          read(u0values(139),*) ti0_0d
          read(u0values(134),*) te0_0d
          read(u0values(111),*) q0_0d
          read(u0values(110),*) q95_0d
          read(u0values(131),*) gtauth_exp
          read(u0values(133),*) gtautot_exp
          read(u0values(159),*) wth_0d
          read(u0values(160),*) wtot_0d
        endif
        if(idatzero.eq.2) then
          read(u0values(37),*) apgasa
          read(u0values(38),*) apgasz
          read(u0values(73),1010) apimpa
          read(u0values(74),1010) apimpz
          read(u0values(39),*) abgasa
          read(u0values(40),*) abgasz
          u0phase = u0values(8)
          read(u0values(25),*) rgeo_0d
          read(u0values(27),*) rmin_0d
          read(u0values(10),*) bt_0d
          read(u0values(11),*) ip_0d
          read(u0values(32),*) kap_0d
          read(u0values(35),*) delt_0d
          read(u0values(64),*) nebar_0d
          read(u0values(43),*) pnbi_0d
          read(u0values(84),*) pich_0d
          read(u0values(91),*) pech_0d
          read(u0values(56),*) ti0_0d
          read(u0values(53),*) te0_0d
          read(u0values(59),*) q0_0d
          read(u0values(60),*) q95_0d
          read(u0values(20),*) gtauth_exp
          read(u0values(21),*) gtautot_exp
          read(u0values(14),*) wth_0d
          read(u0values(15),*) wtot_0d
        endif
        if(idatzero.eq.3) then
          read(u0values(92),*) apgasa
          read(u0values(93),*) apgasz
          read(u0values(97),1010) apimpa
          read(u0values(98),1010) apimpz
          read(u0values(8),*) abgasa
          read(u0values(10),*) abgasz
          u0phase = u0values(94)
          read(u0values(113),*) rgeo_0d
          read(u0values(1),*) rmin_0d
          read(u0values(15),*) bt_0d
          read(u0values(49),*) ip_0d
          read(u0values(54),*) kap_0d
          read(u0values(19),*) delt_0d
          read(u0values(62),*) nebar_0d
          read(u0values(104),*) pnbi_0d
          read(u0values(96),*) pich_0d
          read(u0values(89),*) pech_0d
          read(u0values(137),*) ti0_0d
          read(u0values(132),*) te0_0d
          read(u0values(109),*) q0_0d
          read(u0values(108),*) q95_0d
          read(u0values(129),*) gtauth_exp
          read(u0values(131),*) gtautot_exp
          read(u0values(157),*) wth_0d
          read(u0values(158),*) wtot_0d
        endif
        ip_0d=ip_0d/1.D6
        nebar_0d=nebar_0d/1.D19
        pnbi_0d=pnbi_0d/1.D6
        pich_0d=pich_0d/1.D6
        pech_0d=pech_0d/1.D6
        ti0_0d=ti0_0d/1.D3
        te0_0d=te0_0d/1.D3
        wth_0d=wth_0d/1.D6
        wtot_0d=wtot_0d/1.D6
c
        if (u0phase(1:1).eq.'o'.or.u0phase(1:1).eq.'O') phase='O'
        if (u0phase(1:1).eq.'l'.or.u0phase(1:1).eq.'L') phase='L'
        if (u0phase(1:1).eq.'h'.or.u0phase(1:1).eq.'H') phase='Hnoelm'
        if (u0phase(2:2).eq.'s'.or.u0phase(2:2).eq.'S') phase='Helm'
        if (u0phase(2:2).eq.'g'.or.u0phase(2:2).eq.'G') phase='Helm'
      endif
c
c     Assume the working and beam ions are the same and the impurity
c     is set by pimpa,pimpz. In xpin, the defaults are the following:
c       amassgas_exp=2.0     amassimp_exp=12.0
c       zgas_exp=1.0         zimp_exp=6.0         
c       pimpa=12.0           pimpz=6.0
c     All five variables may also be set in the input file 
c     if the 0D file information is not present or desired.
c
      if (apgasa.lt.1.D-3) apgasa=amassgas_exp
      if (apgasz.lt.1.D-3) apgasz=zgas_exp
      if (abgasa.lt.1.D-2) abgasa=amassgas_exp
      if (abgasz.lt.1.D-2) abgasz=1.D0
      if (apimpa.lt.1.D-2) apimpa=amassimp_exp
      if (apimpz.lt.1.D-2) apimpz=zimp_exp
c
      if (iproc.eq.0) then
        write(*,11)
        if(idatzero.ge.1)then
          write(*,13) rgeo_0d,rmin_0d,bt_0d,ip_0d,nebar_0d,pnbi_0d,
     &                pich_0D,pech_0d,ti0_0d,te0_0d,q0_0d,q95_0d,
     &                kap_0d,delt_0d,wth_0d,wtot_0d
          write(*,12) apgasa,apgasz,abgasa,abgasz,apimpa,apimpz
        else
          write(*,12) apgasa,apgasz,abgasa,abgasz,apimpa,apimpz
        endif
      endif
c
c...PLTH second best option:
c...  Read total input power corrected for cx and orbit losses
c...  Note that here p_glob_exp can still be .le.0
c...  In that case it is assigned a value (according
c...  to 1D data) just before ITER89P is calculated.
c
      if (p_glob_exp.le.0.0 .and. idatzero.eq.0) then
        if (ifail.eq.0) then
          read(u0values(76),1010,err=1012) p_glob_exp
          p_glob_exp = p_glob_exp*1.D-6                  ![MW]
        endif
      endif
c
 1010 format(1pe10.3)
 1011 continue
 1012 continue
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c nj : the size of the vectors printed in this file
          nj_d=nr
c
c nion : the number of ion species          
          nion_d=2
c
c nprim : the number of primary ion species
          nprim_d=1
c
c nimp : the number of impurity ion species          
          nimp_d=1
c
c nneu : the number of neutral ion species
          nneu_d=1
c
c ibion : index of beam species          
          ibion_d=1
c
c namep : name(s) of primary ion species          
c         (namep_d(i),i=1,nprim_d)
c
c namei : name(s) of impurity ion species          
c         (namei_d(i),i=1,nimp_d)
c
c namen : name(s) of neutral ion species         
c         (namen_d(i),i=1,nneu_d)
c              
c time :  time at which data is printed          
          time_d=xp_time
c
c Rgeom : major radius of geometric          
         rgeom_d=u1d(9)                  !meters
c 
c Rmag :  major radius of mag axis, meters        
         rmag_d=u2d(1,15)                    !meters
         if (rmag_d.eq.0.) rmag_d=u1d(9)
c
c R0 : major radius of vaccuum btor ref            
         rmajor_d=u1d(9)                !meters
         amin_d=u1d(1)                  !meters
c
c JAK 960325 warning output added
         if (rmajor_d.eq.0.0)
     >     write(6,*) 'Warning rep_iter: rmajor zero'
c
c kappa : plasma elongation         
         if(u1d(6).gt.0.) then
           kappa_d=u1d(6)
         else
           u1d(6)=1.D0
           kappa_d=1.D0
         endif
c delta : plasma triangularity         
         deltao_d=u1d(3)
c  
c pindent : plasma indentation         
         pindento_d=u1d(4)
c
c volo : plasma volume,meters**3         
         volo_d=u2d(nr,21)                    !meters**3
c
c Btor : vaccuum toroidal field at rmajor, tesla         
         btor_d=u1d(2)                    !tesla
         if(btscale.gt.1.) btor_d=btor_d*btscale
c
c total, ohmic, bootstrap, beam, and rf currents, amps       
c        tocur_d
c        totohm_d
c        totboot_d
c        totbeam_d
c        totrf_d
c
         tocur_d=u1d(5)
         curtot=tocur_d/1.D6
c
c betap : poloidal beta         
c         betap_d
c
c beta : toroidal beta         
c         beta_d
c
c ali : plasma inductance         
         ali_d=u1d(7)
c
c te0 : central electron temperature         
         te0_d=u2d(1,23)/1.D3                        !kev
c
c ti0 : central ion temperature         
         ti0_d=u2d(1,22)/1.D3                        !kev
c
c pohm: total ohmic input power                      !Wa  (JAK)
         pohm_d=u1d(12)
c
c vsurf: loop voltage at the plasma boundary         !V   (JAK)
         vsurf_d=u1d(13)
c
c---psi on rho grid,volt*sec/rad 
c              
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(psir_d(j),j=1,nj_d)      !volt*se/rad        
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c   
c---rho grid, meters
c
c get rho_a from edge kappa
c r_d = rho*a*sqrt(kappa)
c
          do j=1,nj_d
             r_d(j)=rho_d(j)*u1d(1)*sqrt(u1d(6))    !meters
             r_d(j)=r_d(j)*rho_norm_loc
          enddo
          if (r_d(nj_d).eq.0.0 .and. iproc.eq.0) then
            write(6,*) 'Warning rep_iter: r_d zero'
          endif
c        stop
c
c---flux surface average absolute grad rho and
c   flux surface average ( grad rho)**2
c
         do j=1,nj_d
           grho1npsi_d(j)=u2d(j,5)*r_d(nj_d)
           grho2npsi_d(j)=u2d(j,6)*r_d(nj_d)**2
         enddo
c
c---fcap, (ie f(psilim)/f(psi) )
c   place holder for fcap until can get from U-file  GMS 1-12-2007
         do j=1,nj_d
           fcap_d(j) = 1.D0
         enddo
c
c---gcap, (ie <(grad rho)**2*(R0/R)**2> )
c   place holder for gcap until can get from U-file  GMS 1-12-2007
         do j=1,nj_d
           gcap_d(j) = grho2npsi_d(j)/(1.D0-(r_d(j)/rmajor_d)**2)**1.5
         enddo
c       
c---hcap, (surface area/(4*pi*pi*R0*gradrho*rho))
c         (vprime/(4*pi*pi*R0*gradrho*rho*<abs(gradrho)>))
c         with surface area = vprime*<abs(gradrho)> and
c         vprime = dV/drho
c
           rfix=1.D0
           if(istk.eq.88) rfix=4.46
           sfix=1.D0
           if(istk.eq.88) sfix=22.825
           if(istk.eq.88) write(6,*) 'sfix=22.825'
c
          do j=2,nj_d
            hcap_d(j)=sfix*u2d(j,19)/(u2d(j,5)*
     >      r_d(nj_d)/rfix)/(4.D0*pi_m*pi_m*rmajor_d*r_d(j))
c           write(*,100) j, r_d(j), u2d(j,19), hcap_d(j)
          enddo
          hcap_d(1)=hcap_d(2)
c
c---volume of each flux surface,  meters**3
c
        do j=1,nj_d
          psivolp_d(j)=u2d(j,21)
        enddo
c                         
c--- flux surface area,  meters**2 4.*pi*pi*R0*hcap*rho*<abs(grad rho)>
c
c        u2d(2,19)=0.D0
        if ( u2d(2,19).eq.0 ) then
          if(iproc.eq.0)
     >     write(6,'(a33)') 'Computing flux surface area, hcap'
          do j=2,nj_d-1
           sfareanpsi_d(j)=
     >       (psivolp_d(j+1)-psivolp_d(j-1))/2.D0 /
     >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) * grho1npsi_d(j)
           hcap_d(j)=
     >       (psivolp_d(j+1)-psivolp_d(j-1)) /2.D0 /
     >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) /
     >       (4.D0*pi_m*pi_m*rmajor_d*r_d(j))
          enddo
          sfareanpsi_d(1)=1.D-6
          sfareanpsi_d(nj_d)=
     >       (psivolp_d(nj_d)-psivolp_d(nj_d-1)) /
     >       (rho_d(nj_d)-rho_d(nj_d-1)) / r_d(nj_d) * 
     >       grho1npsi_d(nj_d)
          hcap_d(5)=hcap_d(6)
          hcap_d(4)=hcap_d(5)
          hcap_d(3)=hcap_d(4)
          hcap_d(2)=hcap_d(3)
          hcap_d(1)=hcap_d(2)
        else
          do j=1,nj_d
            sfareanpsi_d(j)=sfix*u2d(j,19)
          enddo
        endif
c        write(*,*) 'sfareanpsi = ',sfareanpsi_d(nj_d)
c        write(*,*) 'area = ',3.14*amin_d**2*kappa_d
c        write(*,*) 'kappa = ',kappa_d
c
c       do j=2,nj_d-1
c         write(*,54) j, r_d(j), sfareanpsi_d(j), 
c    >       (psivolp_d(j+1)-psivolp_d(j-1)) /2.D0 /
c    >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) * grho1npsi_d(j),
c    >        hcap_d(j),
c    >       (psivolp_d(j+1)-psivolp_d(j-1)) /2.D0 /
c    >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) /
c    >       (4.D0*pi_m*pi_m*rmajor_d*r_d(j))
c       enddo
c
c---read original experimental data for ne, te and ti ---JAK
c
c ne and errorbar
c
          nnexp_d = nx3(1)
          do j=1,nnexp_d
           xnexp_d(j)=u3dx(j,1)
           ynexp_d(j)=u3dy(j,1)*1.D-19
           xnexpeb_d(j)=u3dx(j,2)
           ynexpeb_d(j)=u3dy(j,2)*1.D-19
          enddo
c
c te and errorbar
c
          ntexp_d = nx3(3)
          do j=1,ntexp_d
           xtexp_d(j)=u3dx(j,3)
           ytexp_d(j)=u3dy(j,3)*1.D-3        !keV
           xtexpeb_d(j)=u3dx(j,4)
           ytexpeb_d(j)=u3dy(j,4)*1.D-3      !keV
          enddo
c
c ti and errorbar
c
          ntixp_d = nx3(5)
          do j=1,ntixp_d
           xtixp_d(j)=u3dx(j,5)
           ytixp_d(j)=u3dy(j,5)*1.D-3        !keV
           xtixpeb_d(j)=u3dx(j,6)
           ytixpeb_d(j)=u3dy(j,6)*1.D-3      !keV
          enddo
c
c---electron temperature, kev
c
          do j=1,nj_d
           te_d(j)=u2d(j,23)/1.D3            !kev
          enddo
c
c---ion temperature, kev
c
          do j=1,nj_d
           ti_d(j)=u2d(j,22)/1.D3          !kev
          enddo
  
       if( iexp_exch.eq.2) then
        do j=1,nj_d
         ti_d(j)=te_d(j)-.001
        enddo
       endif    
c
c---q (ie safety factor) profile
c
          do j=1,nj_d
           q_d(j)=u2d(j,27)
c           write(*,*) j, r_d(j), q_d(j), 'q_d'
          enddo
          q01_exp=q_d(1)
          qa1_exp=q_d(nj_d)
c
c---electron density,#/m**3
c
          do j=1,nj_d
           ene_d(j)=u2d(j,7)           !#/meter**3
c           write(*,54) j, rho_d(j), ene_d(j)
          enddo                         
c      
c---primary ion density,#/m**3,species (NM1,NM2,NM3)
c---impurity ion density,#/m**3,species (NIMP)
c   assuming abgasz=apgasz
             do j=1,nj_d
               en_d(j,1)=u2d(j,10) + u2d(j,42) + u2d(j,43) !#/meter**3
               en_nm1_d(j)=u2d(j,10)
               en_nm2_d(j)=u2d(j,42)
               en_nm3_d(j)=u2d(j,43)
c              if (imodel.lt.0) then
c                en_d(j,1)=u2d(j,10) + u2d(j,42) + 
c    >                     u2d(j,43) + u2d(j,9) ! chi=1 test
c              endif
               en_d(j,2)=u2d(j,9)    !#/meter**3
c               write(*,100) j, rho_d(j),u2d(j,20)
               zeff_t=u2d(j,20)
               if (zeff_t.eq.0.) zeff_t=u1d(11)
               if (zeff_t.eq.0.) then
                 zeff_t=1.
                 if(iproc.eq.0) write(6,*) 'zeff_t=1.'
               endif
c
               if(iproc.eq.0 .and. apimpz.lt.2.D0)
     >           write(6,*) '*** WARNING: apimpz looks fishy ***'
               if(en_d(j,1).eq.0.) 
     >          en_d(j,1)=(apimpz-zeff_t)/(apimpz-apgasz)/apgasz
     >                    *u2d(j,7)-u2d(j,8)
               if(en_d(j,2).eq.0.) 
     >          en_d(j,2)=(zeff_t-apgasz)/(apimpz-apgasz)/apimpz
     >		              *u2d(j,7)
               if(en_d(j,1).lt.0..or.en_d(j,2).lt.0.) then
                 if(iproc.eq.0) write(6,*) 'negative density fixed'
                 en_d(j,1)=dabs(en_d(j,1))
                 en_d(j,2)=dabs(en_d(j,2))
               endif
c               write(*,100) j, rho_d(j), ene_d(j), en_d(j,1), en_d(j,2)
             enddo
c
c---sion   : source due to ionization,        #/(m**3*sec),species
c---srecom : source due to recombination,     #/(m**3*sec),species
c---scx    : source due to cx thermal neut.,  #/(m**3*sec),species
c---sbcx   : sink due to cx with beam neut.   #/(m**3*sec),species
c---s      : total source rate,               #/(m**3*sec),species
c---dudt   : s dot,                           #/(m**3*sec),species
c
       do jj=1,nprim_d       
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(sion_d(j,jj),j=1,nj_d)    !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c              read(niterdb,10)(srecom_d(j,jj),j=1,nj_d)  !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(scx_d(j,jj),j=1,nj_d)     !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(sbcx_d(j,jj),j=1,nj_d)    !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(s_d(j,jj),j=1,nj_d)       !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(dudtsv_d(j,jj),j=1,nj_d)  !#/meter**3/sec
c
            do j=1,nj_d
              sion_d(j,jj)=0.
              srecom_d(j,jj)=0.
              scx_d(j,jj)=0.
              sbcx_d(j,jj)=0.
              s_d(j,jj)=0.
              dudtsv_d(j,jj)=0.
             enddo
          enddo
c
c---fast ion density, #/m**3, species
c
          do j=1,nj_d
             enbeam_d(j)=u2d(j,8)         !#/meter**3
          enddo
c      
c---neutral density,  #/m**3,species
c
       do jn=1,nneu_d
           do j=1,nj_d
              enn_d(j,jn)=0.     !#/meter**3
          enddo
       enddo
c
c---neutral density from  wall source, #/m**3, species
c
       do jn=1,nneu_d        
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(ennw_d(j,jn),j=1,nj_d)     !#/meter**3     
         do j=1,nj_d
          ennw_d(j,jn)=0.
         enddo
       enddo
c
c---neutral density from volume source, #/m**3,species (recomb and beam cx)
c
       do jn=1,nneu_d
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(ennv_d(j,jn),j=1,nj_d)     !#/meter**3
        do j=1,nj_d
         ennv_d(j,jn)=0.
        enddo
       enddo
c
c---volume source of neutrals, #/(m**3*sec),species
c
c       do jn=1,nneu_d   
c               read(niterdb,'(a)')stflg
cx               read(niterdb,10)(volsn_d(j,jn),j=1,nj_d)!#/(m**3*sec)
c               read(niterdb,10)(blank_d(j),j=1,nj_d) 
c       enddo
c
c---sbion : beam electron source,  #/(m**3*sec
c        
c               read(niterdb,'(a)')stflg
cx               read(niterdb,10)(sbion_d(j),j=1,nj_d)    !#/(m**3*sec)
c               read(niterdb,10)(blank_d(j),j=1,nj_d)          
c
c---sbeam : beam thermal ion source,  #/(m**3*sec)
c
c               read(niterdb,'(a)')stflg   
c               read(niterdb,10)(sbeam_d(j),j=1,nj_d)      !#/(m**3*sec)
c
c...Note that SNBII as read contains all sources and sinks due to
c...beam thermal ion source and charge exchange processes; 
c...so: field SNBII = sbeam + scx + sbcx + srecom    (note sbcx<0 sink)
c...Read the source in sbeam_d and let all other source terms be 0.          
c
             do j=1,nj_d
              sbeam_d(j)=u2d(j,18)                       !#/(m**3*sec)
             enddo
c
c...thermal ion particle source due to ionisation; field SWALL
c
             do j=1,nj_d
              sion_d(j,1)=u2d(j,30)                       !#/(m**3*sec)
             enddo
             do j=1,nj_d
               dudtsv_d(j,1)=u2d(j,41)                    ! dn/dt
             enddo
c
c---total current density, amps/m**2
c   (avoids negative currenty density)
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(curden_d(j),j=1,nj_d)        !amps/meter**2
             do j=1,nj_d
              curden_d(j)=DABS(u2d(j,4))
             enddo
c    
c---ohmic current density, amps/m**2
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curohm_d(j),j=1,nj_d)        !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d) 
c          
c---bootstrap current density, amps/m**2
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curboot_d(j),j=1,nj_d)       !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c      
c--- beam driven current density, amps/m**2
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curdbeam_d(j),j=1,nj_d)      !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---rf current density, amps/m**2
c    
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(currf_d(j),j=1,nj_d)         !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c          
c---rho*bp0*fcap*gcap*hcap, tesla*meters
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(rbp_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)         !tesla*meters      
c
c---zeff profile:
c
          do j=1,nj_d
           zeff_d(j)=u2d(j,20)
           if (zeff_d(j).eq.0.) zeff_d(j)=u1d(11)
c          write(*,54) j, rho_d(j), zeff_d(j)
          enddo
c
c---angular rotation speed profile, rad/sec
c
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(angrot_d(j),j=1,nj_d)       !rad/sec
          do j=1,nj_d
            angrot_d(j)=u2d(j,29)
          enddo
c
c---torque density, N/M**2
c Note: TRANSP ufiles have N-m/cm**3
c
          if(itorque.ne.0)then
            if(tok.eq.'jet') conv_torq=1.D0
            do j=1,nj_d
              torque_d(j)=u2d(j,46)*conv_torq
            enddo
          endif
c
c---electron thermal diffusivity, meters*2*/sec on half grid
c
          do j=1,nj_d
           chieinv_d(j)=u2d(j,1)       !meters**2/sec
          enddo
c                 
c---ion thermal diffusivity,meters*2*/sec  on half grid
c     
          do j=1,nj_d
           chiinv_d(j)=u2d(j,2)       !meters**2/sec
          enddo
c         
c---ion neocl.  thermal conductivity, 1/(m*sec) on half grid 
c        
          do j=1,nj_d
           xkineo_d(j)=0.D0      !meters**2/sec
          enddo
c
c---wdot,electrons, watts/m**3 d(electron energy)/dt profile (DWER):
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(dpedtc_d(j),j=1,nj_d)        !watts/meter**3      
          do j=1,nj_d
           dpedtc_d(j)=u2d(j,31)
          enddo
c    
c---wdot,ions,watts/m**3 d(ion energy)/dt profile (DWIR):
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(dpidtc_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           dpidtc_d(j)=u2d(j,32)
          enddo
c
c---electron conduction, watts/m**3
cc      
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qconde_d(j),j=1,nj_d)        !watts/meter**3 
          do j=1,nj_d
           qconde_d(j)=0.D0
          enddo
c   
c---ion conduction, watts/m**3
c
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qcondi_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qcondi_d(j)=0.D0
          enddo
c
c---electron convection,watts/m**3 
c     
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qconve_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qconve_d(j)=0.D0
          enddo
c
c---ion convection, watts/m**3 
c     
c         read(niterdb,'(a)')stflg
c         read(niterdb,10)(qconvi_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qconvi_d(j)=0.D0
          enddo
c
c---power to elec.from beam, watts/m**3 
c     
           do j=1,nj_d
             qbeame_d(j)=u2d(j,11)        !watts/meter**3
           enddo
c          
c---qdelt : electron-ion equilibration, watts/m**3
c
         do j=1,nj_d
          qdelt_d(j)=u2d(j,28)
         enddo
c
c NRL formulary :
c
       if(qdelt_d(nj_d).eq.0.or.iexp_exch.ge.1) then       
         do j=1,nj_d
           zbrac=en_d(j,1)/ene_d(j)+apgasa/
     >         apimpa*en_d(j,2)/ene_d(j)*apimpz**2
           alamda=24.D0-dlog((dabs(ene_d(j))*1.D-6)**zhalf/
     >         (dabs(te_d(j))*1.D3))
           tau_e=3.44D5*(dabs(te_d(j))*1.D3)**1.5D0/
     >         (ene_d(j)*1.D-6)/alamda
           qdelt_d(j)=-10**6*1.6022D-19*zbrac*3.D0/(1836.D0*apgasa)/
     >         tau_e*ene_d(j)*1.D-6*(te_d(j)-ti_d(j))*1.D3
         enddo                    
       endif     !watts/meter**3
       if(iexp_exch.eq.-1) then
        do j=1,nj_d
         qdelt_d(j)=0.D0
        enddo
       endif
c
c---power to ions from beam, watts/m**3
c   
           do j=1,nj_d
             qbeami_d(j)=u2d(j,12)        !watts/meter**3
           enddo
c    
c---qrfe, rf electron heating, watts/m**3 (QECHE + QICRHE)
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qrfe_d(j),j=1,nj_d)          !watts/meter**3 
           do j=1,nj_d
            qrfe_d(j)=u2d(j,33)*echconv + u2d(j,34)
           enddo
c    
c---qrfi, rf ion heating, watts/m**3 (QECHI + QICRHI)
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qrfi_d(j),j=1,nj_d)          !watts/meter**3 
           do j=1,nj_d
            qrfi_d(j)=u2d(j,35) + u2d(j,36)
           enddo
c
c... calculate integrated RF heating
c
       prfesum=0.D0
       prfisum=0.D0
       do j=2,nj_d
         drm=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*
     >      hcap_d(j-1))/2.D0*drm
         prfesum=prfesum+
     >     zhalf*(dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*drm
         prfisum=prfisum+
     >     zhalf*(dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*drm
       enddo
c      if(iproc.eq.0) then
c        write(6,'(a27,2F10.3)') 'Total integrated Prfe,Prfi [MW]:',
c    >        prfesum*1.D-6,prfisum*1.D-6
c      endif
c
c
c---qlhe, lower hybrid heating to electrons, watts/m**3
c   
        do j=1,nj_d
          qlhe_d(j)=u2d(j,45)
        enddo
c
c---qione, recombination and impact ionization, watts/m**3
c   use QWALLE = ''thermal electron heat loss due
c   to the ionization of wall neutrals W/m3''
c
        do j=1,nj_d
          qione_d(j)=u2d(j,37)
        enddo
c      
c---qioni, recombination and impact ionization, watts/m**3
c   use QWALLI = ''main thermal ion heat loss due to
c   ionization and charge exchange with wall neutrals in W/m3''
c   SO: QWALLI=qioni_d + qcx_d
c   
        do j=1,nj_d
          qioni_d(j)=u2d(j,38)
        enddo
c      
c---qcx,  neutral-ion charge exchange, watts/m**3
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qcx_d(j),j=1,nj_d)           !watts/meter**3
        do j=1,nj_d
          qcx_d(j)=0.D0
        enddo
c     
c---2d MHD equil. electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qe2d_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---2d MHD equil.ion heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qi2d_d(j),j=1,nj_d)         !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---fusion electron heating, watts/m**3
c
        do j=1,nj_d
          qfuse_d(j)=u2d(j,39)
        enddo
c               
c---fusion ion heating, watts/m**3
c
        do j=1,nj_d
          qfusi_d(j)=u2d(j,40)
        enddo    
c         
c---beam fusion electron heating,watts/m**3
c   
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfuse_d(j),j=1,nj_d)       !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c                
c---beam fusion ion heating profile
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfusi_d(j),j=1,nj_d)       !watts/meter**3
c         read(niterdb,10)(blank_d(j),j=1,nj_d)
c                 
c---qmag electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qmag_d(j),j=1,nj_d)         !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---sawtooth electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawe_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---sawtooth ion  heating, watts/m**3
c       
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawi_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---radiated power density, watts/m**3
c       
        do j=1,nj_d
          qrad_d(j)=u2d(j,14)
        enddo
c
c---(electron) ohmic power density,watts/m**3
c    
        do j=1,nj_d
          qohm_d(j)=u2d(j,13)
        enddo
c
c--- Total plasma pressure, Kev/m**3
c    Pa = N/m**2, keV/m**3=Pa/1.602e-16
c    PTOTR has units of Pascals
c    Use thermal pressure if total pressure not provided
c
       if(u2d(1,44) .eq. 0.0) then
         if(iproc.eq.0) then
          write(6,'(a31)') 'Error: total pressure not found'
          write(6,'(a32)') 'Computing total thermal pressure'
         endif
         do j=1,nj_d
            ptot_d(j) = ene_d(j)*(te_d(j) + ti_d(j))*1.6022D-16   ! Pascals
            pfast_d(j) = 0.0
c            write(*,53) j, r_d(j)/r_d(nj_d),ptot_d(j)
         enddo
       else
         do j=1,nj_d
            ptot_d(j)=u2d(j,44)    
            pfast_d(j)=ptot_d(j) 
     >       - ene_d(j)*(te_d(j)+ti_d(j))*1.6022D-16
c     >            (ene_d(j) - (zimp_exp-1.D0)*
c     >            en_d(j,nprim_d+1))*ti_d(j))
            if(ipfst.eq.1) pfast_d(j)=0.D0
c           write(*,53) j, r_d(j)/r_d(nj_d),ptot_d(j),pfast_d(j)
         enddo
       endif
c
c... calculate integrated radiation loss
c
       pradsum=0.D0
       do j=2,nj_d
         drm=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*
     >            2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >            2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*
     >      hcap_d(j-1))/2.D0*drm
         pradsum=pradsum+
     >     zhalf*(dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*drm
       enddo
c      if(iproc.eq.0) write(6,'(a27,F10.3)') 
c    >   'Total integrated Prad [MW]:', pradsum*1.D-6
c
c... calculate integrated ohmic heating power
c
       if (pohm_d.eq.0.) then
         pohsum=0.D0
         do j=2,nj_d
           drm=r_d(j)-r_d(j-1)
           dvoldr_p=2.D0*pi_m*r_d(j)*
     >              2.D0*pi_m*rmajor_d*hcap_d(j)
           dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >              2.D0*pi_m*rmajor_d*hcap_d(j-1)
           pohsum=pohsum+
     >      zhalf*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
         enddo
       if(iproc.eq.0) write(6,'(a27,F10.3)') 
     >   'Total integrated Pohm [MW]:', pohsum*1.D-6
       else
       if(iproc.eq.0) write(6,'(a30,F7.3)') 
     >   'Total Pohm [MW] from 1D ufile:',pohm_d*1.D-6
       endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c...Calculate q-profile if not available
c   using qa1_exp and q01_exp values
c   alfj_exp=(qa1_exp/q01_exp-1.) corresponds to ja=0 and shata_exp=2.
c
       qa1_exp=u1d(8)
       q01_exp=q0_exp
       alfj1_exp=(qa1_exp/q01_exp-1.D0)
      if(q_d(nj_d).eq.0.or.iexp_q.eq.-1) then
        if(iproc.eq.0) write(6,*) 'Calculating q(r) ...'
c
      do j=2,nj_d
       xxq=rho_xp(j)**2
       q_d(j)=1.D0/(1.D0/(q01_exp*xxq*(alfj1_exp+1.D0))*
     >               (1.D0-(1.D0+1.D-10-xxq)**(alfj1_exp+1)))
      enddo
c
      if(qmin_exp.ne.0.) then
      do j=2,nj_d
       xxq=rho_xp(j)**2
       q_d(j)=(rho_xp(j)-rho_qm_exp)**2/(1.D0-rho_qm_exp)**2
     >   *(qa_exp-qmin_exp)+qmin_exp
      enddo 
      endif
         q_d(1)=q_d(2)
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c...Calculate J(r) from q-profile if not available
c
           if (curden_d(1).EQ.0) then 
c
             if(iproc.eq.0) 
     &          write(6,*) 'Calculating j(r) from q-profile'
c
c...         Check if total current is given
c
c...         Normalize to total current; if not available then use
c...         d3D formula to relate Ip (in MegAmps) to q_95, etc., (obtained from
c...         Jim Deboo in Feb. 1995., who says it usually agrees with efit
c...         to within 5-10%):
c...         1.e6 convert totohm_d from MA to A,
c
             if (tocur_d.EQ.0) then
               if(iproc.eq.0) write(6,*) 'Total 
     &            current calculated from q(r)'
               tocur_d=(5.D0*r_d(nj_d)**2*btor_d/rmajor_d/q_d(nj_d))
     &         *(1.D0+1.5D0*(r_d(nj_d)/rmajor_d)**2)
     &         *(1.D0+kappa_d**2)/2.D0*1.D6
             endif
c
c...         Assume j(r)=j(0)*(1-rhohat**2)**(alpha)
c...         Note that by approximation (circular geometry)
c...         alpha = qa/qa0-1
c
             alpha = qa1_exp/q01_exp-1.D0
             sum=0.D0
             do j=2,nj_d
               drm=r_d(j)-r_d(j-1)
               dvoldr_p=2.D0*pi_m*r_d(j)*
     >                  2.D0*pi_m*rmajor_d*hcap_d(j)
               dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >                  2.D0*pi_m*rmajor_d*hcap_d(j-1)
               ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+
     >            r_d(j-1)*hcap_d(j-1))/2.D0*drm
               curden_d(j) = (1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
               sum = sum + curden_d(j)*ds
c    >           zhalf*(dvoldr_p*curden_d(j)+dvoldr_m*curden_d(j-1))*drm
c    >           /(2.D0*pi_m*rmajor_d) 
             enddo
             cur0 = dabs(tocur_d)/sum
c
             do j=1,nj_d
               curden_d(j)=cur0*(1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
             enddo
           endif 
c
         if (qohm_d(1).EQ.0) then
c
           if(iproc.eq.0) write(6,*)
     >       'Warning rep_iter: no ohmic power profile; making one up'
c
c...       calculate qohm_d from j(r) and loop voltage
c...       by: qohm = E.j
c...       Also calculate volume-integrated power 
c
           cursum=0.D0
           pohsum=0.D0
           do j=2,nj_d
             drm=r_d(j)-r_d(j-1)
             dvoldr_p=2.D0*pi_m*r_d(j)*
     >                2.D0*pi_m*rmajor_d*hcap_d(j)
             dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >                2.D0*pi_m*rmajor_d*hcap_d(j-1)
             ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+
     >          r_d(j-1)*hcap_d(j-1))/2.D0*drm
             cursum = cursum + curden_d(j)*ds
crew error fixed 7.1.96
crew             qohm_d(j)=ABS(vsurf_d*curden_d(j))/(2.*pi_m*rmajor_d)
             qohm_d(j)=(vsurf_d*curden_d(j))/(2.D0*pi_m*rmajor_d)
             pohsum=pohsum+
     >         zhalf*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
           enddo
           qohm_d(1)=(vsurf_d*curden_d(1))/(2.D0*pi_m*rmajor_d) 
crew....vsurf_d or curden_d may not be reliable...
crew     better renorm by pohm_d
crew thus in the end we use only a shape of current profile
crew and assume dE/dr=0 ie steady state current profile
           if(pohm_d.gt.0.) then
            do j=1,nj_d
              qohm_d(j)=pohm_d/pohsum*qohm_d(j)
            enddo
           endif
crew why calc?           pohm_d = pohsum
         if(iproc.eq.0) write(6,'(a27,F10.3,a10,F7.3)') 
     >     'Total calculated Pohm [MW]:', pohsum*1.D-6,
     >     'Vloop = ',vsurf_d
         if(iproc.eq.0) write(6,'(a27,F10.3)') 
     >     'Total integrated Ip  [MA]:', cursum*1.D-6
c
         end if  !end making up ohmic profile
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c---avg major radius of each flux surface, meters, at elevation of mag. axis

c          read(niterdb,'(a)')stflg
          do j=1,nj_d
             rmajavnpsi_d(j) = u2d(j,15)   !meters
          enddo
c
c---avg minor radius of each flux surface,meters, at elevation of mag. axis
c
          do j=1,nj_d
             rminavnpsi_d(j) = u2d(j,16)   !meters
             if(rminavnpsi_d(1).lt.0) rminavnpsi_d(1)=1.D-6
          enddo  
c--- elongation of each flux surface
c
          do j=1,nj_d
           elongx_d(j)=u2d(j,26)
          enddo
c
c use 1D value for NTCC benchmarking if itest_ntcc.gt.0
c
         if( elongx_d(nj_d).lt.1.e-3) then
          if(itest_ntcc.gt.0) then
            write(*,32)
            do j=1,nj_d
              elongx_d(j)=u1d(6)
            enddo
          else
            do j=1,nj_d
              elongx_d(j)= (r_d(j)/rminavnpsi_d(j))**2
            enddo
          endif
         endif
c
c--- triangularity at each flux surface
          do j=1,nj_d
           deltax_d(j)=u2d(j,24)
          enddo 
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 999  return
c
 89   continue
      if(iproc.eq.0) 
     >   write(*,*) 'Error reading iter namelist'
      close(50)
      goto 999
c     
 90   continue
      if(iproc.eq.0) 
     >   write(*,*) 'Error reading the 0d file'
      stop
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 11   format(/,'0D quantities:',/)
 12   format(/,2x,'A-plasma = ',0p1f10.3,4x,'Z-plasma = ',0p1f10.3,/,
     >         2x,'A-beam   = ',0p1f10.3,4x,'Z-beam   = ',0p1f10.3,/,
     >         2x,'A-imp    = ',0p1f10.3,4x,'Z-imp    = ',0p1f10.3,/)
 13   format(2x,'R0  = ',0p1f10.3,4x,'a   = ',0p1f10.3,/,
     >       2x,'BT  = ',0p1f10.3,4x,'Ip  = ',0p1f10.3,/,
     >       2x,'ne  = ',0p1f10.3,4x,'Pnb = ',0p1f10.3,/,
     >       2x,'Pic = ',0p1f10.3,4x,'Pec = ',0p1f10.3,/,
     >       2x,'Ti0 = ',0p1f10.3,4x,'Te0 = ',0p1f10.3,/,
     >       2x,'q0  = ',0p1f10.3,4x,'q95 = ',0p1f10.3,/,
     >       2x,'kappa = ',0p1f8.3,4x,'delta = ',0p1f8.3,/,
     >       2x,'Wth = ',0p1f10.3,4x,'Wtot= ',0p1f10.3)
 32   format(' Using 1D ufile value for elongation')
 35   format(' Warning: xp_time <= min exp. time, setting xp_time = ',
     >         0pf10.5)
 50   format(i2,2x,'rho',5x,'powibeam',5x,'powiion',6x,'powicx',7x,
     &       'powei',8x,'powiwdot',5x,'powitot')
 53   format(i2,2x,0p1f4.2,0p6e13.5)
 54   format(i2,2x,0p1f5.3,1p6e13.5)
 55   format(i2,2x,'rho',5x,'powebeam',5x,'poweohm',6x,'poweion',6x,
     &       'powerad',6x,'powei',8x,'powewdot',5x,'powetot')
 57   format('   ierr = ',i3,',ufile=',a50)
 58   format(2x,i2,2x,0p1f10.6,1p6e15.7)
 60   format(' Reading time-dependent 2D data ...')
 65   format(' Reading time-dependent 1D data ...')
 91   format(' Smoothing electron density profile ...',0p1f4.2)
 92   format(' Smoothing ion density profile ...',0p1f4.2)
 93   format(' Smoothing Zeff profile ...',0p1f4.2)
 94   format(' Smoothing q-profile ...',0p1f4.2)
 95   format(' Smoothing Qnb profile ...',0p1f4.2)
 96   format(' Smoothing fast ion density profile ...',0p1f4.2)
 97   format(' Smoothing vrot profile ...',0p1f4.2)
 98   format(' Smoothing torque profile ...',0p1f4.2)
101   format(' Smoothing grho profile ...',0p1f4.2)
100   format(i2,2x,0p1f4.2,1p8e13.5)
150   format(2x,i2,2x,0p1f10.6,1p6e15.7)
c
      end
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      integer function linblnk(string,max)
      character string*(*)
      integer i,j,max
      i=0
      do j=1,max
         if ( string(j:j)  .ne. ' ' ) then
            i = i + 1
         else
            go to 10
         endif
      enddo
 10   continue 
      linblnk=i
      return
      end

c@u3read.f   17 Nov 1994   Glenn Bateman, PPPL
c Adapted by J.A.Konings to read experimental data without extrapolation
c to 51 pts grid 7/07/95
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u3read ( cdir, cfile, kufile,
     & ppxu, kkxu, ppyu, pdata, kerr )
c
c   This routine reads 2-D ASCII U-files.
c
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  kkxu      = number of elements in x array                    (output)
c  pxu(jx)   = 1-D independent variable x array from U-file     (output)
c  ppyu       = time of desired profile				(input)
c  parray(jx,jy) = 2-D array from U-file                        (output)
c  
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the U-file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer jmaxm,kmaxm
      parameter (jmaxm=301,kmaxm=2000)
      real*8 parray(jmaxm,kmaxm), pxu(jmaxm), pyu(kmaxm), zu(6)
      real*8 pdata(jmaxm), ppxu(jmaxm), ppyu
      real*8 pxu1, pxu0, parray1, parray0
c
      integer kxdim, kxu, kydim, kyu, kufile, kdchar, kfchar, kerr
c
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, jy, ju, i, ii, kkxu
c
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c
c  idchar  = length of directory name cdir
c  ifchar  = length of U-file name cfile
c  iufile  = unit number for reading UFILE
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
      character cdir*(*), cfile*(*), cdfile*130, cline*72
c
c
c..set defaults
c
      ierr   = kerr
c
      kerr   = 0
c
      iprint = 0
c
      ichmax = 60
c
c kxdim is the maximum number of radial points in the 2D ufiles.
c
      kxdim=jmaxm
c
c kydim is the maximum number of time points in the 2D ufiles.
c
      kydim=kmaxm
c
c..clear arrays
c
      do jx=1,kxdim
        pxu(jx)    = 0.D0
      enddo
c
      do jy=1,kydim
        pyu(jy)    = 0.D0
      enddo
c
      do jy=1,kydim
        do jx=1,kxdim
          parray(jx,jy) = 0.D0
        enddo
      enddo
c
c..determine the length of the directory and file names
c  Note:  These names must have no imbedded blanks.
c
      idchar = 0
      do j=1,ichmax
        if ( cdir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 10
        endif
      enddo
  10  continue
c
      ifchar = 0
      do j=1,ichmax
        if ( cfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 11
        endif
      enddo
  11  continue
c
c..protect against no file name given
c
      if ( ifchar .lt. 1 ) then
        write (*,*)
        write (*,*) 'No file name is given for 2-D U-file'
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u3read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u3read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( idchar .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
        ichar  = ifchar
        cdfile = cfile
c
      else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cdir(idchar:idchar) .ne. '/'
     &    .and. cdir(idchar:idchar) .ne. ']' ) then
c
          ichar = idchar + 1 + ifchar
          cdfile = cdir(1:idchar) // '/' // cfile(1:ifchar)
c
        else
c
          ichar  = idchar + ifchar
c
          cdfile = cdir(1:idchar) // cfile(1:ifchar)
c
        endif
c
      endif
c
      if ( iprint .gt. 0 ) then
        write (*,*) cdfile(1:ichar),' = cdfile'
        write (*,*) ichar,' = idchar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u3read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u3read'
          kerr = 1
          if ( ierr .lt. 0 ) stop
          return
        endif
c
c..open UFILE
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) iufile = 70
c
      open ( iufile, file=cdfile(1:ichar), status='old', err=90 )
c
c
c..read header
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
 100  format (a)
c
c..read number of scalars
c
      read (iufile,102) iscalar
 102  format (i4)
c
c..diagnostic printout
c
      if ( iprint .gt. 8 ) then
        write (*,*)
        write (*,*) cdfile(1:ichar),' = name of dir/file'
        write (*,*) cfile,' = name of UFILE 2'
        write (*,*) iscalar,' = iscalar'
      endif
c
c..read lines of scalars
c
      if ( iscalar .gt. 0 ) then
        do j=1,iscalar*2
          read (iufile,100,err=92) cline
        enddo
      endif
c
c..read variable labels
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94) kxu, kyu
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kxu,' = kxu'
        write (*,*) kyu,' = kyu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kyu .gt. kydim ) then
        write (*,*) 'Too many elements in 2-D U-file arrays'
        write (*,*) kxu, kxdim, kyu, kydim
     &   ,' = kxu, kxdim, kyu, kydim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn u3read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x and y 1-D arrays
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
      enddo
c
      iu = 7
      do jy=1,kyu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kyu + 1 - jy )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pyu(jy) = zu(iu)
      enddo
c
c
c..read 2-D profile
c
 
      iu = 7
      do jy=1,kyu
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kyu*kxu + 1 - (jy-1)*kxu - jx )
            read ( iufile, 111, err=92 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx,jy) = zu(iu)
        enddo
      enddo
  
c
c
c..end of routine
c
      close ( iufile, err=90 )
c
c
c select timeslice in case more measurements are available
c In case of only one timeslice: this one is selected.
c
	kkxu=kxu
		
	i=0
	do ii=1,kyu
	  i=i+1
	  if(pyu(ii).ge.ppyu) goto 95
	enddo
  95    continue
  
         do j=1,kkxu
          ppxu(j)=pxu(j)
          pdata(j)=parray(j,i)
         enddo
        
      return
c
c
c..error conditions
c
  90  continue
c
c
c...JAK more compact output
crew      write (6,'(a22,a20)') 'cannot find 2d UFILE:',cfile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn u3read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) kxu,' = kxu'
      write (*,*) kyu,' = kyu'
      write (*,*)
     & 'error reading kxu or kyu in UFILE in sbrtn u3read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u3read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
c
 111  format (1x,6e13.6)
 114  format (1p6e12.4)
 116  format (1p6e13.6)
c
      end

c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u2tread ( cdir, cfile, kufile,
     & kxu, kyu, pxu, pyu, parray, nsmooth, it, kerr)
c
c   This routine reads 2-D ASCII U-files.
c
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  kxu       = number of elements in x array from U-file        (output)
c  kyu       = number of elements in y array from U-file        (output)
c  pxu(jx)   = 1-D independent variable x array from U-file     (output)
c  pyu(jy)   = 1-D independent variable y array from U-file     (output)
c  parray(jx,jy) = 2-D array from U-file                         
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the U-file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer nrmax, ntmax,it, nsmooth
      parameter(nrmax=51) !max number of radial points in experimental grid
c      parameter(ntmax=801) !max number of time points in experimental grid
      parameter(ntmax=3801) !max number of time points in experimental grid
c     parameter(nrmax=301) !max number of radial points in experimental grid
c     parameter(ntmax=2500) !max number of time points in experimental grid
      real*8 parray(nrmax,ntmax), pxu(nrmax), pyu(ntmax), zu(6)
      integer kxdim, kxu, kydim, kyu, kufile, kdchar, kfchar, kerr
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, jy, ju, i, ii, kkxu
c
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c
c  idchar  = length of directory name cdir
c  ifchar  = length of U-file name cfile
c  iufile  = unit number for reading UFILE
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
      character cdir*(*), cfile*(*), cdfile*130, cline*72
c
c
c..set defaults
c
      ierr   = kerr
c
      kerr   = 0
c
      iprint = 0
c
c      if(it.eq.27) iprint = 10  ! for diagnosti  purposes on a given 2d file
c
      ichmax = 60
c
c kxdim is the maximum number of radial points in the 2D ufiles.
c
      kxdim=nrmax
c
c kydim is the maximum number of time points in the 2D ufiles.
c
      kydim=ntmax
c
c..clear arrays
c
      do jx=1,kxdim
        pxu(jx)    = 0.D0
      enddo
c
      do jy=1,kydim
        pyu(jy)    = 0.D0
      enddo
c
      do jy=1,kydim
        do jx=1,kxdim
          parray(jx,jy) = 0.D0
        enddo
      enddo
c
c..determine the length of the directory and file names
c  Note:  These names must have no imbedded blanks.
c
      idchar = 0
      do j=1,ichmax
        if ( cdir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 10
        endif
      enddo
  10  continue
c
      ifchar = 0
      do j=1,ichmax
        if ( cfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 11
        endif
      enddo
  11  continue
c
c..protect against no file name given
c
      if ( ifchar .lt. 1 ) then
        write (*,*)
        write (*,*) 'No file name is given for 2-D U-file'
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u2read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u2read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( idchar .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
        ichar  = ifchar
        cdfile = cfile
c
      else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cdir(idchar:idchar) .ne. '/'
     &    .and. cdir(idchar:idchar) .ne. ']' ) then
c
          ichar = idchar + 1 + ifchar
          cdfile = cdir(1:idchar) // '/' // cfile(1:ifchar)
c
        else
c
          ichar  = idchar + ifchar
c
          cdfile = cdir(1:idchar) // cfile(1:ifchar)
c
        endif
c
      endif
c
      if ( iprint .gt. 0 ) then
        write (*,*) cdfile(1:ichar),' = cdfile'
        write (*,*) ichar,' = idchar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u2read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u2read'
          kerr = 1
          if ( ierr .lt. 0 ) stop
          return
        endif
c
c..open UFILE
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) iufile = 70
c
      open ( iufile, file=cdfile(1:ichar), status='old', err=90 )
c
c
c..read header
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
 100  format (a)
c
c..read number of scalars
c
      read (iufile,102) iscalar
 102  format (i4)
c
c..diagnostic printout
c
      if ( iprint .gt. 8 ) then
        write (*,*)
        write (*,*) cdfile(1:ichar),' = name of dir/file'
        write (*,*) cfile,' = name of UFILE 2'
        write (*,*) iscalar,' = iscalar'
      endif
c
c..read lines of scalars
c
      if ( iscalar .gt. 0 ) then
        do j=1,iscalar*2
          read (iufile,100,err=92) cline
        enddo
      endif
c
c..read variable labels
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94) kxu, kyu
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kxu,' = kxu'
        write (*,*) kyu,' = kyu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kyu .gt. kydim ) then
        write (*,*) 'Too many elements in 2-D U-file arrays'
        write (*,*) kxu, kxdim, kyu, kydim
     &   ,' = kxu, kxdim, kyu, kydim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn u2read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x and y 1-D arrays
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
        if(iprint.ge.10) write(*,*) jx, pxu(jx),' pxu'
      enddo
c
      iu = 7
      do jy=1,kyu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kyu + 1 - jy )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pyu(jy) = zu(iu)
        if(iprint.ge.10) write(*,*) jy, pyu(jy),' pyu'
      enddo
c
c
c..read 2-D profile
c
 
      iu = 7
      do jy=1,kyu
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kyu*kxu + 1 - (jy-1)*kxu - jx )
            read ( iufile, 111, err=92 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx,jy) = zu(iu)
        enddo
      enddo
      if(nsmooth.ne.0)then
        do jy=1,nsmooth
c smooth in space
        do i=1,kyu
          do j=2,kxu-1
            parray(j,i) = 2.D0*parray(j,i)/3.D0 +
     >      (parray(j-1,i)+parray(j+1,i))/6.D0
          enddo
        enddo
c smooth in time
        do j=1,kxu
          do i=2,kyu-1
            parray(j,i) = 2.D0*parray(j,i)/3.D0 +
     >      (parray(j,i-1)+parray(j,i+1))/6.D0
          enddo
        enddo
        enddo
      endif       
c
c
c..end of routine
c
      close ( iufile, err=90 )
      return
c
c
c..error conditions
c
  90  continue
c
c
c...JAK more compact output
crew      write (6,'(a22,a20)') 'cannot find 2d UFILE:',cfile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn u2read'
      kerr = 1
c
c...  JAK 96/04/04 close statement added
c
      close(iufile)
c
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) kxu,' = kxu'
      write (*,*) kyu,' = kyu'
      write (*,*)
     & 'error reading kxu or kyu in UFILE in sbrtn u2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
c
 111  format (1x,6e13.6)
 114  format (1p6e12.4)
 116  format (1p6e13.6)
c
      end
c@u2read.f   17 Nov 1994   Glenn Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u2read ( idat, cdir, cfile, kufile,
     & ppxu, kkxu, ppyu, pdata, kerr, mxtime,
     & xtime, ytime, nty, time_flag, ishift, nsmooth, islice,it,iproc)
c
c   This routine reads 2-D ASCII U-files.
c
c JEK Use either format 111 or 116 for reading ufiles on line 309
c
c  idat      = flag for ufile format                            (input)
c              0 for 111 format, 1 for 116 format
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  pxu(jx)   = 1-D independent variable x array from U-file     (output)
c  kkxu      = number of elements in interpolated x array       (output)
c  mxtime    = max number of time pts                           (input)
c  ppyu       = time of desired profile				(input)
c  kyu       = number of elements in y array from U-file        (output)
c  parray(jx,jy) = 2-D array from U-file                        
c  pdata(jjx) = interpolated radial profile from parray 
c  at time index nt_i 						(output)
c  
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the U-file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c...adapted by J.A.Konings 95/10/25 to also supply time information 
c...added parameters;
c
c   ytime(51,mxtime) data array with max. mxtime timeslices
c   nty            number of actually read timeslices
c   xtime(mxtime)    array containing the time values
c   time_flag      = 0 don t read all timeslices
c                  = 1 fill arrays xtime, ytime 
c   ishift         number of sample points shifted with respect to ppyu
c...JAK 96/05/16 islice added
c   islice         index of array xtime giving the time of interest
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer it, jmaxm, idat, iproc
      integer mxtime, kkxu
      parameter (jmaxm=301)
      real*8 parray(jmaxm,mxtime), pxu(jmaxm), pyu(mxtime), zu(6)
      real*8 pdata(kkxu), ppxu(kkxu), ppyu, der(jmaxm)
      real*8 pxu1, pxu0, parray1, parray0
c
      real*8 xtime(mxtime), ytime(kkxu,mxtime)
      real*8 ytmp(kkxu), ydata(kkxu), sum
      integer time_flag,nty, islice, ishift, nsmooth, ism
c
      integer kxdim, kxu, kydim, kyu, kufile, kdchar, kfchar, kerr
c
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, jy, ju, i, ii
c
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c
c  idchar  = length of directory name cdir
c  ifchar  = length of U-file name cfile
c  iufile  = unit number for reading UFILE
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
      character cdir*(*), cfile*(*), cdfile*130, cline*72
c
c      write(*,*) 'inside u2read ...'
c      write(*,*) 'u2read-ufile = ',cfile
c
c..set defaults
c
      ierr   = kerr
c
      kerr   = 0
c
      iprint = 0
c      iprint=10
c
c      write(*,*) 'mxtime = ',mxtime
c      if(it.eq.29) iprint = 10  ! for diagnostic purposes on a given 2d file
c     7 = NE, 18 = SNBII, 23 = Te
c
      ichmax = 60
c
c kxdim is the maximum number of radial points in the 2D ufiles.
c
      kxdim=jmaxm
c
c kydim is the maximum number of time points in the 2D ufiles.
c
      kydim=mxtime
c
c..clear arrays
c
      do jx=1,kxdim
        pxu(jx)    = 0.D0
      enddo
c
      do jy=1,kydim
        pyu(jy)    = 0.D0
      enddo
c
      do jy=1,kydim
        do jx=1,kxdim
          parray(jx,jy) = 0.D0
        enddo
      enddo
c
c..determine the length of the directory and file names
c  Note:  These names must have no imbedded blanks.
c
      idchar = 0
      do j=1,ichmax
        if ( cdir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 10
        endif
      enddo
  10  continue
c
      ifchar = 0
      do j=1,ichmax
        if ( cfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 11
        endif
      enddo
  11  continue
c
c..protect against no file name given
c
      if ( ifchar .lt. 1 ) then
        write (*,*)
        write (*,*) 'No file name is given for 2-D U-file'
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u2read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u2read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( idchar .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
        ichar  = ifchar
        cdfile = cfile
c
      else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cdir(idchar:idchar) .ne. '/'
     &    .and. cdir(idchar:idchar) .ne. ']' ) then
c
          ichar = idchar + 1 + ifchar
          cdfile = cdir(1:idchar) // '/' // cfile(1:ifchar)
c
        else
c
          ichar  = idchar + ifchar
c
          cdfile = cdir(1:idchar) // cfile(1:ifchar)
c
        endif
c
      endif
c
      if ( iprint .gt. 0 ) then
        write (*,*) cdfile(1:ichar),' = cdfile',it
        write (*,*) ichar,' = idchar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u2read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u2read'
          kerr = 1
          if ( ierr .lt. 0 ) stop
          return
        endif
c
c..open UFILE
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) iufile = 70
c
      open ( iufile, file=cdfile(1:ichar), status='old', err=90 )
c
c
c..read header
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
 100  format (a)
c
c..read number of scalars
c
      read (iufile,102) iscalar
 102  format (i4)
c
c..diagnostic printout
c
      if ( iprint .gt. 8 ) then
        write (*,*)
        write (*,*) cdfile(1:ichar),' = name of dir/file'
        write (*,*) cfile,' = name of UFILE 2'
        write (*,*) iscalar,' = iscalar'
      endif
c
c..read lines of scalars
c
      if ( iscalar .gt. 0 ) then
        do j=1,iscalar*2
          read (iufile,100,err=92) cline
        enddo
      endif
c
c..read variable labels
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94) kxu, kyu
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kkxu,' = kkxu'
        write (*,*) kxu,' = kxu'
        write (*,*) kyu,' = kyu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kyu .gt. kydim ) then
        write (*,*) 'Too many elements in 2-D U-file arrays'
        write (*,*) kxu, kxdim, kyu, kydim
     &   ,' = kxu, kxdim, kyu, kydim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn u2read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x and y 1-D arrays
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
        if(iprint.ge.10) write(*,*) jx, pxu(jx),' pxu'
      enddo
c
      iu = 7
      do jy=1,kyu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kyu + 1 - jy )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pyu(jy) = zu(iu)
        if(iprint.ge.11) write(*,*) jy, pyu(jy),' pyu'
      enddo
c
c
c..read 2-D profile
c
      iu = 7
      do jy=1,kyu
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kyu*kxu + 1 - (jy-1)*kxu - jx )
            if (idat.eq.1) then
               read(iufile,116,err=92) ( zu(ju),ju=1,iumax )
            else
               read(iufile,111,err=92) ( zu(ju),ju=1,iumax )
            endif
            iu = 1
          endif
        parray(jx,jy) = zu(iu)
        if(iprint.ge.12) write(*,*) jy, jx, parray(jx,jy),' parray'
        enddo
      enddo
c
c
c..end of routine
c
      close ( iufile, err=90 )
c
c
c interpolate the iter data to a kxdim point grid and return 
c only the interpolated quantities.
c
c
c   kkxu=51
c

	do i=1,kkxu
	  ppxu(i)=dfloat(i-1)/dfloat(kkxu-1)
	enddo
	ppxu(1)=ppxu(1)+1.D-10
	ppxu(kkxu)=ppxu(kkxu)-1.D-10
	if(iprint.ge.10) then
	  write(*,*) 'time_flag = ',time_flag
	  write(*,*) 'kkxu = ',kkxu
	  do j=1,kkxu
	    write(*,*) j, ppxu(j), ' ppxu'
	  enddo
	endif
c
	i=0
	do ii=1,kyu
	  i=i+1
          if(iprint.ge.10) write(*,*) 'ppyu = ',ppyu,pyu(ii)
	  if(pyu(ii).ge.ppyu) goto 95
	enddo
  95    continue
        if(pyu(ii).gt.ppyu) then
          i=i-1
          pyu(ii)=pyu(ii-1) ! use time-slice where xp_time >= u2d-time
        endif
        if(iprint.ge.10) write(*,*) 'using data at t = ',pyu(ii)
c
c...JAK shift this timeslice and keep thi i-value of the actual timeslice
c...    Do not shift if timeslice is not available
c...JAK dont use ishift
c
ccc     i = i + ishift
ccc     if ((i.lt.1).or.(i.gt.kyu)) i = i - ishift
        islice = i  
c
c...JAK smooth data (only if no time data is read further on)
c
        if (time_flag.eq.0) then 
          nsmooth = min(nsmooth, min(i-1,kyu-i))
          do j=1,kxu
            sum=0.D0
            do ism=-nsmooth,nsmooth
              sum = sum + parray(j,i+ism)
            end do
            parray(j,i) = sum/(2.D0*dfloat(nsmooth)+1.D0)
          end do
        endif
c
c...JAK structure changed: don t return but continue for reading timedaata
c
        if(kxu.ge.kkxu) then
         do j=1,kkxu
          ppxu(j)=pxu(j)
          pdata(j)=parray(j,i)
         enddo
        if(iprint.ge.10) then
          write(*,*)
          do j=1,kxu
            write(*,150) j, pxu(j), parray(j,i)
          enddo
        endif
c
c could be problems here if pxu(1) not 0 or pxu(kxu) not 1.0
cjak     return
        endif
c 
      if(kxu.lt.kkxu) then
c-->
      if(iprint.ge.10) then
        write(*,*)
        do j=1,kxu
          write(*,150) j, pxu(j), parray(j,i)
        enddo
      endif
c
c      pxu0=0.D0     ! jek: disabled 5/22/00
      pxu0=pxu(1)
      if (pxu(1).gt.0.025) pxu0=1.D-3 ! set to zero if pxu(1) too large
      parray0=parray(1,i)+(parray(1,i)-parray(2,i))
     >        /(pxu(1)-pxu(2))*(pxu0-pxu(1))
c      pxu(1)=pxu0   ! jek: disabled 5/22/00
      if (pxu(1).gt.0.025) pxu(1)=pxu0
      parray(1,i)=parray0
  
c     pxu1=1.D0
      pxu1=pxu(kxu)
      parray1=parray(kxu,i)+(parray(kxu,i)-parray(kxu-1,i))
     >        /(pxu(kxu)-pxu(kxu-1))*(pxu1-pxu(kxu))
c      pxu(kxu)=pxu1   ! jek: disabled 5/22/00
      parray(kxu,i)=parray1
      if(iprint.ge.10) then
        do j=1,kxu
          write(*,155) j, pxu(j), parray(j,i)
        enddo
      endif
c
      call INTER_CSPL(kxu,pxu,parray(1,i),kkxu,ppxu,pdata,der)
c
      if(iprint.ge.10) then
        write(*,*)
        do j=1,kkxu
          write(*,66) j, ppxu(j), der(j), pdata(j)
        enddo
      endif
c      write(*,*) 'u2, okay so far-11c...'
c<--
      endif
  66  format(i2,2x,1p3e12.4,2x,'interpolated pdata')
c--------------------------------------------  
c
      if (time_flag.eq.1) then   !fill timearrays
c
c...    reset arrays
c
c     do jy=1,6
c       write(6,*) (parray(jx,jy),jx=1,6)
c     enddo

      do jy=1,mxtime
        xtime(jy) = 0.D0
        do jx=1,kkxu
          ytime(jx,jy) = 0.D0
        enddo
      enddo
c
        nty=kyu  
c
        if(iproc.eq.0) write(6,*) 'Number of timesamples:', kyu
c
        do i=1,kyu           !start time loop
c
          xtime(i)=pyu(i)
c
c
c...      note that ppxu is already calculated at the time ppyu of
c...      interest.
c
          if(kxu.lt.kkxu) then
c
c...        don t adjust the boundary condition for timeslice that
c...        has already been done
c
            if (i.ne.islice) then
              pxu0=0.D0
              parray0=parray(1,i)+(parray(1,i)-parray(2,i))
     >                /(pxu(1)-pxu(2))*(pxu0-pxu(1))
              pxu(1)=pxu0
              parray(1,i)=parray0
 
              pxu1=1.D0
              parray1=parray(kxu,i)+(parray(kxu,i)-parray(kxu-1,i))
     >                /(pxu(kxu)-pxu(kxu-1))*(pxu1-pxu(kxu))
              pxu(kxu)=pxu1
              parray(kxu,i)=parray1
            endif
c
c...        copy to local array
c
            do j=1,kxu
              ytmp(j) = parray(j,i)
            end do
            call INTER_CSPL(kxu,pxu,ytmp(1),kkxu,ppxu,ydata,der)
            do j=1,kkxu
              ytime(j,i) = ydata(j)
            end do
          else
            do j=1,kkxu
              ytime(j,i) = parray(j,i)
            end do
          endif 
        end do
c
      endif
c--------------------------------------------  
c
      return
c
c
c..error conditions
c
  90  continue
c
c
c...JAK more compact output
crew      write (6,'(a22,a20)') 'cannot find 2d UFILE:',cfile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write(*,*) 'zu =',zu(ju)
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn u2read'
      kerr = 1
c
c...  JAK 96/04/04 close statement added
c
      close(iufile)
c
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) kxu,' = kxu'
      write (*,*) kyu,' = kyu'
      write (*,*)
     & 'error reading kxu or kyu in UFILE in sbrtn u2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
c
 111  format (1x,6e13.6)
 114  format (1p6e12.4)
 115  format (1pe14.6,1p5e13.6)
 116  format (1p6e13.6)
 150  format (2x,i2,2x,1p2e13.6,2x'parray-before')
 155  format (2x,i2,2x,1p2e13.6,2x'parray-after')
c
      end
c@u1tread.f   derived from u1read.f by G.M.Staebler 6/28/2000
c reads full time seris from 1d u-files.
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u1tread ( cdir, cfile, kufile, 
     &  kxu ,pxu, parray, nsmooth, kerr )
c
c   This routine reads 1-D ASCII U-files.
c
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  kxdim     = maximum number of elements allowed in arrays     (input)
c
c  kxu       = number of time points            		(input)
c  parray     = data time series         			(output)
c
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the U-file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
      integer nrmax, ntmax, nsmooth
      parameter(nrmax=51) !max number of radial points in experimental grid
c      parameter(ntmax=801) !max number of time points in experimental grid
      parameter(ntmax=3801) !max number of time points in experimental grid
c
      real*8 parray(ntmax), pxu(ntmax), zu(6)
c
      integer kxdim, kxu, kufile, kerr
c
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, ju, i, ii
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c
c  idchar  = length of directory name cdir
c  ifchar  = length of U-file name cfile
c  iufile  = unit number for reading UFILE
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
      character cdir*(*), cfile*(*), cdfile*130, cline*72
c
c
c..set defaults
c
      ierr   = kerr
c
      kerr   = 0
c
      iprint = 0
c
      ichmax = 60
c
c kxdim is the maximum number of time slices we can read
c
      kxdim=ntmax

c
c..clear arrays
c
      do j=1,kxdim
        parray(j) = 0.0
        pxu(j)    = 0.0
      enddo
c
c..determine the length of the directory and file names
c  Note:  These names must have no imbedded blanks.
c
      idchar = 0
      do j=1,ichmax
        if ( cdir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 10
        endif
      enddo
  10  continue
c
      ifchar = 0
      do j=1,ichmax
        if ( cfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 11
        endif
      enddo
  11  continue
c
c..protect against no file name given
c
      if ( ifchar .lt. 1 ) then
        write (*,*)
        write (*,*) 'No file name is given for 1-D U-file'
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u1read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u1read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( idchar .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
        ichar  = ifchar
        cdfile = cfile
c
      else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cdir(idchar:idchar) .ne. '/'
     &    .and. cdir(idchar:idchar) .ne. ']' ) then
c
          ichar = idchar + 1 + ifchar
          cdfile = cdir(1:idchar) // '/' // cfile(1:ifchar)
c
        else
c
          ichar  = idchar + ifchar
c
          cdfile = cdir(1:idchar) // cfile(1:ifchar)
c
        endif
c
      endif
c
      if ( iprint .gt. 0 ) then
        write (*,*) cdfile(1:ichar),' = cdfile'
        write (*,*) ichar,' = ichar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u1read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u1read'
          kerr = 1
          if ( ierr .lt. 0 ) stop
          return
        endif
c
c..open UFILE
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) iufile = 70
c
c         write(*,*) 'file = ',cdfile(1:ichar)
      open ( iufile, file=cdfile(1:ichar), status='old', err=90 )
c
c..read header
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
 100  format (a)
c
c..read number of scalars
c..JAK 96/03/25 error check added
c
      read (iufile,102,err=92) iscalar
 102  format (i4)
c
c..diagnostic printout
c
      if ( iprint .gt. 8 ) then
        write (*,*)
        write (*,*) cdfile(1:ichar),' = name of dir/file'
        write (*,*) cfile,' = name of UFILE 1'
        write (*,*) iscalar,' = iscalar'
      endif
c
c..read lines of scalars
c
      if ( iscalar .gt. 0 ) then
        do j=1,iscalar*2
          read (iufile,100,err=92) cline
        enddo
      endif
c
c..read variable labels
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94) kxu
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kxu,' = kxu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kxu .lt. 1 ) then
        write (*,*) 'Too many elements in 1-D U-file array'
        write (*,*) kxu, kxdim
     &   ,' = kxu, kxdim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kxu .lt. 1 in sbrtn u1read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x  1-D array
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
c       write(*,*) jx, pxu(jx), ' pxu-u1tread'
      enddo
c
c
c..read 1-D profile
c
  
cfix      if(cfile(1:3).ne.'d3d') then 

      iu = 7
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kxu + 1 - jx )
            read ( iufile, 111, err=92 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx) = zu(iu)
        enddo

cfix      else 
cfix
cfix      iu = 7
cfix        do jx=1,kxu
cfix          iu = iu + 1
cfix          if ( iu .gt. 6 ) then
cfix            iumax = min ( 6, kxu + 1 - jx )
cfix            read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
cfix            iu = 1
cfix          endif
cfix        parray(jx) = zu(iu)
cfix        enddo
  

cfix      endif
      if(nsmooth.ne.0)then
c smooth in time
        do i=2,kxu-1
            parray(i) = 2.D0*parray(i)/3.D0 +
     >      (parray(i-1)+parray(i+1))/6.D0
        enddo
      endif       
c
c
c
c..end of routine
c
      close ( iufile, err=90 )
c
      return
c
c
c..error conditions
c
  90  continue
c
c...JAK more compact output
c      write (6,'(a22,a24)') 'cannot find 1d UFILE:',cfile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) kxu,' = kxu'
      write (*,*)
     & 'error reading kxu in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
c
 111  format (1x,6e13.6)
 114  format (1p6e12.4)
 116  format (1p6e13.6)
c
      end
c@u1read.f   17 Nov 1994   Glenn Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u1read ( cdir, cfile, kufile, 
     & ppxu, pdata, kerr )
c
c   This routine reads 1-D ASCII U-files.
c
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  kxdim     = maximum gnumber of elements allowed in arrays     (input)
c
c  ppxu      = time of interest from the experiment		(input)
c  pdata     = data at time of interest 			(output)
c
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the U-file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer ntmax
      parameter(ntmax=3801) !max number of time points in experimental grid
      real*8 parray(ntmax), pxu(ntmax), zu(6), ppxu, pdata
c
      integer kxdim, kxu, kufile, kerr
c
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, ju, i, ii
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c
c  idchar  = length of directory name cdir
c  ifchar  = length of U-file name cfile
c  iufile  = unit number for reading UFILE
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
      character cdir*(*), cfile*(*), cdfile*130, cline*72
c
c
c..set defaults
c
      ierr   = kerr
c
      kerr   = 0
c
      iprint = 0
c
      ichmax = 60
c
c kxdim is the maximum number of time slices we can read
c
      kxdim=ntmax

c
c..clear arrays
c
      do j=1,kxdim
        parray(j) = 0.0
        pxu(j)    = 0.0
      enddo
c
c..determine the length of the directory and file names
c  Note:  These names must have no imbedded blanks.
c
      idchar = 0
      do j=1,ichmax
        if ( cdir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 10
        endif
      enddo
  10  continue
c
      ifchar = 0
      do j=1,ichmax
        if ( cfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 11
        endif
      enddo
  11  continue
c
c..protect against no file name given
c
      if ( ifchar .lt. 1 ) then
        write (*,*)
        write (*,*) 'No file name is given for 1-D U-file'
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u1read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u1read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( idchar .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
        ichar  = ifchar
        cdfile = cfile
c
      else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cdir(idchar:idchar) .ne. '/'
     &    .and. cdir(idchar:idchar) .ne. ']' ) then
c
          ichar = idchar + 1 + ifchar
          cdfile = cdir(1:idchar) // '/' // cfile(1:ifchar)
c
        else
c
          ichar  = idchar + ifchar
c
          cdfile = cdir(1:idchar) // cfile(1:ifchar)
c
        endif
c
      endif
c
      if ( iprint .gt. 0 ) then
        write (*,*) cdfile(1:ichar),' = cdfile'
        write (*,*) ichar,' = idchar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u1read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u1read'
          kerr = 1
          if ( ierr .lt. 0 ) stop
          return
        endif
c
c..open UFILE
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) iufile = 70
c
      open ( iufile, file=cdfile(1:ichar), status='old', err=90 )
c
c..read header
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
 100  format (a)
c
c..read number of scalars
c..JAK 96/03/25 error check added
c
      read (iufile,102,err=92) iscalar
 102  format (i4)
c
c..diagnostic printout
c
      if ( iprint .gt. 8 ) then
        write (*,*)
        write (*,*) cdfile(1:ichar),' = name of dir/file'
        write (*,*) cfile,' = name of UFILE 1'
        write (*,*) iscalar,' = iscalar'
      endif
c
c..read lines of scalars
c
      if ( iscalar .gt. 0 ) then
        do j=1,iscalar*2
          read (iufile,100,err=92) cline
        enddo
      endif
c
c..read variable labels
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94) kxu
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kxu,' = kxu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kxu .lt. 1 ) then
        write (*,*) 'Too many elements in 1-D U-file array'
        write (*,*) kxu, kxdim
     &   ,' = kxu, kxdim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kxu .lt. 1 in sbrtn u1read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x and y 1-D arrays
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
      enddo
c
c
c..read 1-D profile
c
  
cfix      if(cfile(1:3).ne.'d3d') then 

      iu = 7
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kxu + 1 - jx )
            read ( iufile, 111, err=92 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx) = zu(iu)
        enddo

cfix      else 
cfix
cfix      iu = 7
cfix        do jx=1,kxu
cfix          iu = iu + 1
cfix          if ( iu .gt. 6 ) then
cfix            iumax = min ( 6, kxu + 1 - jx )
cfix            read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
cfix            iu = 1
cfix          endif
cfix        parray(jx) = zu(iu)
cfix        enddo
  

cfix      endif
c
	i=0
	do ii=1,kxu
	  i=i+1
	  if(pxu(ii).ge.ppxu) goto 95
	enddo
  95    continue
	pdata=parray(i)	
c
c..end of routine
c
      close ( iufile, err=90 )
c
      return
c
c
c..error conditions
c
  90  continue
c
c...JAK more compact output
crew      write (6,'(a22,a20)') 'cannot find 1d UFILE:',cfile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) kxu,' = kxu'
      write (*,*)
     & 'error reading kxu in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
c
 111  format (1x,6e13.6)
 114  format (1p6e12.4)
 116  format (1p6e13.6)
c
      end
c@u0read.f
c  Read in data in ITER/ITPA PDB 0d ufiles
c  idat = 0 ITER PDB 0D ufile format
c  idat = 1 new IPTA PDB format (163 entries) 3/27/07
c           Note: must remove any slashes
c  from show0d.f from Stan Attenberger, ORNL, June 95.
c  adapted by Joop Konings 04/95, J. Kinsey 9/06
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u0read(idat,udfile,name,value,ifail)
c
      implicit none
c
      integer ifail, nun, mxshot, mxrow, mxcol, nnames, i, ishot
      integer jdone, inames, nline, jcol
      integer idat, mrow, ju, imax
      character*100 udfile
      character*67 name(200)
      character*11 value(200)
      character*80 line
      character*132 line2
      data nun/10/, mrow/25/
c     data mxshot/50/,mxrow/7/
      data mxshot/1/,mxrow/16/,mxcol/12/,nnames/200/
c      include 'mpif.h'
      include 'glf.m'
c
c     open(unit=nun,file='readnames_0d.dat',status='old',err=888)
c     read(nun,'(i5)') nnames
c     read(unit=nun,fmt='(a67)',end=100,err=888) (name(i),i=1,nnames)
c     close(nun)
      write(*,500) udfile
      open(unit=nun,file=udfile,status='old',err=888)
c
      if(idat.eq.1.or. idat.eq.3) then
      imax=163
      do ishot=1,mxshot
        do inames=1,mrow
          read(nun, *,end=101 ) (value(ju),ju=1,imax)
c          write(*,*) ( value(ju),ju=1,imax )
        enddo
      enddo
c
      elseif(idat.eq.2) then
c      write(*,*) 'mxshot = ',mxshot,mrow
      imax=164
      do ishot=1,mxshot
        do inames=1,mrow
          read(nun, *,end=101 ) (value(ju),ju=1,imax)
c          write(*,*) ( value(ju),ju=1,imax )
        enddo
      enddo
c
      else
c
      do ishot=1,mxshot
        jdone=0
c       do inames=1,nnames,mxrow
        do inames=1,mxrow
          read(unit=nun,fmt='(a80)',end=101) line
          nline=80
          do jcol=1,mxcol
            jdone=jdone+1
            call part(line,nline,value(jdone))
c            write(*,*) inames,jdone,'value = ',value(jdone)
            if(jdone.eq.nnames)go to 50
          enddo
        enddo
 50     continue
      enddo
      endif
c
c*****normal exit
c
 101  return
c
c*****error exit
c
 888  continue
      ifail=1
      return
c
 100  write(*,*) 'fatal: expected ',nnames,' lines in names_0d'
      goto 888
 500  format('0D file =',a50)
      end
      subroutine inter_cspl(n,r,data,m,x,ds)
c
c... INTERPOLATE DATA TO NEW GRID
c
      implicit none
c
      integer j, n, m, ifail
      real*8 r(n), data(n), x(m), ds(m), br(n)
      logical sk
c
      sk = .TRUE.
c 
      ifail=0
c     call e01bee(n,r,data,br,ifail)
c     call e01bfe(n,r,data,br,m,x,ds,ifail)
c
      call dpchim(n,r,data,br,1,ifail)
      call dpchfe(n,r,data,br,1,sk,m,x,ds,ifail)
c
      return
      end
c@datavg.f
c jek 19-Jan-11 version 2.0
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Smoothes data via 7-pt averaging if(ismooth_data.ne.0)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine datavg
c
      implicit none
      include 'data_d.m'
      include 'input.m'
c
      integer j, jj, jn
c
       call average7_1d(fcap_d,nj_d)
       call average7_1d(fcap_d,nj_d)
       call average7_1d(gcap_d,nj_d)
       call average7_1d(gcap_d,nj_d)
       call average7_1d(hcap_d,nj_d)
       call average7_1d(hcap_d,nj_d)
c
       call average7_1d(te_d,nj_d)
       call average7_1d(ti_d,nj_d)
c
       call average7_1d(q_d,nj_d)
       call average7_1d(q_d,nj_d)
c
       call average7_1d(ene_d,nj_d)
       do jj=1,nion_d
         call average7_1d(en_d(1,jj),nj_d)
       enddo
       do jj=1,nprim_d
          call average7_1d(sion_d(1,jj),nj_d)
          call average7_1d(srecom_d(1,jj),nj_d)
          call average7_1d(scx_d(1,jj),nj_d)
          call average7_1d(sbcx_d(1,jj),nj_d)
          call average7_1d(s_d(1,jj),nj_d)
          call average7_1d(dudtsv_d(1,jj),nj_d)
       enddo
       do jn=1,nneu_d
          call average7_1d(enn_d(1,jn),nj_d)
          call average7_1d(ennw_d(1,jn),nj_d)
          call average7_1d(ennv_d(1,jn),nj_d)
c          call average7_1d(volsn_d(1,jn),nj_d)
       enddo
       call average7_1d(enbeam_d,nj_d)
       call average7_1d(enbeam_d,nj_d)
c
       call average7_1d(sbion_d,nj_d)
       call average7_1d(sbeam_d,nj_d)
c
       call average7_1d(curden_d,nj_d)
       call average7_1d(curohm_d,nj_d)
       call average7_1d(curboot_d,nj_d)
       call average7_1d(curbeam_d,nj_d)
       call average7_1d(currf_d,nj_d)
c
       call average7_1d(rbp_d,nj_d)
c
       call average7_1d(zeff_d,nj_d)
       call average7_1d(zeff_d,nj_d)
c
       call average7_1d(angrot_d,nj_d)
c
       if(ismooth_data.ge.2) call average7_1d(dpedtc_d,nj_d)
       if(ismooth_data.ge.2) call average7_1d(dpidtc_d,nj_d)
c
       call average7_1d(qdelt_d,nj_d)
       call average7_1d(qbeame_d,nj_d)
       call average7_1d(qbeami_d,nj_d)
       call average7_1d(qrfe_d,nj_d)
       call average7_1d(qrfi_d,nj_d)
       call average7_1d(qfuse_d,nj_d)
       call average7_1d(qfusi_d,nj_d)
       if(ismooth_data.ge.2) call average7_1d(qrad_d,nj_d)
       if(ismooth_data.ge.2) call average7_1d(qohm_d,nj_d)
       if(ismooth_data.ge.2) call average7_1d(qione_d,nj_d)
       if(ismooth_data.ge.2) call average7_1d(qioni_d,nj_d)
       if(ismooth_data.ge.2) call average7_1d(qcx_d,nj_d)
c
       call average7_1d(rmajavnpsi_d,nj_d)
       call average7_1d(rmajavnpsi_d,nj_d)
c
       call average7_1d(rminavnpsi_d,nj_d)
       call average7_1d(rminavnpsi_d,nj_d)
c
       call average7_1d(psivolp_d,nj_d)
       call average7_1d(psivolp_d,nj_d)
c
       call average7_1d(elongx_d,nj_d)
       call average7_1d(elongx_d,nj_d)
c
       call average7_1d(deltax_d,nj_d)
       call average7_1d(deltax_d,nj_d)
c
       call average7_1d(sfareanpsi_d,nj_d)
       call average7_1d(sfareanpsi_d,nj_d)
c
       call average7_1d(cxareanpsi_d,nj_d)
       call average7_1d(cxareanpsi_d,nj_d)
c      
       call average7_1d(grho1npsi_d,nj_d)
       call average7_1d(grho1npsi_d,nj_d)
c
       call average7_1d(grho2npsi_d,nj_d)
       call average7_1d(grho2npsi_d,nj_d)
c
       if(itorque.ne.0 .or. iptotr.eq.2)then
         call average7_1d(torque_d,nj_d)
       endif
c
       call average7_1d(ptot_d,nj_d)
       if(iptotr.eq.2) then
         call average7_1d(pfast_d,nj_d)
       endif
c
       if(ncl_flag.eq.1) then
         call average7_1d(xb2_d,nj_d)
         call average7_1d(xbm2_d,nj_d)
         call average7_1d(xb2_d,nj_d)
         call average7_1d(xbm2_d,nj_d)
         call average7_1d(xngrth_d,nj_d)
         call average7_1d(xgrbm2_d,nj_d)
         call average7_1d(fm1_d,nj_d)
         call average7_1d(fm2_d,nj_d)
         call average7_1d(fm3_d,nj_d)
         call average7_1d(fhat_d,nj_d)
       endif
c
      end
      subroutine average7_1d(f,nk)
c*****************************************************
c
c performs a two pass seven point average 
c
c*****************************************************
      implicit none
c
      integer nk,k,i,m
      real*8 f(nk),g(nk)
c
c check if grid is too small
      if(nk.lt.7)return
c
        g(1)=f(1)
        g(nk)=f(nk)
        k=2
        g(k)=(f(k)+f(k+1)+f(k-1))/3.D0
        k=3
        g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.D0
        k=nk-1
        g(k)=(f(k)+f(k+1)+f(k-1))/3.D0
        k=nk-2
        g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.D0
        do k=4,nk-3
          g(k)=(f(k)+f(k-1)+f(k+1)+f(k-2)+f(k+2)
     >          + f(k+3)+f(k-3))/7.D0
        enddo
        k=2
        f(k)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=3
        f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        k=nk-1
        f(k)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=nk-2
        f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        do k=4,nk-3
          f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2)
     >            + g(k+3)+g(k-3))/7.D0
        enddo
      return
      end
c******************************************************
      subroutine part(work,largs,strnxt)
c******************************************************
c*****part finds the next item in the string.
c*****Input:
c*****work-a character variable.
c*****largs-the number of characters to be parsed
c*****Output:
c*****strnxt-the next item in the list,
c*****       delimited by matched blanks, "", or {}
c*****Note: On return, work and largs have been modified to reflect
c*****      the removal of the first item.
c*****last revision: 12/94 s.e.attenberger, w.a.houlberg, ornl
c*******************************************************
      implicit none
c
      integer largs, ltst, lenst, ntst
      character*(*) work,  strnxt
      character*1 find,    tab,    wtst
      character*1 space,      quote,      lbr,     rbr
      data        space/' '/, quote/'"'/, lbr/'{'/, rbr/'}'/
      tab=char(9)
c
      strnxt=' '
      if(largs.eq.0) return
c*****look for start of item (non-blank character)
      ltst=0
100   ltst=ltst+1
      if    ((work(ltst:ltst).eq.space.or.work(ltst:ltst).eq.tab)
     >        .and. ltst.lt.largs) then
        go to 100
      elseif((work(ltst:ltst).eq.space.or.work(ltst:ltst).eq.tab)
     >        .and. ltst.eq.largs) then
c*****  no items in this list
        strnxt=' '
        work=' '
        largs=0
        return
      elseif(work(ltst:ltst).eq.quote) then
        find=quote
        ltst=ltst+1
      elseif(work(ltst:ltst).eq.lbr) then
        find=rbr
        ltst=ltst+1
      else
        find=space
      endif
      lenst=1
      strnxt(lenst:lenst)=work(ltst:ltst)
c*****Start of item is character ltst.  Now search for end of item.
      ntst=ltst
200   ntst=ntst+1
c*****Treat tab like space for testing, but dont replace tab in work.
      wtst=work(ntst:ntst)
      if(wtst.eq.tab) wtst=space
      if(ntst.gt.largs.and.find.ne.space) then
        write(*,*)' fatal error, missing ///', find,'/// near'
        write(*,*)work(1:largs)
        stop
      elseif(ntst.le.largs.and.wtst.ne.find) then
        lenst=lenst+1
        strnxt(lenst:lenst)=work(ntst:ntst)
        go to 200
      elseif(ntst.gt.largs.and.find.eq.space) then
c*****  successful exit, end of string.
c*****  (no space delimiter is required at the end of a string)
        work(1:largs)=' '
        largs=0
      elseif(wtst.eq.find)then
c*****  successful exit
        work(1:largs-ntst)=work(ntst+1:largs)
        work(largs-ntst+1:largs)=' '
        largs=largs-ntst
      endif
c
      return
      end
!
      SUBROUTINE W_LIN_INTERP(n1,x1,y1,n2,x2,y2,iflag,message)
!***********************************************************************
!W_LIN_INTERP does a linear interpolation to obtain n2 values for the
!  target array y2 on the grid x2 from the n1 values of y1 on the grid
!  x1
!References:
!  W.A.Houlberg 3/2000
!Input:
!  n1-number of abscissas and values in the source arrays
!  x1-array of source abscissas
!  y1-array of source values
!  n2-number of target values to be found
!  x2-array of target abscissas
!Output:
!  y2-array of target values
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        n1,                      n2
      REAL           x1(*),                   y1(*),
     &               x2(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           y2(*)
!Declaration of local variables
      INTEGER        i,                       il
      IF(n1.lt.2) THEN
        iflag=1
        message='W_LIN_INTERP/ERROR:less than 2 points in source array'
      ELSE
        il=1
        DO i=1,n2
   10     IF(x2(i).lt.x1(1)) THEN
!           Use innermost data value
            y2(i)=y1(1)
            iflag=-1
            message='W_LIN_INTERP(1)/WARNING:x<x(1), use end point'
          ELSEIF(x2(i).eq.x1(1)) THEN
            y2(i)=y1(1)
          ELSEIF(x2(i).gt.x1(il+1)) THEN
            IF(il.lt.n1-1) THEN
!             Step x1 grid forward and loop
              il=il+1
              GOTO 10
            ELSE
!             Set to last value
              y2(i)=y1(n1)
              iflag=-1
              message='W_LIN_INTERP(2)/WARNING:x>x(n1), use end point'
            ENDIF
          ELSE
!           Interpolate
            y2(i)=y1(il)
     &            +(y1(il+1)-y1(il))*(x2(i)-x1(il))/(x1(il+1)-x1(il))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
   
