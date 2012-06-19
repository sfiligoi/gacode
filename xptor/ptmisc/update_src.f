cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_src
c
c Update experimental sources
c   upsrc(1) : q-profile
c   upsrc(2) : power flows

      implicit none
c
      include '../inc/ptor.m'
      include '../inc/vtrans.m'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/glf.m'
c
c nrmax = max number of radial grid points in experimental data
c nr   = actual number of radial grid points in experimental data
c xp_time = time of desired profile (input)
c ismooth = number of time over which profiles are smoothed (input)
c islice_d = index of array xtime_d giving the time of interest
c
      integer nrmax, n2d
c     parameter(nrmax=51) ! max number of radial pts in experimental grid
      parameter(n2d=44)   ! max number of 2D arrays
      character*7 fields2d(n2d)
      character*50 ufile
      integer j, k, ifd, ierr, ig, is, i2d, linblnk, nr
      integer time_flag, itime_shift, ierr_ei, ierr_rf, kgrid
      real*8 xtime_d, ytime_d, time_ptor
      real*8 zhalf, drho, dvoldr_m, dvoldr_p, zbrac
      real*8 alamda, tau_e
      real*8 rho_xp(mxgrid+1), u2d(mxgrid+1,n2d)
      real*8 powe_p, powe_mm, powi_p, powi_mm, volintk
c
      data fields2d/'CHIE   ','CHII   ','CURNBI ','CURTOT ',
     .	'GRHO1  ','GRHO2  ','NE     ','NFAST  ','NIMP   ','NM1    ',
     .	'QNBIE  ','QNBII  ','QOHM   ','QRAD   ','RMAJOR ','RMINOR ',
     .	'SNBIE  ','SNBII  ','SURF   ','ZEFFR  ','VOLUME ','TI     ',
     .	'TE     ','DELTAR ','INDENTR','KAPPAR ','Q      ','QEI    ',
     .	'VROT   ','SWALL  ','DWER   ','DWIR   ','QECHE  ','QICRHE ',
     .	'QECHI  ','QICRHI ','QWALLE ','QWALLI ','QFUSE  ','QFUSI  ',
     .	'DNER','NM2 ','NM3 ','PTOTR'/
c
      if(i_proc.eq.0) write(*,*) 'Updating sources ...'
      nrmax=mxgrid+1
      jmaxm=mxgrid
c      ig=linblnk(tok,10)
c      is=linblnk(shot,10)
      ig=linblnk(tok)
      is=linblnk(shot)
      time_flag=0  ! don't read all timeslices
      time_ptor=xp_time+(time_upd-time)
      zhalf = 1.D0/2.D0
      ierr_ei=0  ! don't warn about missing exchange data
      ierr_rf=0  ! don't warn about missing RF data
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... q (i.e. safety factor) and shear profiles
c
      if (upsrc(1) .ge. 1) then
        ifd = 27
c        i2d=linblnk(fields2d(ifd),7)
        i2d=linblnk(fields2d(ifd))
        ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields2d(ifd)(1:i2d)
c        write(*,*) xp_time, time, time_ptor
c
c...read data
c
        ierr = 0
        nr=mxgrid+1
        call u2read(cudir,ufile,88,rho_xp,nr,time_ptor,
     &              u2d(1,ifd),ierr,xtime_d, ytime_d, ntime_d,
     &              time_flag,itime_shift,ismooth, islice_d)
        if(ierr.ge.1) then
            if(i_proc.eq.0) write(*,*) '  Error reading ufiles for ',
     &                      fields2d(ifd)(1:i2d),ierr
        endif
c
      do j=1,nrmax
        q_d(j)=u2d(j,27)
      enddo
      do j=1,jmaxm
        q_exp(j-1)=q_d(j)
      enddo
      do j=1,ngrid
        qpro(j)=q_exp(j)
      enddo
      qa=q_exp(jmaxm)
      q0=q_exp(0)
c
      do j=1,jmaxm
        shat_exp(j)=(rho(j)+rho(j-1))/(q_exp(j)+q_exp(j-1))*
     >              (q_exp(j)-q_exp(j-1))/(rho(j)-rho(j-1))
c       write(*,100) j, rho(j), q_exp(j), shat_exp(j)
      enddo
c
      endif
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Electron and ion power flows
c    Electron: Beam, RF, Ohmic, Radiation
c    Ion: Beam, RF
c
      if (upsrc(2) .ge. 1) then
c
c---power to electrons and ions from beam, watts/m**3 

        do j=11,12
          ierr = 0
          ifd = j
          nr=mxgrid+1
c          i2d=linblnk(fields2d(ifd),7)
          i2d=linblnk(fields2d(ifd))
          ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields2d(ifd)(1:i2d)
          call u2read(cudir,ufile,88,rho_xp,nr,time_ptor,u2d(1,ifd),
     &                ierr,xtime_d, ytime_d, ntime_d, time_flag,
     &                itime_shift, ismooth, islice_d)
          if(ierr.ge.1 .and. i_proc.eq.0) then
             write(*,*) '  Error reading ufiles for ',
     &                  fields2d(ifd)(1:i2d),ierr
          endif
        enddo
c
        do j=1,nrmax
          qbeame_d(j)=u2d(j,11)        !watts/meter**3
          qbeami_d(j)=u2d(j,12)        !watts/meter**3
        enddo
c
c---qrfe, rf electron heating, watts/m**3
c
        do j=33,36
          ifd = j
          nr=mxgrid+1
c          i2d=linblnk(fields2d(ifd),7)
          i2d=linblnk(fields2d(ifd))
          ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields2d(ifd)(1:i2d)
          call u2read(cudir,ufile,88,rho_xp,nr,time_ptor,u2d(1,ifd),
     &                ierr,xtime_d, ytime_d, ntime_d, time_flag,
     &                itime_shift, ismooth, islice_d)
          if(ierr.ge.1 .and. ierr_rf.ne.0 .and. i_proc.eq.0) then
             write(*,*) '  Error reading ufiles for ',
     &                  fields2d(ifd)(1:i2d),ierr
          endif
        enddo
c
        do j=1,nrmax
          qrfe_d(j)=u2d(j,33) + u2d(j,34)  !watts/meter**3
          qrfi_d(j)=u2d(j,35) + u2d(j,36)
        enddo
c
c---(electron) ohmic power density, watts/m**3
c---radiated power density, watts/m**3
c   
        do j=13,14
          ierr = 0
          ifd = j
          nr=mxgrid+1
c          i2d=linblnk(fields2d(ifd),7)
          i2d=linblnk(fields2d(ifd))
          ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields2d(ifd)(1:i2d)
          call u2read(cudir,ufile,88,rho_xp,nr,time_ptor,u2d(1,ifd),
     &                ierr,xtime_d, ytime_d, ntime_d, time_flag,
     &                itime_shift, ismooth, islice_d)
          if(ierr.ge.1 .and. i_proc.eq.0) then
             write(*,*) '  Error reading ufiles for ',
     &                  fields2d(ifd)(1:i2d),ierr
          endif
        enddo
c
        do j=1,nrmax
          qohm_d(j)=u2d(j,13)        !watts/meter**3
          qrad_d(j)=u2d(j,14)        !watts/meter**3
        enddo
c
c---qdelt : electron-ion equilibration, watts/m**3
c
        ifd = 28
c        i2d=linblnk(fields2d(ifd),7)
        i2d=linblnk(fields2d(ifd))
        ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields2d(ifd)(1:i2d)
        ierr = 0
        nr=mxgrid+1
        call u2read(cudir,ufile,88,rho_xp,nr,time_ptor,u2d(1,ifd),
     &              ierr,xtime_d, ytime_d, ntime_d, time_flag,
     &              itime_shift, ismooth, islice_d)
        if(ierr.ge.1 .and. ierr_ei.ne.0 .and. i_proc.eq.0) then
            write(*,*) '  Error reading ufiles for ',
     &                 fields2d(ifd)(1:i2d),ierr
        endif
c
        do j=1,nrmax
          qdelt_d(j)=u2d(j,28)
        enddo
        if(qdelt_d(nrmax).eq.0.or.iexp_exch.ge.1) then
          do j=1,nrmax
           zbrac=en_d(j,1)/ene_d(j)+apgasa/
     >         apimpa*en_d(j,2)/ene_d(j)*apimpz**2
           alamda=24.D0-dlog((dabs(ene_d(j))*1.D-6)**zhalf/
     >         (dabs(te_d(j))*1.D3))
           tau_e=3.44D5*(dabs(te_d(j))*1.D3)**1.5D0/
     >         (ene_d(j)*1.D-6)/alamda
           qdelt_d(j)=-10**6*1.6022D-19*zbrac*3.D0/(1836.D0*apgasa)/
     >         tau_e*ene_d(j)*1.D-6*(te_d(j)-ti_d(j))*1.D3
         enddo                    
        if(iexp_exch.eq.-1) then
          do j=1,nrmax
            qdelt_d(j)=0.D0
          enddo
        endif
        endif
c
c---integrated power densities
c
        do j=2,nrmax
          drho=r_d(j)-r_d(j-1)
          dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
          dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
          powe_beam_exp(j-1)=powe_beam_exp(j-2) + 1.D-6*zhalf*
     &         (dvoldr_p*qbeame_d(j)+dvoldr_m*qbeame_d(j-1))*drho
          powi_beam_exp(j-1)=powi_beam_exp(j-2)+ 1.D-6*zhalf*
     &         (dvoldr_p*qbeami_d(j)+dvoldr_m*qbeami_d(j-1))*drho
          powe_rf_exp(j-1)=powe_rf_exp(j-2)+ 1.D-6*zhalf*
     &         (dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*drho
          powi_rf_exp(j-1)=powi_rf_exp(j-2)+ 1.D-6*zhalf*
     &         (dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*drho
          powe_oh_exp(j-1)=powe_oh_exp(j-2)+ 1.D-6*zhalf*
     &         (dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drho
          powe_rad_exp(j-1)=powe_rad_exp(j-2)+ 1.D-6*zhalf*
     &         (dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*drho
c
          pow_ei_exp(j-1)=pow_ei_exp(j-2)-1.D-6*zhalf*
     &         (dvoldr_p*qdelt_d(j)+dvoldr_m*qdelt_d(j-1))*drho
c
          powe_exp(j-1)=powe_beam_exp(j-1)+
     &      powe_rf_exp(j-1)+xoh_exp*powe_oh_exp(j-1)
     &      -xrad_exp*powe_rad_exp(j-1)-powe_ion_exp(j-1)-
     &      xwdot*powe_wdot_exp(j-1)
     &      -pow_ei_exp(j-1)+xfus_exp*powe_fus_exp(j-1)
          powi_exp(j-1)=powi_beam_exp(j-1)+
     &      powi_rf_exp(j-1)
     &      -powi_ion_exp(j-1)+powi_cx_exp(j-1)-
     &      xwdot*powi_wdot_exp(j-1)
     &      +pow_ei_exp(j-1)+xfus_exp*powi_fus_exp(j-1)
c        write(*,100) j, rho(j-1),powe_exp(j-1),powi_exp(j-1)
        enddo
c
c Need to add wdot, ionization, cx, flows
c
      do k=1,ngrid
         powe_p=powe_exp(k)
         powe_mm=powe_exp(k-1)
         powi_p=powi_exp(k)
         powi_mm=powi_exp(k-1)
         if (iexch.ge.1) then
           powe_p=powe_p+pow_ei_exp(k)
           powe_mm=powe_mm+pow_ei_exp(k-1)
           powi_p=powi_p-pow_ei_exp(k)
           powi_mm=powi_mm-pow_ei_exp(k-1)
         endif
         if (iohm.ge.1) then
           powe_p=powe_p-powe_oh_exp(k)
           powe_mm=powe_mm-powe_oh_exp(k-1)
         endif
         Peaux(k)=(powe_p-powe_mm)/vprime(k,1)/dr(k,1) 
         Piaux(k)=(powi_p-powi_mm)/vprime(k,1)/dr(k,1)
      enddo
c
c...compute psume and psumi needed by trcoef
c
      do k=1,ngrid
        kgrid=k
        psume(k)=volintk(Peaux,kgrid)
        psumi(k)=volintk(Piaux,kgrid)
      enddo
c
c Now add stuff from vtr_pack ...
c
      endif
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
100   format(i2,2x,0p1f4.2,0p6f10.5)
c
      return
      end

