!
!--------------------------------------------------------------
!
      SUBROUTINE get_multiscale_spectrum
!
!***********************************************************************
!  questions  should be addressed to
!  Gary Staebler 858-455-3466  or email: gary.staebler@gat.com
!***********************************************************************
!  TGLF multiscale saturation rule: sat_rule_in = 1 option
!
!***********************************************************************
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      USE tglf_xgrid
      IMPLICIT NONE
      !
      LOGICAL :: USE_MIX=.TRUE.
      LOGICAL :: USE_TTF=.FALSE.
      INTEGER :: i,is,k,j,j1,jmax1,jmax2
      REAL :: test,testmax1,testmax2
      REAL :: gammamax1,kymax1,gammamax2,kymax2,ky0,ky1,ky2
      REAL :: f0,f1,f2,a,b,c,x0,x02,dky,xmax
      REAL :: gamma0,gammaeff
      REAL :: cnorm, phinorm, kylow, czf, cz1, cz2, kyetg
      REAL :: cky,sqcky,delta,ax,ay,kx
      REAL :: mix1,mix2,mixnorm,gamma_ave
      REAL :: vzf,dvzf,vzf1,vzf2,bz1,bz2
      REAL,DIMENSION(nkym) :: gamma_net=0.0
      REAL,DIMENSION(nkym) :: gamma=0.0
      REAL,DIMENSION(nkym) :: gamma_mix=0.0
      REAL,PARAMETER :: small=1.0E-10
      !
      ! model fit parameters
      ! need to set alpha_zf_in = 1.0
      ! Miller geometry values igeo=1
      czf = alpha_zf_in
      bz1=0.0
      bz2=0.0
      kyetg=1.28
      cnorm=14.21
!      cz1=0.48*czf
      cz1 = 0.48
      cz2=1.0*czf
      ax=0.0
      ay=0.0
      if(alpha_quench_in.eq.0.0)then
      !spectral shift model parameters
        ax = 1.15
        ay = 0.56
      endif
      !   write(*,*)" ax= ",ax," ay= ",ay
      do i=1,nky
         kx=spectral_shift_out(i)
         gamma_net(i) = eigenvalue_spectrum_out(1,i,1)/(1.0 + (ax*kx)**4)
      !   write(*,*)i,"gamma_net = ",gamma_net(i)
      enddo
      if(USE_MIX)then
        kyetg=1.9
        cky=3.0
        sqcky=SQRT(cky)
        cnorm = 14.29
!        cz1=0.48*czf
        cz1 = 0.48
        cz2=0.92*czf  
        if(USE_TTF)then
           bz1=0.12
           bz2=0.12
           cz2=1.2*czf
           kyetg=1.0
           cky=2.0
          sqcky=SQRT(cky)
        endif  
      endif    
      kyetg = kyetg*ABS(zs(2))/SQRT(taus(2)*mass(2))
      if(igeo.eq.0)then ! s-alpha 
       cnorm=14.63
       cz1=0.90*czf
       cz2=1.0*czf
      endif
      !
      ! renormalize the fluxes and intensities to the phi-norm from the v-norm
      do j=1,nky
      !  write(*,*)"spectal_shift_out(",j,") = ",spectral_shift_out(j)
         do i=1,nmodes_in
            phinorm=1.0
            if(ABS(field_spectrum_out(2,j,i)).gt.small)phinorm=field_spectrum_out(2,j,i)
            do is=1,ns
               intensity_spectrum_out(1,is,j,i) = intensity_spectrum_out(1,is,j,i)/phinorm
               intensity_spectrum_out(2,is,j,i) = intensity_spectrum_out(2,is,j,i)/phinorm
               do k=1,3
                  flux_spectrum_out(1,is,k,j,i) = flux_spectrum_out(1,is,k,j,i)/phinorm
                  flux_spectrum_out(2,is,k,j,i) = flux_spectrum_out(2,is,k,j,i)/phinorm
                  flux_spectrum_out(3,is,k,j,i) = flux_spectrum_out(3,is,k,j,i)/phinorm
                  flux_spectrum_out(4,is,k,j,i) = flux_spectrum_out(4,is,k,j,i)/phinorm
                  flux_spectrum_out(5,is,k,j,i) = flux_spectrum_out(5,is,k,j,i)/phinorm
              enddo
           enddo
        enddo
      enddo
      ! find the maximum of gamma/ky 
      gammamax1= gamma_net(1)
      kymax1 = ky_spectrum(1)
      testmax2 = gammamax1/kymax1
      testmax1 = testmax2
      jmax1=1
      jmax2=nky-2
      kylow=MIN(0.8/SQRT(taus_in(2)/mass_in(2)),ky_spectrum(nky-1))
      j1=0
      do j=2,nky
         ky0 = ky_spectrum(j)
         if(ky0 .lt. kylow)j1=j1+1
         test = gamma_net(j)/ky0
         if(test .gt. testmax2)then
           testmax2 = test
           jmax2=j
           if(ky0 .lt. kylow)then
             testmax1=testmax2 
             jmax1=jmax2
           endif      
         endif        
      enddo
!      if(jmax1.lt.nky)then
!        test=gamma_net(jmax1+1)/ky_spectrum(jmax1+1)
!        if(testmax1.le.test)then
          ! there is no low-k peak 
!          write(*,*)" gammamax1/kymax1 = ",gammamax1/kymax1, test
!          jmax1=jmax2  
!        endif
!      endif
      gammamax2 = gamma_net(jmax2)
      kymax2 = ky_spectrum(jmax2)
      gammamax1 = gamma_net(jmax1)
      kymax1 = ky_spectrum(jmax1)
!      write(*,*)" jmax1 = ",jmax1," jmax2= ",jmax2
!      write(*,*)" g/k 1 = ",gammamax1/kymax1," g/k 2 = ",gammamax2/kymax2
      !interpolate to find a more accurate low-k maximum gamma/ky 
      ! this is cut of at j1 since a maximum may not exist in the low-k range
      if(jmax1.gt.1.and.jmax1.lt.j1)then
         f0 =  gamma_net(jmax1-1)/ky_spectrum(jmax1-1)
         f1 =  gamma_net(jmax1)/ky_spectrum(jmax1)
         f2 =  gamma_net(jmax1+1)/ky_spectrum(jmax1+1)
         dky = (ky_spectrum(jmax1+1)-ky_spectrum(jmax1-1))
         x0 = (ky_spectrum(jmax1)-ky_spectrum(jmax1-1))/dky
         a = f0
         x02 = x0*x0
         b = (f1 - f0*(1-x02)-f2*x02)/(x0-x02)
         c = f2 - f0 - b
         xmax = -b/(2.0*c)
         if(xmax .ge. 1.0)then
           kymax1 = ky_spectrum(jmax1+1)
           gammamax1 = f2*kymax1
         elseif(xmax.lt.0.0)then
           kymax1 = ky_spectrum(jmax1-1)
           gammamax1 = f0*kymax1
         else
           kymax1 = ky_spectrum(jmax1-1)+dky*xmax
           gammamax1 = (a+b*xmax+c*xmax*xmax)*kymax1
         endif     
      endif
!      write(*,*)"gammamax1 = ",gammamax1," kymax1 = ",kymax1," kylow = ",kylow
      ! compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
      ! note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
      vzf1 = gammamax1/kymax1
      vzf2 = gammamax2/kymax2
      dvzf = MAX(vzf2-vzf1,0.0)
!      write(*,*)"dvzf = ",dvzf," vzf1 = ",vzf1," vzf2 = ",vzf2
      vzf = vzf1 + bz1*dvzf
      do j=1,nky
! include zonal flow effects on growthrate model
!          gamma=0.0
          gamma0 = gamma_net(j)
          ky0=ky_spectrum(j)
          if(ky0.lt.kymax1)then
            gamma(j) = Max(gamma0 + bz2*dvzf*kymax1  - cz1*(kymax1 - ky0)*vzf,0.0)
          elseif(USE_TTF)then
            gamma(j) = gammamax1 + bz2*dvzf*kymax1  +  Max(gamma0 - cz2*vzf*ky0,0.0)
          else
            gamma(j) = cz2*gammamax1  +  Max(gamma0 - cz2*vzf*ky0,0.0)          
          endif 
          gamma_mix(j) = gamma(j)
      enddo
    if(USE_MIX)then
      !mix over ky > kymax with integration weight = sqcky*ky0**2/(ky0**2 + cky*(ky-ky0)**2)
      do j=jmax1+2,nky
        gamma_ave = 0.0
        ky0 = ky_spectrum(j)
        mixnorm = ky0*(ATAN(sqcky*(ky_spectrum(nky)/ky0-1.0))-  &
                  ATAN(sqcky*(ky_spectrum(jmax1+1)/ky0-1.0)))
        do i=jmax1+1,nky-1
          ky1 = ky_spectrum(i)
          ky2 = ky_spectrum(i+1)
          mix1 = ky0*(ATAN(sqcky*(ky2/ky0-1.0))- ATAN(sqcky*(ky1/ky0-1.0)))
          delta = (gamma(i+1)-gamma(i))/(ky2-ky1)
          mix2 = ky0*mix1 + (ky0*ky0/(2.0*sqcky))*(LOG(cky*(ky2-ky0)**2+ky0**2)- &
                 LOG(cky*(ky1-ky0)**2+ky0**2))
          gamma_ave = gamma_ave + (gamma(i)-ky1*delta)*mix1 + delta*mix2
        enddo  
        gamma_mix(j) = gamma_ave/mixnorm  
!        write(*,*)j,ky0,gamma(j),gamma_mix(j)
      enddo  
    endif      
! intensity model
      do j=1,nky
        gamma0 = eigenvalue_spectrum_out(1,j,1)
        ky0 = ky_spectrum(j)
        kx = spectral_shift_out(j)
        do i=1,nmodes_in
          gammaeff = 0.0
          if(gamma0.gt.small)gammaeff = gamma_mix(j)*(eigenvalue_spectrum_out(1,j,i)/gamma0)**2
          if(ky0.gt.kyetg)gammaeff = gammaeff*SQRT(ky0/kyetg)
          field_spectrum_out(2,j,i) = (cnorm*gammaeff*gammaeff/ky0**4)/(1.0+ay*kx**2)**2
        enddo
     enddo
     ! recompute the intensity and flux spectra
      do j=1,nky
         do i=1,nmodes_in
            do is=1,ns
               intensity_spectrum_out(1,is,j,i) = intensity_spectrum_out(1,is,j,i)*field_spectrum_out(2,j,i)
               intensity_spectrum_out(2,is,j,i) = intensity_spectrum_out(2,is,j,i)*field_spectrum_out(2,j,i)
               do k=1,3
                  flux_spectrum_out(1,is,k,j,i) = flux_spectrum_out(1,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(2,is,k,j,i) = flux_spectrum_out(2,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(3,is,k,j,i) = flux_spectrum_out(3,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(4,is,k,j,i) = flux_spectrum_out(4,is,k,j,i)*field_spectrum_out(2,j,i)
                  flux_spectrum_out(5,is,k,j,i) = flux_spectrum_out(5,is,k,j,i)*field_spectrum_out(2,j,i)
              enddo
           enddo
        enddo
      enddo    

      END SUBROUTINE get_multiscale_spectrum
      
