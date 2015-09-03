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
      INTEGER :: i,is,k,j,j1,jmax1
      REAL :: test1,testmax1
      REAL :: gammamax1,kymax1,ky0
      REAL :: f0,f1,f2,a,b,c,x0,x02,dky,xmax
      REAL :: gamma0,gamma,gammaeff
      REAL :: cnorm, phinorm, kylow, czf, cz1,cz2
      REAL,PARAMETER :: small=1.0E-10
      !
      ! model fit parameters
      ! need to set alpha_zf_in = 1.0
      ! Miller geometry values igeo=1
      czf = alpha_zf_in
      cnorm = 13.51
      cnorm = 12.84
      cz1 = 0.35*czf
      cz2=0.60*czf
      if(igeo.eq.0)then ! s-alpha 
       cnorm=12.14
       cz1=0.61*czf
       cz2=0.60*czf
      endif
      !
      ! renormalize the fluxes and intensities to the phi-norm from the v-norm
      do j=1,nky
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
      gammamax1= eigenvalue_spectrum_out(1,1,1)
      kymax1 = ky_spectrum(1)
      testmax1 = gammamax1/kymax1
      jmax1=1
      kylow=0.8/SQRT(taus_in(2)/mass_in(2))
      j1=0
      do j=2,nky
         ky0 = ky_spectrum(j)
         if(ky0 .lt. kylow)then
           j1=j1+1
           test1 = eigenvalue_spectrum_out(1,j,1)/ky0
           if(test1 .gt. testmax1)then
            testmax1 = test1
            jmax1=j
           endif        
         endif        
      enddo
      gammamax1 = eigenvalue_spectrum_out(1,jmax1,1)
      kymax1 = ky_spectrum(jmax1)
      !interpolate to find a more accurate low-k maximum gamma/ky 
      ! this is cut of at j1 since a maximum may not exist in the low-k range
      if(jmax1.gt.1.and.jmax1.lt.j1)then
         f0 =  eigenvalue_spectrum_out(1,jmax1-1,1)/ky_spectrum(jmax1-1)
         f1 =  eigenvalue_spectrum_out(1,jmax1,1)/ky_spectrum(jmax1)
         f2 =  eigenvalue_spectrum_out(1,jmax1+1,1)/ky_spectrum(jmax1+1)
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
      do j=1,nky
! efective growthrate model
          gamma=0.0
          gamma0 = eigenvalue_spectrum_out(1,j,1)
          ky0=ky_spectrum(j)
          if(ky0.lt.kymax1)then
            gamma = Max(gamma0 - cz1*(1.0 - ky0/kymax1)*gammamax1,0.0)
          else
            gamma = cz2*gammamax1 +  Max(gamma0 - cz2*gammamax1*ky0/kymax1,0.0)
          endif   
! intensity model
        do i=1,nmodes_in
          gammaeff = 0.0
          if(gamma0.gt.small)gammaeff = gamma*(eigenvalue_spectrum_out(1,j,i)/gamma0)**2
          field_spectrum_out(2,j,i) = cnorm*gammaeff*gammaeff/ky0**4
!
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
      
