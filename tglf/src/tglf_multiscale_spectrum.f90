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
      INTEGER i,is,k,j,jmax0
      REAL test,testmax,gammamax,kymax,ky0
      REAL f0,f1,f2,a,b,c,dky,xmax
      REAL gamma0,gamma1,gamma2
      REAL cnorm
      !
      cnorm = 12.66
      if(igeo.eq.0)cnorm = 8.95
      ! renormalize the fluxes and intensities to the phi-norm from the v-norm
      do j=1,nky
         do i=1,nmodes_in
            if(field_spectrum_out(2,j,i).ne.0.0)then
               do is=1,ns
                 intensity_spectrum_out(1,is,j,i) = intensity_spectrum_out(1,is,j,i)/field_spectrum_out(2,j,i)
                 intensity_spectrum_out(2,is,j,i) = intensity_spectrum_out(2,is,j,i)/field_spectrum_out(2,j,i)
                 do k=1,3
                    flux_spectrum_out(1,is,k,j,i) = flux_spectrum_out(1,is,k,j,i)/field_spectrum_out(2,j,i)
                    flux_spectrum_out(2,is,k,j,i) = flux_spectrum_out(2,is,k,j,i)/field_spectrum_out(2,j,i)
                    flux_spectrum_out(3,is,k,j,i) = flux_spectrum_out(3,is,k,j,i)/field_spectrum_out(2,j,i)
                    flux_spectrum_out(4,is,k,j,i) = flux_spectrum_out(4,is,k,j,i)/field_spectrum_out(2,j,i)
                    flux_spectrum_out(5,is,k,j,i) = flux_spectrum_out(5,is,k,j,i)/field_spectrum_out(2,j,i)
                enddo
             enddo
          endif
        enddo
      enddo
      ! find the maximum of gamma/ky for ky < 1
      testmax=0
      jmax0=0
      do j=1,nky
         ky0 = ky_spectrum(j)
         if(ky0.le.1.0)then
           test = eigenvalue_spectrum_out(1,j,1)/ky0
           if(test .ge. testmax)then
              testmax = test
              jmax0=j
           endif
         endif        
      enddo
      gammamax = eigenvalue_spectrum_out(1,jmax0,1)
      kymax = ky_spectrum(jmax0)
      !interpolate to find a more accurate maximum gamma/ky
      if(jmax0.gt.1.and.jmax0.lt.nky)then
         f0 =  eigenvalue_spectrum_out(1,jmax0-1,1)/ky_spectrum(jmax0-1)
         f1 =  eigenvalue_spectrum_out(1,jmax0,1)/ky_spectrum(jmax0)
         f2 =  eigenvalue_spectrum_out(1,jmax0+1,1)/ky_spectrum(jmax0+1)
         a = f0
         b = 2.0*f1 - 0.5*f2 - 1.5*f0
         c = 0.5*f0 + 0.5*f2 - f1
         dky=(ky_spectrum(jmax0+1)-ky_spectrum(jmax0-1))/2.0
         xmax = -b/(2.0*c)
         kymax = ky_spectrum(jmax0-1)+dky*xmax
         gammamax = (a-0.25*b*b/c)*kymax       
      endif
      ! compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
      ! note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
      do j=1,nky
        gamma0 = 0.0
        gamma1 = eigenvalue_spectrum_out(1,j,1)
        gamma0 = gammamax
        ky0=ky_spectrum(j)
        if(ky0.le.kymax)gamma0 = gammamax*(ky0/kymax)**2.5
        if(gamma1 .ge. ky0*gammamax/kymax)gamma0 = gammamax + gamma1 - ky0*gammamax/kymax
        do i=1,nmodes_in
          gamma2 = 0.0
          if(gamma1.ne.0.0)gamma2 = gamma0*(eigenvalue_spectrum_out(1,j,i)/gamma1)**2
          field_spectrum_out(2,j,i) = cnorm*gamma2*gamma2/ky0**4
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
      
