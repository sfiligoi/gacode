!
!--------------------------------------------------------------
!
      SUBROUTINE get_multiscale_spectrum
!
!***********************************************************************
!  questions  should be addressed to
!  Gary Staebler 858-455-3466  or email: gary.staebler@gat.com
!***********************************************************************
!
!  TGLF multiscale saturation rule: sat_rule_in = 1 option
!  April 8, 2016: G.M. Staebler, J. Candy, N. T. Howard, and C. Holland
!                  Physics of Plasmas, 23 (2016) 062518
!  June 22, 2017: Retuned after coding error to Laplacian terms in Ampere 
!                 and Poisson equations was fixed
!  TGLF multiscale saturation rule: sat_rule_in = 2 option
!  Aug. 6, 2020 based on papers submitted to PPCF & NF
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
      LOGICAL :: first_pass = .TRUE.
      LOGICAL :: USE_SUB1=.FALSE.
      INTEGER :: i,is,k,j,j1,j2,jmax1,jmax2
      
      INTEGER :: expsub=2, exp_ax, jmax_mix
      REAL :: test,testmax1,testmax2
      REAL :: gammamax1,kymax1,gammamax2,kymax2,ky0,ky1,ky2
!      REAL :: f0,f1,f2,a,b,c,x0,x02,dky,xmax
      REAL :: gamma0,gamma1,gammaeff,dgamma
      REAL :: cnorm, phinorm, czf, cz1, cz2, kyetg
      REAL :: kycut
      REAL :: cky,sqcky,delta,ax,ay,kx
      REAL :: mix1,mix2,mixnorm,gamma_ave
      REAL :: vzf,dvzf,vzf1,vzf2
      REAL :: vzf_mix, kymax_mix
      REAL :: etg_streamer
      REAL :: kxzf, sat_geo_factor
      REAL :: b0,b1,b2,b3
      REAL :: d1,d2,Gq,kx_width,dlnpdr,ptot
      REAL :: measure
      REAL,DIMENSION(nkym) :: gamma_net=0.0
      REAL,DIMENSION(nkym) :: gamma=0.0
      REAL,DIMENSION(nkym) :: gamma_kymix=0.0
      REAL,PARAMETER :: small=1.0E-10
      !
       
      if(jmax_out.eq.0)then
         first_pass = .TRUE.
      else
        first_pass = .FALSE.
      endif
      if(first_pass)then  ! first pass for spectral shift model or only pass for quench rule
        gamma_net(:) = eigenvalue_spectrum_out(1,:,1)
        CALL get_zonal_mixing(nky,ky_spectrum,gamma_net,vzf_mix,kymax_mix,jmax_mix)
        vzf_out = vzf_mix
        kymax_out = kymax_mix
        jmax_out = jmax_mix
        gamma_net(:) = eigenvalue_spectrum_out(1,:,1)
!        write(*,*)"FIRST PASS: vzf_out = ",vzf_out,"  jmax_out = ",jmax_out
!        write(*,*)"FIRST PASS: kymax_out = ",kymax_out,"  gammamax_out = ",kymax_out*vzf_out
      endif
      ! model fit parameters
      ! need to set alpha_zf_in = 1.0
      ! Miller geometry values igeo=1
      if(rlnp_cutoff_in.gt.0.0)then
         dlnpdr = 0.0
         ptot = 0.0
 !        do is=1,nstotal_in   ! include all species even non-kinetic ones like fast ions
         do is=1,ns
           ptot = ptot + as(is)*taus(is)
           dlnpdr = dlnpdr + as(is)*taus(is)*(rlns(is)+rlts(is))
         enddo
         dlnpdr = rmaj_input*dlnpdr/MAX(ptot,0.01)
         if(dlnpdr .ge. rlnp_cutoff_in)dlnpdr = rlnp_cutoff_in
         if(dlnpdr .lt. 4.0)dlnpdr = 4.0
      else
         dlnpdr = 12.0
      endif
!         write(*,*)"dlnpdr = ",dlnpdr
!
      czf = ABS(alpha_zf_in)
!
! coefficients for SAT_RULE = 1
      if(sat_rule_in.eq.1) then
         kyetg=1.28
         cnorm=14.21
        if(USE_SUB1)then
          cnorm=12.12
          expsub=1
        endif
        cz1=0.48*czf
        cz2=1.0*czf
        cnorm = 14.29
        measure = SQRT(taus(1)*mass(2))
!
!        if(USE_X3)then
!         cnorm = 12.94  ! note this is normed to GASTD CGYRO units
!         cz1 = 0.0
!         cz2=1.4*czf
!         etg_streamer = 1.0
!         kyetg=etg_streamer*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
!         cky=3.0
!         sqcky=SQRT(cky)
!        endif
        if(USE_MIX)then
!original        kyetg=1.9*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
          cky=3.0
          sqcky=SQRT(cky)
          if(USE_SUB1)cnorm=12.12
          cz1=0.48*czf
!original         cz2=0.92*czf
! retuned June 22,2017
          cz2 = 1.0*czf
          etg_streamer=1.05
          if(alpha_quench_in .ne. 0.0)etg_streamer=2.1
!          kyetg=etg_streamer*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
          kyetg=etg_streamer/rho_ion  ! fixed to ion gyroradius
        endif
        if(igeo.eq.0)then ! s-alpha
           cnorm=14.63
           cz1=0.90*czf
           cz2=1.0*czf
        endif
     endif
! coefficents for SAT_RULE = 2
     if(sat_rule_in.eq.2)then
       ! SAT2 fit for CGYRO linear modes NF 2021 paper
!        b0 = 0.72
        b0 = 0.76
        b1 = 1.22
!        b2 = 24.52  ! note this is b2**2 in PPCF paper 2020
!        b2 = 11.21
!        b2 = 8.44
!        b2 = 2.13
        b2 = 3.74
        if(nmodes_in.gt.1)b2 = 3.55
!        b3 = 0.88
!        b3 = 2.4
        b3 = 1.0
!        d1 = 0.5475*((Bt0_out/B_geo0_out)**4)/grad_r0_out**2
        d1 = (Bt0_out/B_geo0_out)**4    ! PPCF paper 2020
        d1 = d1/grad_r0_out
        Gq = B_geo0_out/grad_r0_out
        d2 = b3/Gq**2
        cnorm = b2*(12.0/dlnpdr)
        kyetg = 1000.0   ! does not impact SAT2
        cky=3.0
        sqcky=SQRT(cky)
        kycut = b0*kymax_out
        cz1=0.0
        cz2 = 1.05*czf
        measure = 1.0/kymax_out     ! changed from alpha_i on Aug 8, 2021
      endif
!      write(*,*)"Bt0 = ",Bt0_out," B0 = ",B_geo0_out," gradr0 = ",grad_r0_out
!      write(*,*)"G1 = ",sat_geo1_out,"  G2 = ",sat_geo2_out,"  sat_geo0 = ",sat_geo0_out
!      write(*,*)"d1 = ",d1,"  d2 = ",d2
! coefficients for spectral shift model for ExB shear
      ax=0.0
      ay=0.0
      exp_ax = 1
      if(alpha_quench_in.eq.0.0)then
      !spectral shift model parameters
        ax = 1.15
        ay = 0.56
        exp_ax = 4
       if(sat_rule_in.eq.2)then
         ax = 1.21
         ay = 1.0
         exp_ax = 2
       endif
      endif
!    write(*,*)" ax= ",ax," ay= ",ay," kycut = ",kycut
!  apply the spectral shift temporal suppression factor
!
     if(first_pass .eqv. .FALSE.)then  ! second pass for spectral shift model
        do i=1,nky
          kx = spectral_shift_out(i)
          if(sat_rule_in.eq.2)then
            ky0=ky_spectrum(i)
            if(ky0.lt.kycut)then
              kx_width = kycut/grad_r0_out
            else
              kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
            endif
            kx = kx*ky0/kx_width
          endif
          gamma_net(i) = eigenvalue_spectrum_out(1,i,1)/(1.0 + ABS(ax*kx)**exp_ax)
!         write(*,*)ky0,kx,(ax*kx)**4
!          write(*,*)i,"gamma_net = ",gamma_net(i)
        enddo
!     write(*,*)"gammamax_out = ",vzf_out*kymax_out, gamma_net(jmax_out),"  vexb_shear = ",vexb_shear_s!!
! find the maximum of gamma_net/ky
        if(sat_rule_in.eq.1)then
          CALL get_zonal_mixing(nky,ky_spectrum,gamma_net,vzf_mix,kymax_mix,jmax_mix)
          vzf_out = vzf_mix
          kymax_out = kymax_mix
          jmax_out = jmax_mix
        else
          vzf_out = vzf_out*gamma_net(jmax_out)/MAX(eigenvalue_spectrum_out(1,jmax_out,1),small)
        endif
!        write(*,*)"2nd PASS: vzf_out = ",vzf_out,"  jmax_out = ",jmax_out
!        write(*,*)"2nd PASS: kymax_out = ",kymax_out,"  gammamax_out = ",kymax_out*vzf_out
!        write(*,*)"2nd PASS: gamma_net(jmax) = ",gamma_net(jmax_out)
      endif   ! second pass complete


! compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
      ! note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
      gammamax1= vzf_out*kymax_out
      kymax1 = kymax_out
      jmax1 = jmax_out
      vzf1 = vzf_out
!      vzf2 = gammamax2/kymax2
!      dvzf = MAX(vzf2-vzf1,0.0)
!      dgamma= dvzf*kymax1
!      gamma1=gammamax1
!      write(*,*)"dvzf = ",dvzf," vzf1 = ",vzf1," vzf2 = ",vzf2
!      write(*,*)"gammamax1 = ",gammamax1
!      write(*,*)"kymax1 = ",kymax1,"  kycut = ",kycut
      do j=1,nky
! include zonal flow effects on growthrate model
!          gamma=0.0
          gamma0 = gamma_net(j)
          ky0=ky_spectrum(j)
          if(sat_rule_in.eq.1)then
            if(ky0.lt.kymax1)then
              gamma(j) = MAX(gamma0  - cz1*(kymax1 - ky0)*vzf1,0.0)
            else
              gamma(j) = cz2*gammamax1  +  Max(gamma0 - cz2*vzf1*ky0,0.0)          
            endif
          elseif(sat_rule_in.eq.2)then
! note that the gamma mixing is for gamma_model not G*gamma_model
            if(ky0.lt.kymax1)then
              gamma(j) = gamma0
            else
              gamma(j) = gammamax1 +Max(gamma0 - cz2*vzf1*ky0,0.0)
            endif
          endif
          gamma_kymix(j) = gamma(j)     ! used if USE_MIX=F 
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
        gamma_kymix(j) = gamma_ave/mixnorm
!        write(*,*)j,ky0,gamma(j),gamma_kymix(j)
      enddo  
    endif      
! intensity model
       do j=1,nky
        gamma0 = eigenvalue_spectrum_out(1,j,1)
        ky0 = ky_spectrum(j)
        kx = spectral_shift_out(j)
        if(sat_rule_in.eq.1)then
          sat_geo_factor = SAT_geo0_out
          kx_width = ky0
        endif
        if(sat_rule_in.eq.2)then
          if(ky0.lt. kycut)then
            kx_width = kycut/grad_r0_out
            sat_geo_factor = SAT_geo0_out*d1*SAT_geo1_out
           else
             kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
             sat_geo_factor = SAT_geo0_out*(d1*SAT_geo1_out*kycut +  &
                              (ky0 - kycut)*d2*SAT_geo2_out)/ky0
           endif
           kx = kx*ky0/kx_width
    !       write(*,*)"sat2 ",B_geo0_out,grad_r0_out,SAT_geo1_out,SAT_geo2_out,sat_geo_factor
        endif
        do i=1,nmodes_in
          gammaeff = 0.0
          if(gamma0.gt.small)gammaeff = &
               gamma_kymix(j)*(eigenvalue_spectrum_out(1,j,i)/gamma0)**expsub
          if(ky0.gt.kyetg)gammaeff = gammaeff*SQRT(ky0/kyetg)
          field_spectrum_out(2,j,i) = measure*cnorm*((gammaeff/(kx_width*ky0))/(1.0+ay*kx**2))**2
          if(units_in.ne.'GYRO')field_spectrum_out(2,j,i) = sat_geo_factor*field_spectrum_out(2,j,i)
        enddo
     enddo
     ! recompute the intensity and flux spectra
      do j=1,nky
         do i=1,nmodes_in
            phinorm=field_spectrum_out(2,j,i) 
            field_spectrum_out(1,j,i) = phinorm
            field_spectrum_out(3,j,i) = QL_field_spectrum_out(3,j,i)*phinorm
            field_spectrum_out(4,j,i) = QL_field_spectrum_out(4,j,i)*phinorm
            do is=1,ns
               intensity_spectrum_out(1,is,j,i) = QL_intensity_spectrum_out(1,is,j,i)*phinorm
               intensity_spectrum_out(2,is,j,i) = QL_intensity_spectrum_out(2,is,j,i)*phinorm
               intensity_spectrum_out(3,is,j,i) = QL_intensity_spectrum_out(3,is,j,i)*phinorm
               intensity_spectrum_out(4,is,j,i) = QL_intensity_spectrum_out(4,is,j,i)*phinorm
               do k=1,3
                  flux_spectrum_out(1,is,k,j,i) = QL_flux_spectrum_out(1,is,k,j,i)*phinorm
                  flux_spectrum_out(2,is,k,j,i) = QL_flux_spectrum_out(2,is,k,j,i)*phinorm
                  flux_spectrum_out(3,is,k,j,i) = QL_flux_spectrum_out(3,is,k,j,i)*phinorm
                  flux_spectrum_out(4,is,k,j,i) = QL_flux_spectrum_out(4,is,k,j,i)*phinorm
                  flux_spectrum_out(5,is,k,j,i) = QL_flux_spectrum_out(5,is,k,j,i)*phinorm
              enddo
           enddo
        enddo
      enddo    

      END SUBROUTINE get_multiscale_spectrum
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE get_zonal_mixing(nmix,ky_mix,gamma_mix,vzf_mix,kymax_mix,jmax_mix)
!
!  finds the maximum of gamma/ky spectrum vzf_out and kymax_out
!
    USE tglf_dimensions
    USE tglf_global
    USE tglf_species
    IMPLICIT NONE
    LOGICAL :: use_kymin = .false.
    INTEGER :: nmix,i,k,j,j1,j2,jmax1,jmin
    REAL :: test,testmax,peakmax
    REAL :: kx, kyhigh, kycut, kymin
    REAL :: gammamax1,kymax1,testmax1,ky0,ky1,ky2
    REAL :: f0,f1,f2,a,b,c,x1,deltaky,xmax,xmin
    REAL :: vzf1, vzf2
    REAL :: kymax_mix, vzf_mix
    INTEGER :: jmax_mix, down
    REAL, DIMENSION(nkym) :: gamma_mix, ky_mix
!  initialize output of subroutine
    vzf_mix = 0.0
    kymax_mix = 0.0
    jmax_mix = 1
    xmin = 0.0
    if(alpha_zf_in.lt.0.0)use_kymin = .true. 
!
! find the local maximum of gamma_mix/ky_mix with the largest gamma_mix/ky_mix^2
!
      gammamax1= gamma_mix(1)
      kymax1 = ky_mix(1)
      testmax = 0.0
      peakmax=0.0
      jmax_mix=1
      kycut=0.8/rho_ion
      kymin = 0.0
      jmin = 0
      if(use_kymin)kymin = 0.173*sqrt(2.0)/rho_ion
      if(sat_rule_in.eq.2)then
        kycut = grad_r0_out*kycut
        kymin = grad_r0_out*kymin
      endif
!      write(*,*)" kycut = ",kycut," kymin = ",kymin
      ! find the low and high ky peaks of gamma/ky
      do j=1,nky-1
       if(ky_mix(j).lt.kymin)jmin = j
       if((ky_mix(j+1).ge.kymin).and.(ky_mix(j).le.kycut))then
         j1=j
         kymax1 = ky_mix(j)
         testmax1 = gamma_mix(j)/kymax1
!         write(*,*)"j=",j,"ky = ",ky0," gamma_net = ",gamma_net(j)
         if(testmax1.gt.testmax)then
           testmax = testmax1
           jmax_mix = j
          endif
         endif
      enddo
      if(testmax.eq.0.0)jmax_mix=j1
! no unstable modes in range set kymax index to end of range
! this is cut of at j1 since a maximum may not exist in the low-k range
      kymax1 = ky_mix(jmax_mix)
      gammamax1 = gamma_mix(jmax_mix)
      if(kymax1.lt.kymin)then
          kymax1 = kymin
!interpolate to find a more accurate low-k maximum gamma/ky
         gammamax1 = gamma_mix(1)   &
          +(gamma_mix(2)-gamma_mix(1))*(kymin-ky_mix(1))/(ky_mix(2)-ky_mix(1))
      endif
!        write(*,*)" jmax_mix = ",jmax_mix,"  gammamax1 = ",gammamax1," kymax1 = ",kymax1
      if(jmax_mix .gt. 1 .and. jmax_mix .lt. j1)then
 !        write(*,*)"refining low-k maximum"
! determine kymax1 and gammamax1 bounded by the tree points f0,f1,f2
! use a quadratic fit: f = a + b x + c x^2  to f = gamma/ky centered at jmax1
         jmax1 = jmax_mix
         f0 =  gamma_mix(jmax1-1)/ky_mix(jmax1-1)
         f1 =  gamma_mix(jmax1)/ky_mix(jmax1)
         f2 =  gamma_mix(jmax1+1)/ky_mix(jmax1+1)
         deltaky = ky_mix(jmax1+1)-ky_mix(jmax1-1)
         x1 = (ky_mix(jmax1)-ky_mix(jmax1-1))/deltaky
         a = f0
         b = (f1 - f0*(1-x1*x1)-f2*x1*x1)/(x1-x1*x1)
         c = f2 - f0 - b
!         write(*,*)"f0 = ",f0,"  f1 = ",f1,"  f2 = ",f2,"  x1 = ",x1
         if(f0 .ge.f1)then
! if f0>f1 then f1 is not a local maximum
             kymax1 = ky_mix(jmax1-1)
             gammamax1 = f0*kymax1
             if(kymax1.lt.kymin)then
!interpolate to find the value of gammamax1 at kymin
               kymax1 = kymin
               xmin = (kymin - ky_mix(jmax1-1))/deltaky
               gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
             endif
        endif
        if(f0.lt.f1  )then
!        if f0<f1 then f1>f2 due to the maximum search
! use the quadratic fit to refine the local maximum:
             xmax = -b/(2.0*c)
             xmin = 0.0
             if(ky_mix(jmax1-1).lt.kymin)then
               xmin = (kymin - ky_mix(jmax1-1))/deltaky
             endif
             if(xmax .ge. 1.0)then
! if xmax >= 1  use f2 as the maximum
               kymax1 = ky_mix(jmax1+1)
               gammamax1 = f2*kymax1
             elseif(xmax.le.xmin)then
               if(xmin .gt. 0.0)then
                 kymax1 = kymin
! use the quadratic fit to determine gammamax1 at kymin
                 gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
               elseif(xmax .le. 0.0)then
! if xmax<=0 use f0 as the maximum
                 kymax1 = ky_mix(jmax1-1)
                 gammamax1 = f0*kymax1
               endif
             else
! the conditions f0<f1<f2 and xmin<xmax<1 are satisfied
! use the quadratic fit to determine gammamax1 and kymax1
               kymax1 = ky_mix(jmax1-1)+deltaky*xmax
               gammamax1 = (a+b*xmax+c*xmax*xmax)*kymax1
             endif !xmax tests
           endif !f0 < f1
       endif  ! jmax_mix > 1
       vzf_mix = gammamax1/kymax1
       kymax_mix = kymax1
!      jmax_mix = jmax1
!      write(*,*)"get_zonal_mxing"
!      write(*,*)"xmax = ",xmax, "  xmin = ",xmin
!      write(*,*)"gammamax1 = ",gammamax1," kymax1 = ",kymax1, "  kymin = ",kymin
 
 END SUBROUTINE get_zonal_mixing
