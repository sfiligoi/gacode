!
!--------------------------------------------------------------
!
      SUBROUTINE gftm_LS
!
!***********************************************************************
! questions  should be addressed to
!  Gary Staebler 858-455-3466  or email: gary.staebler@gat.com
!***********************************************************************
!  gftm Linear Stability driver computes all quasilinear quantities
!  for a single ky
!
!***********************************************************************
!
      USE gftm_dimensions
      USE gftm_global
      USE gftm_species
      USE gftm_eigen
      USE gftm_GFS
      USE gftm_weight
!
      IMPLICIT NONE
      REAL,PARAMETER :: epsilon1 = 1.E-12
      INTEGER :: info
      INTEGER :: j1, j, i, iroot, jmax(maxmodes)
      INTEGER :: imax,is
      INTEGER :: mi,me
      INTEGER :: ip,jp
      REAL :: zgamax
      REAL :: phi2_bar,kyi
      REAL :: sum_v_bar, sum_modB_bar
      REAL :: one, henorm
      COMPLEX :: testhe
!  ZGESV storage
      REAL :: small = 1.0E-10
      REAL :: get_gamma_net
!
!      cputime0=MPI_WTIME()
!
! set ky in units of k_theta*rho_s  rho_s=C_s/omega_s
! C_s=sqrt(Te/mi), omega_s=eB/(mi c)
!
      ky = ky_s
!      write(*,*)"ky = ",ky
      !
      ! check co-dependencies
      !
      if(new_geometry)new_width=.TRUE.
      if(new_width)new_matrix=.TRUE.
!      write(*,*)"new_start=",new_start
!      write(*,*)"new_geometry=",new_geometry
!      write(*,*)"new_width=",new_width
!      write(*,*)"new_matrix=",new_matrix
!      write(*,*)"width_in=",width_in
!      write(*,*)"nbasis_max_in=",nbasis_max_in
!      write(*,*)"new_eikonal_in=",new_eikonal_in
!      write(*,*)"betae_in=",betae_in
!      write(*,*)"debye_in=",debye_in
!      write(*,*)"xnue_in=",xnue_in
!      write(*,*)"zeff_in=",zeff_in
!      write(*,*)"filter_in=",filter_in
!      write(*,*)"park_in=",park_in
!      write(*,*)"ghat_in=",ghat_in
!      write(*,*)"wd_zero_in=",wd_zero_in
!      write(*,*)"Linsker_factor_in=",Linsker_factor_in
!      write(*,*)"gradB_factor_in=",gradB_factor_in
!      write(*,*)"xnu_factor_in=",xnu_factor_in
!      write(*,*)"debye_factor_in=",debye_factor_in
!      write(*,*)"theta_trapped_in=",theta_trapped_in
!      write(*,*)"alpha_p_in=",alpha_p_in
!      write(*,*)"alpha_e_in=",alpha_e_in
!
      ky = ky_s
!
      if(new_eikonal_in)then
        if(new_geometry)then
          if(igeo.eq.1)then
! set up MILLER geometry
            trace_path(4)=1
            call miller_geo
          elseif(igeo.eq.2)then
! set up FOURIER geometry
            trace_path(5)=1
            call fourier_geo
          elseif(igeo.eq.3)then
! set up ELITE geometry
            trace_path(8)=1
            call ELITE_geo
          endif
!        write(*,*)"geometry done"
!
! compute the eikonal functions for general geometry (igeo>0)
!
          if(igeo.gt.0)call mercier_luc 
!
        endif ! new_geometry
        new_geometry=.FALSE.
!
!  load the x-grid eikonal functions wdx,b0x
!
        if(new_width)then
          trace_path(6)=1
          call get_xgrid_functions
        endif
!       write(*,*)"ft=",ft
      endif  !new_eikonal_in
!      write(*,*)"eikonal done"
!
if(new_matrix)then
        trace_path(7)=1
        call get_matrix
      endif
!      write(*,*)"matrix done"
!
! allocate eigenvalues, eigenvectors
!
!      write(*,*)"eigen allocation done"
!
!  solver for linear eigenmodes of gftm equations
!
      call gftm_eigensolver
!      write(*,*)"eigensolver done"
!
!      initalize output to zero
!
      do j1=1,maxmodes
        jmax(j1) = 0
        gamma_out(j1) = 0.0
        freq_out(j1) = 0.0
        v_QL_out(j1) = 0.0
        a_par_QL_out(j1) = 0.0
        b_par_QL_out(j1) = 0.0
        phi_bar_out(j1) = 0.0
        a_par_bar_out(j1) = 0.0
        b_par_bar_out(j1) = 0.0
        v_bar_out(j1) = 0.0
        do is=1,nsm
          particle_QL_out(j1,is) = 0.0
          energy_QL_out(j1,is) = 0.0
          stress_par_QL_out(j1,is) = 0.0
          stress_tor_QL_out(j1,is) = 0.0
          exchange_QL_out(j1,is) = 0.0
          N_QL_out(j1,is) = 0.0
          T_QL_out(j1,is) = 0.0
          U_QL_out(j1,is) = 0.0
          Q_QL_out(j1,is) = 0.0
          N_bar_out(j1,is) = 0.0
          T_bar_out(j1,is) = 0.0
          U_bar_out(j1,is) = 0.0
          Q_bar_out(j1,is) = 0.0
          Ns_Ts_phase_out(j1,is) = 0.0
        enddo
        ne_te_phase_out(j1) = 0.0
      enddo
!
      if(ibranch_in.eq.0)then
!       write(*,*)"nmodes_in=",nmodes_in
      nmodes_out = nmodes_in
!
! sort the unstable modes into electron and ion frequencies
      mi = 0
      me = 0
      do j1=1,ntot
        if(ri(j1).gt.epsilon1)then
         if(rr(j1).lt.0.0)then
           mi = mi+1
           di(mi)=j1
         else
           me = me+1
           de(me)=j1
         endif
!          write(*,*)"debug sort",j1,rr(j1),ri(j1),me,mi
        endif
      enddo
! find the most unstable mode for each branch
      if(me.gt.0)then
        zgamax = 0.0
        do iroot=1,me
         if(ri(de(iroot)).gt.zgamax)then
           zgamax = ri(de(iroot))
           jmax(1) = de(iroot)
         endif
        enddo
        gamma_out(1)=ri(jmax(1))
        freq_out(1)=rr(jmax(1))
      endif
      if(mi.gt.0)then
        zgamax = 0.0
        do iroot=1,mi
          if(ri(di(iroot)).gt.zgamax)then
           zgamax = ri(di(iroot))
           jmax(2) = di(iroot)
         endif
        enddo
        gamma_out(2)=ri(jmax(2))
        freq_out(2)=rr(jmax(2))
      endif
!        write(*,*)"debug sort", ky
!        write(*,*)"nmodes_out = ",nmodes_out,nmodes_in
!        write(*,*)"me = ",me,"mi = ",mi
!        write(*,*)"gamma_out=",gamma_out(1),gamma_out(2)
!        write(*,*)"freq_out=",freq_out(1),freq_out(2)
      endif ! ibranch_in .eq.0 
      if(ibranch_in.eq.-1)then
!
!  find the top nmodes most unstable modes
       CALL sort_eigenvalues(nmodes_in,jmax)
        nmodes_out = 0
         do j1=1,nmodes_in 
          if(jmax(j1).ne.0)then
            nmodes_out = nmodes_out + 1
            gamma_out(j1)=ri(jmax(j1))
            freq_out(j1)=rr(jmax(j1))
          endif
        enddo
!       write(*,*)"gamma_out = ",gamma_out(1)," freq_out = ", freq_out(1)
!      write(*,*)"debug jmax =",jmax(1),jmax(2)
!      write(*,*)"debug nmodes_out = ",nmodes_out
      endif ! ibranch_in .eq. -1
!
      if(alpha_quench_in.ne.0.0)then
! apply quench rule
        do j1=1,nmodes_in
          gamma_out(j1) = get_gamma_net(gamma_out(j1))
        enddo
      elseif(vexb_shear_s.ne.0.0)then
! use spectral shift model for second pass
       do j1=1,nmodes_in
          gamma_out(j1) = gamma_reference_kx0(j1)
          freq_out(j1) = freq_reference_kx0(j1)
        enddo
      endif
!      
!  get the fluxes for the most unstable modes
!
         do imax=1,nmodes_out
         if(jmax(imax).gt.0)then
!  alpha/beta=-(frequency+xi*growthrate)
           we = alpha(jmax(imax))/beta(jmax(imax))
!          write(*,*)"eigenvalue=",we,jmax(imax)
          if(iflux_in)then
           he(:) = hetot(:,jmax(imax))
          else
            do i=1,ntot
            he(i) = small
            do j=1,ntot
              emat(i,j) = we*bmat(i,j) - amat(i,j) + small*one(i,j)
            enddo
           enddo
!  find the eigenvector he correponding to we eigenvalue
          call zgesv(ntot,1,emat,ntot,ipive,he,ntot,info)

          ! Check for clean exit from ZGESV
!          write(*,*)" info = ",info
!          write(*,*)" amat(1,1) = ",amat(1,1)," bmat(1,1) = ",bmat(1,1), "zmat(1,1) = ",zmat(1,1)
          if (info /= 0)CALL gftm_error(1,"1st call to ZGESV failed in gftm_LS")
!
          endif  ! iflux_in
! debug
          go to 10
          henorm = 0.0
          do i=1,ntot
            henorm = henorm + REAL(CONJG(he(i))*he(i))
          enddo
          henorm = SQRT(henorm)
          write(*,*)" henorm = ",henorm
          do i=1,ntot
             he(i) = he(i)/henorm
             write(*,*)i,"  he = ",he(i)
           enddo
           do i=1,ntot
           testhe = 0.0
           do j=1,ntot
              testhe = testhe + (we*bmat(i,j) - amat(i,j))*he(j)
           enddo
            write(*,*)i,"  test he = ", testhe
           enddo
     10   continue
!
!  compute propagator gmatinv for the eigenvalue we
!  separate in to species block diagonal form before inverting
!
         
          do is=1,ns
            do ip=1,nphase
            do jp=1,nphase
              gmatinv(is,ip,jp) = 0.0
              i = ip + nphase*(is-1)
              j = jp + nphase*(is-1)
               gmat(ip,jp) = one(ip,jp)
               zmat(ip,jp) = we*one(i,j) - mateq(i,j)
!            write(*,*)i,j,"  zmat = ",zmat(i,j)
            enddo
            enddo
            call zgesv(nphase,nphase,zmat,nphase,ipiv,gmat,nphase,info)
! Check for clean exit from ZGESV
            if (info /= 0)CALL gftm_error(1,"2nd call to ZGESV failed in gftm_LS")
!
            do ip=1,nphase
            do jp=1,nphase
              gmatinv(is,ip,jp) = gmat(ip,jp)
            enddo
            enddo
          enddo  ! is
! debug
!          do i=1,ntot
!          do j=1,ntot
!            write(*,*)i,j,"  gmatinv = ",gmatinv(i,j)
!          enddo
!          enddo
!
!  compute the response matrix and QL weights for this eigenvalue
!
          call get_QL_weights
!
!  save the results for output
!
          wd_bar_out(imax)=wd_bar
          b0_bar_out(imax)=b0_bar
          modB_bar_out(imax)=modB_bar
          v_QL_out(imax)=v_weight
          a_par_QL_out(imax)=a_par_weight
          b_par_QL_out(imax)=b_par_weight
          kx_bar_out(imax)=kx_bar
          kpar_bar_out(imax)=kpar_bar/(R_unit*q_unit*width_in)
          do i=1,nbasis
            do j=1,3
             field_weight_out(imax,j,i)=field_weight_QL_out(j,i)
            enddo
          enddo
!
          do is=ns0,ns
            particle_QL_out(imax,is)=flux_out(is,1)
            energy_QL_out(imax,is)=taus(is)*(flux_out(is,2)/sqrt_two+flux_out(is,3)+1.5*flux_out(is,1))
!            stress_par_QL_out(imax,is)=stress_par_weight(is)
            stress_tor_QL_out(imax,is) = mass(is)*taus(is)*flux_out(is,4)
            exchange_QL_out(imax,is) = delta_out(is)
            N_QL_out(imax,is)=N_weight(is)
            T_QL_out(imax,is)=T_weight(is)
            U_QL_out(imax,is)=U_weight(is)
            Q_QL_out(imax,is)=Q_weight(is)
            Ns_Ts_phase_out(imax,is)=Ns_Ts_phase(is)
          enddo
          ne_te_phase_out(imax) = Ne_Te_phase
          phi_bar_out(imax) = phi2_bar
          a_par_bar_out(imax) = phi2_bar*a_par_QL_out(imax)
          b_par_bar_out(imax) = phi2_bar*b_par_QL_out(imax)
          do is=ns0,ns
            N_bar_out(imax,is)=phi2_bar*N_QL_out(imax,is)
            T_bar_out(imax,is)=phi2_bar*T_QL_out(imax,is)
            U_bar_out(imax,is)=phi2_bar*U_QL_out(imax,is)
            Q_bar_out(imax,is)=phi2_bar*Q_QL_out(imax,is)
          enddo
         endif
        enddo !imax
!
!
      END SUBROUTINE gftm_LS
!
!--------------------------------------------------------------
!
      REAL FUNCTION get_gamma_net(gp)
!
      USE gftm_global
      USE gftm_species
!
      REAL,INTENT(IN) :: gp
      REAL :: alpha_exb
!
      alpha_exb = 0.3
      if(igeo.eq.1)alpha_exb=0.3*SQRT(kappa_loc)
      get_gamma_net =  MAX(gp - ABS(alpha_exb*alpha_quench_in*vexb_shear_s),0.0)
!
      END FUNCTION get_gamma_net
!
!--------------------------------------------------------------
!
      REAL FUNCTION get_GAM_freq()
!
      USE gftm_global
      USE gftm_species
!
      REAL :: elong=1.0
!
      if(igeo.eq.1)elong = kappa_loc
!      get_GAM_freq = SQRT(3.5*taus(2)+2.0*taus(1))/R_unit
      get_GAM_freq = (2.0/(1.0+elong))*SQRT(taus(2)+taus(1))/R_unit
!
      END FUNCTION get_GAM_freq
!
!--------------------------------------------------------------
!
!
      SUBROUTINE sort_eigenvalues(nsorted,sorted_index)
!
      USE gftm_dimensions
      USE gftm_eigen
!
      IMPLICIT NONE
      INTEGER :: i, j, ngrow,imid, iend, it, nsorted
      INTEGER :: sorted_index(maxmodes)
      REAL,PARAMETER ::  epsilon1=1.0E-12
!
!  find all of the unstable modes and the most unstable
!
      ngrow = 0
!      write(*,*)ntot,"eigenvalues"
      do i=1,ntot
!        growthrate(i)=0.0
!        frequency(i)=0.0
        grow_index(i)=0
        if(ri(i).gt.epsilon1)then
          ngrow = ngrow + 1
          grow_index(ngrow) = i
!          write(*,*)
!     >    ngrow,grow_index(ngrow),ri(grow_index(ngrow))
        endif 
      enddo
!
      if(ngrow.le.1)go to 30
!
! put the unstable modes in ascending order by growthrate
! using a heapsort algorithm 
!
      imid = ngrow/2 + 1
      iend = ngrow
!      write(*,*)"imid = ",imid,"iend = ",iend
 10   continue
        if(imid.gt.1)then
          imid = imid-1
          it = grow_index(imid)
        else
          it = grow_index(iend)
          grow_index(iend) = grow_index(1)
          iend = iend-1
          if(iend.eq.1)then
            grow_index(1)=it
          do i=1,ngrow
!           write(*,*)i,grow_index(i),ri(grow_index(i))
          enddo
            go to 30
         endif
        endif
        i=imid
        j = imid+imid
!        write(*,*)"i = ",i,"j = ",j
 20     if(j.le.iend)then
          if(j.lt.iend)then
           if(ri(grow_index(j)).gt.ri(grow_index(j+1)))j=j+1
          endif
          if(ri(it).gt.ri(grow_index(j)))then
            grow_index(i)=grow_index(j)
            i=j
            j=j+j
          else
            j = iend+1
          endif
        go to 20
        endif
        grow_index(i)=it
      go to 10
 30   continue
      do i=1,nsorted
        sorted_index(i)=0
        if(i.le.ngrow)sorted_index(i)=grow_index(i)
      enddo
! debug
!      do i=1,ngrow
!        growthrate(i)=ri(grow_index(i))
!        frequency(i)=rr(grow_index(i))
!      enddo
!      do i=1,ngrow
!      write(*,*)i,growthrate(i),frequency(i),grow_index(i)
!      enddo
!
      END SUBROUTINE sort_eigenvalues
!
      SUBROUTINE get_QL_weights
! **************************************************************
!
! compute the quasilinear weights for a single eigenmode
! with eigenvector v. All of the QL weights are normalized to phi_norm
!
!
!***************************************************************
      USE gftm_dimensions
      USE gftm_global
      USE gftm_species
!      USE gftm_hermite
      USE gftm_eigen
      USE gftm_coeff
      USE gftm_weight
      USE gftm_sgrid
      USE gftm_gyro_average
      USE gftm_GFS
! 
      IMPLICIT NONE
!
      INTEGER :: ib,jb,is,iu,ie,k1,k2,kp
      INTEGER :: i,j,k, ip,jp
!
      COMPLEX :: n(nsm,nb)
      COMPLEX :: u_par(nsm,nb)
      COMPLEX :: t_par(nsm,nb)
      COMPLEX :: t_per(nsm,nb)
      COMPLEX :: q_par(nsm,nb)
      COMPLEX :: q_per(nsm,nb)
      COMPLEX :: temp(nsm,nb)
      COMPLEX :: stress_par(nsm,nb,3),stress_per(nsm,nb,3)
      COMPLEX :: phi_wd_phi,wd_phi
      COMPLEX :: phi_b0_phi,b0_phi
      COMPLEX :: phi_modB_phi,modB_phi
      COMPLEX :: phi_kx_phi,kx_phi
      COMPLEX :: phi_kpar_phi,kpar_phi
      REAL :: betae_psi,betae_sig
      REAL :: phi_norm,psi_norm,sig_norm,vnorm
      REAL :: epsilon1
      REAL :: stress_correction,wp
      REAL :: one
!
!      xi=(0.0,1.0)
      epsilon1 = 1.E-12
!
!  compute the electromagnetic potentials
!
      betae_psi = 0.0
      if(use_bper_in)betae_psi = betae_s
      betae_sig = 0.0
      if(use_bpar_in)betae_sig = betae_s
      do ib=1,nbasis
        phi(ib)=0.0
        psi(ib)=0.0
        sig(ib)=0.0
        do k1=1,ntot
          phi(ib) = phi(ib) + phib(ib,k1)*he(k1)
          psi(ib) = psi(ib) + psib(ib,k1)*he(k1)
          sig(ib) = sig(ib) + sigb(ib,k1)*he(k1)
        enddo
      enddo
!
!  compute phi_norm, psi_norm, sig_norm
!
      phi_norm = 0.0
      psi_norm = 0.0
      sig_norm = 0.0
      do ib=1,nbasis
        phi_norm = phi_norm + REAL(phi(ib)*CONJG(phi(ib)))
        psi_norm = psi_norm + REAL(psi(ib)*CONJG(psi(ib)))
        sig_norm = sig_norm + REAL(sig(ib)*CONJG(sig(ib)))
      enddo
      if(phi_norm.lt.epsilon1)phi_norm = epsilon1
!      write(*,*)"phi_norm =",phi_norm
!       write(*,*)" phinorm = ",phi_norm
!       write(*,*)" psinorm = ",psi_norm
!       write(*,*)" signorm = ",sig_norm
!
! save the field weights
      do i=1,nbasis
        field_weight_QL_out(1,i) = xi*phi(i)/SQRT(phi_norm)
        field_weight_QL_out(2,i) = xi*psi(i)/SQRT(phi_norm)
        field_weight_QL_out(3,i) = xi*sig(i)/SQRT(phi_norm)
!        write(*,*)i,"field_weight=",field_weight_QL_out(1,i),field_weight_QL_out(2,i),field_weight_QL_out(3,i)
       enddo
!
! compute <phi|*|phi> averages
      phi_wd_phi = 0.0
      phi_b0_phi = 0.0
      phi_modB_phi = 0.0
      phi_kx_phi = 0.0
      phi_kpar_phi = 0.0
      do i=1,nbasis
         wd_phi = 0.0
         b0_phi = 0.0
         modB_phi = 0.0
         kx_phi = 0.0
         kpar_phi = 0.0
         do j=1,nbasis
           wd_phi = wd_phi +ave_wdpar(i,j)*phi(j)
           b0_phi = b0_phi +ave_b0(i,j)*phi(j)
           modB_phi = modB_phi +ave_c_par_par(i,j)*phi(j)
           kx_phi = kx_phi +ave_kx(i,j)*phi(j)
           kpar_phi = kpar_phi +xi*ave_kpar(i,j)*phi(j)
         enddo
         phi_wd_phi = phi_wd_phi + CONJG(phi(i))*wd_phi
         phi_b0_phi = phi_b0_phi + CONJG(phi(i))*b0_phi
         phi_modB_phi = phi_modB_phi + CONJG(phi(i))*modB_phi
         phi_kx_phi = phi_kx_phi + CONJG(phi(i))*kx_phi
         phi_kpar_phi = phi_kpar_phi + CONJG(phi(i))*kpar_phi
      enddo
      wd_bar = REAL(phi_wd_phi)/phi_norm
      b0_bar = REAL(phi_b0_phi)/phi_norm
      modB_bar = ABS(REAL(phi_modB_phi)/phi_norm)
      kx_bar = REAL(phi_kx_phi)/phi_norm
      kpar_bar = REAL(phi_kpar_phi)/phi_norm
!      write(*,*)"wd_bar = ",wd_bar
!      write(*,*)"b0_bar = ",b0_bar
!       write(*,*)"modB_bar = ",modB_bar
!      write(*,*)"kx_bar = ",kx_bar
!      write(*,*)"kpar_bar = ",kpar_bar
!
! fill the moments
!
      do is=1,ns
      do ie=1,ne
      do ib=1,nbasis
        pe1j0phi(is,ie,ib)=0.0
        pe1j0psi(is,ie,ib)=0.0
        pe1j1sig(is,ie,ib)=0.0
        pe2j0phi(is,ie,ib)=0.0
        pe2j0psi(is,ie,ib)=0.0
        pe2j1sig(is,ie,ib)=0.0
        do jb=1,nbasis
          pe1j0phi(is,ie,ib) = pe1j0phi(is,ie,ib) + mat_pe1j0s(is,ie,ib,jb)*phi(jb)
          pe1j0psi(is,ie,ib) = pe1j0psi(is,ie,ib) + mat_pe1j0s(is,ie,ib,jb)*psi(jb)
          pe1j1sig(is,ie,ib) = pe1j1sig(is,ie,ib) + mat_pe1j1s(is,ie,ib,jb)*sig(jb)
          pe2j0phi(is,ie,ib) = pe2j0phi(is,ie,ib) + mat_pe2j0s(is,ie,ib,jb)*phi(jb)
          pe2j0psi(is,ie,ib) = pe2j0psi(is,ie,ib) + mat_pe2j0s(is,ie,ib,jb)*psi(jb)
          pe2j1sig(is,ie,ib) = pe2j1sig(is,ie,ib) + mat_pe2j1s(is,ie,ib,jb)*sig(jb)
        enddo
      enddo
      enddo
      enddo
!
! now summ the gyro-average Maxwell fields
!
     do is = 1,ns
     do ie = 1,ne
     do iu = 1,nu
     do ib = 1,nbasis
         kp = ib+nbasis*(iu-1)+nbasis*nu*(ie-1)
         k1 = kp + nphase*(is-1)
         hes(is,kp) = he(k1)
         pwns(is,kp) = one(1,iu)*pe1j0phi(is,ie,ib)                             &
               - vs(is)*one(2,iu)*pe1j0psi(is,ie,ib)                            &
               + 2.0*(taus(is)/zs(is))*one(1,iu)*pe1j1sig(is,ie,ib)
         pwtpars(is,kp) = one(3,iu)*pe1j0phi(is,ie,ib)                          &
               - vs(is)*sqrt_two*matu(3,iu)*pe1j0psi(is,ie,ib)                  &
               + 2.0*(taus(is)/zs(is))*one(3,iu)*pe1j1sig(is,ie,ib)
         pwtpers(is,kp) = one(1,iu)*pe2j0phi(is,ie,ib)                          &
               - vs(is)*one(2,iu)*pe2j0psi(is,ie,ib)                            &
               + 2.0*(taus(is)/zs(is))*one(1,iu)*pe2j1sig(is,ie,ib)
         pwps(is,kp) = one(2,iu)*pe1j0phi(is,ie,ib)                             &
               - vs(is)*sqrt_two*matu(2,iu)*pe1j0psi(is,ie,ib)                  &
               + 2.0*(taus(is)/zs(is))*one(2,iu)*pe1j1sig(is,ie,ib)
     enddo
     enddo
     enddo
     enddo
     do is=1,ns
     do ip=1,nphase
       hwns(is,ip) = 0.0
       hwtpars(is,ip) = 0.0
       hwtpers(is,ip) = 0.0
       hwes(is,ip) = 0.0
       hwps(is,ip) = 0.0
       do jp=1,nphase
         hwns(is,ip) = hwns(is,ip) + gmatinv(is,ip,jp)*pwns(is,jp)
         hwtpars(is,ip) = hwtpars(is,ip) + gmatinv(is,ip,jp)*pwtpars(is,jp)
         hwtpers(is,ip) = hwtpers(is,ip) + gmatinv(is,ip,jp)*pwtpers(is,jp)
         hwps(is,ip) = hwps(is,ip) + gmatinv(is,ip,jp)*pwps(is,jp)
       enddo
       hwes(is,ip) = hwns(is,ip)
    enddo
    enddo
!
! conjugates of field perterbation capital PSI for each species
!
      do is=1,ns
        do kp = 1, nphase
          PSIpu1pe1(is,kp) = CONJG(pwns(is,kp))
          PSIpu2pe1(is,kp) = CONJG(pwps(is,kp))
          PSIpu3pe1(is,kp) = CONJG(pwtpars(is,kp))
          PSIpu1pe2(is,kp) = CONJG(pwtpers(is,kp))
!          write(*,*)is,kp,"  phipu2pe1j0 = ",PHIpu2pe1j0(is,kp)
        enddo
      enddo
!
      do is = 1,ns
        do j=1,4
          rwn(is,j) = 0.0
          rwp(is,j) = 0.0
          rwtpar(is,j) = 0.0
          rwtper(is,j) = 0.0
          fluxe(is,j) = 0.0
        enddo
        do kp=1,nphase
        rwn(is,1) = rwn(is,1) + PSIpu1pe1(is,kp)*hwns(is,kp)
        rwtpar(is,1) = rwtpar(is,1) + PSIpu1pe1(is,kp)*hwtpars(is,kp)
        rwtper(is,1) = rwtper(is,1) + PSIpu1pe1(is,kp)*hwtpers(is,kp)
        rwp(is,1) = rwp(is,1) + PSIpu1pe1(is,kp)*hwps(is,kp)
        fluxe(is,1) = fluxe(is,1) + PSIpu1pe1(is,kp)*hes(is,kp)
        rwn(is,2) = rwn(is,2) + PSIpu3pe1(is,kp)*hwns(is,kp)
        rwtpar(is,2) = rwtpar(is,2) + PSIpu3pe1(is,kp)*hwtpars(is,kp)
        rwtper(is,2) = rwtper(is,2) + PSIpu3pe1(is,kp)*hwtpers(is,kp)
        rwp(is,2) = rwp(is,2) + PSIpu3pe1(is,kp)*hwps(is,kp)
        fluxe(is,2) = fluxe(is,2) + PSIpu3pe1(is,kp)*hes(is,kp)
        rwn(is,3) = rwn(is,3) + PSIpu1pe2(is,kp)*hwns(is,kp)
        rwtpar(is,3) = rwtpar(is,3) + PSIpu1pe2(is,kp)*hwtpars(is,kp)
        rwtper(is,3) = rwtper(is,3) + PSIpu1pe2(is,kp)*hwtpers(is,kp)
        rwp(is,3) = rwp(is,3) + PSIpu1pe2(is,kp)*hwps(is,kp)
        fluxe(is,3) = fluxe(is,3) + PSIpu1pe2(is,kp)*hes(is,kp)
        rwn(is,4) = rwn(is,4) + PSIpu2pe1(is,kp)*hwns(is,kp)
        rwtpar(is,4) = rwtpar(is,4) + PSIpu2pe1(is,kp)*hwtpars(is,kp)
        rwtper(is,4) = rwtper(is,4) + PSIpu2pe1(is,kp)*hwtpers(is,kp)
        rwp(is,4) = rwp(is,4) + PSIpu2pe1(is,kp)*hwps(is,kp)
        fluxe(is,4) = fluxe(is,4) + PSIpu2pe1(is,kp)*hes(is,kp)
        enddo
      enddo
!
      do is = 1,ns  ! species index
      do j = 1,4    ! gradient drive index
        diff_out(is,j,1) = as(is)*REAL(xi*ky*rwn(is,j))/phi_norm
        diff_out(is,j,2) = as(is)*REAL(xi*ky*rwtpar(is,j))/phi_norm
        diff_out(is,j,3) = as(is)*REAL(xi*ky*rwtper(is,j))/phi_norm
        diff_out(is,j,4) = as(is)*REAL(xi*ky*rwp(is,j))/phi_norm
        conv_out(is,j) = as(is)*REAL(xi*ky*we*rwn(is,j))/phi_norm
        flux_out(is,j) = as(is)*REAL(xi*ky*fluxe(is,j)/phi_norm)
      enddo
        delta_out(is) = -zs(is)*as(is)*REAL(we)*REAL(xi*fluxe(is,1))
      enddo
!
! debug
!      write(*,*)"fluxe(1) = ",fluxe(1,1),fluxe(2,1)
!      do is=1,ns
!        write(*,*)is,"diff"
!      do i=1,4
!        write(*,*)(diff_out(is,i,j),j=1,4)
!      enddo
!        write(*,*)is,"conv"
!        write(*,*)(conv_out(is,j),j=1,4)
!        write(*,*)is,"flux"
!        write(*,*)(flux_out(is,j),j=1,4)
!      enddo
! check summ of coefficients
       do is=1,ns
       do j=1,4
         write(*,*)is,j,"  flux = ",   &
         diff_out(is,j,1)*ky*rlns(is)+(diff_out(is,j,2)/sqrt_two+diff_out(is,j,3))*ky*rlts(is) &
         +diff_out(is,j,4)*ky*vpar_shear_s(is)/sqrt_two+conv_out(is,j)*zs(is)/taus(is),flux_out(is,j)
       enddo
       enddo
!
      END SUBROUTINE get_QL_weights
!
!    
!***********************
! start of get_wavefunction
      SUBROUTINE get_wavefunction
!
      USE gftm_dimensions
      USE gftm_global
      USE gftm_sgrid
!
      IMPLICIT NONE
!
      INTEGER :: n,i,j,k,np,npi,j0,imax
      REAL :: dx,hp0
      REAL :: hp(nb,max_plot)
      REAL :: xp(max_plot)
!
! set up the theta-grid
! npi is the number of pi intervals for the plot from -npi Pi to +npi Pi
     npi=9  !npi <= 9 limited by max_plot = 2*npi*ms/8+1
     np = ms/8   ! number of points per 1/2 period = 16 for ms=12
     if(igeo.eq.0)then
       dx = REAL(npi)*2.0*pi/REAL(max_plot-1)
       do i=1,max_plot
         xp(i) = -REAL(npi)*pi + REAL(i-1)*dx
         plot_angle_out(i) = xp(i)        
       enddo
     else
! general geometry case 0<y<Ly one loop counterclockwise 
! y is the the straight field line coordinant of the Hermite basis
! t_s is the mapping of the original theta coordinate 0 < t_s < -2Pi with t_s(0)=0, t_s(ms)=-2 Pi
! this will be used for plotting, Note that t_s has the opposite sign to y.
       dx = 2.0*pi/(y(ms)*width_in)
       j0 = npi*np+1 ! the index of the midpoint
       xp(j0) = 0.0
       plot_angle_out(j0) = 0.0
       j = 0  ! j is the local index for one circuit poloidally
       k = 0  ! counts the number of 2 pi  intervals
       imax=np*npi
       do i=1,imax
           j = j+1
           if(j.gt.2*np)then
             j = j - 2*np
             k = k + 1
           endif
           xp(j0+i) = (REAL(k)*y(ms) + y(4*j))*dx  ! remember y is positive 0 <= y <= Ly
           xp(j0-i) = -(REAL(k+1)*y(ms) - y(ms-4*j))*dx ! ok for up/down assymetric cases
           plot_angle_out(j0+i) = -(REAL(k)*t_s(ms) + t_s(4*j))  ! remember t_s is negative  0 <= t_s <= -2 Pi
           plot_angle_out(j0-i) = REAL(k+1)*t_s(ms) - t_s(ms-4*j)
       enddo
     endif
 !     write(*,*)"t_s = ",(t_s(i),i=0,ms)
 !   do i=1,imax
 !     write(*,*)i,"xp=",xp(i),"tp=",plot_angle_out(i)
 !   enddo
! compute the hermite polynomials on the theta-grid using recursion
      hp0 = sqrt_two/pi**0.25
      do i=1,max_plot
       hp(1,i) = hp0*EXP(-xp(i)*xp(i)/2.0)
       hp(2,i) = xp(i)*sqrt_two*hp0*EXP(-xp(i)*xp(i)/2.0)
       if(nbasis.gt.2)then
         do j=3,nbasis
          hp(j,i) = xp(i)*SQRT(2.0/REAL(j-1))*hp(j-1,i) &
           - SQRT(REAL(j-2)/REAL(j-1))*hp(j-2,i)
         enddo
       endif
      enddo
!      compute the fields on the theta-grid
      do n=1,nmodes_out
       do k=1,3
!        write(*,*)"field_weight_out=",(field_weight_out(n,k,j),j=1,nbasis)
        do i=1,max_plot
         plot_field_out(n,k,i) = 0.0
         if((k.eq.1).or.(k.eq.2.and.use_bper_in).or.(k.eq.3.and.use_bpar_in))then
          do j=1,nbasis
           plot_field_out(n,k,i) = plot_field_out(n,k,i) &
           + field_weight_out(n,k,j)*hp(j,i)
          enddo
         endif
        enddo
       enddo
      enddo
!
      END SUBROUTINE get_wavefunction
!
