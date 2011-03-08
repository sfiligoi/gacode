!
      SUBROUTINE gstotal
!
      USE gks_var
      IMPLICIT NONE
!
      INCLUDE 'input.m'
!
      INTEGER iblock,save_nstep,save_icontinue,nstep_tot,i
      REAL save_delt,test
!
      if(igks_model.lt.0.or.igks_model.gt.1)then
        write(*,*)"use igks_model=0 for original gks, =1 for tglf "
      endif
      if(igks_model.eq.0)then
         save_delt=delt
         save_nstep=nstep
         save_icontinue=icontinue
! adjust timestep and nstep with the total diamagnetic drift frequency
         test = ABS(fprim1 + tprim1)
         test = MAX(test,ABS(fprim1*fprim3 + tprim3))
         if(ncspec2.ne.0)test=MAX(test,ABS(fprim1*fprim2 + tprim2))
         test = MAX(test,epsl)
         nstep = save_nstep + INT(10.0*test)
         delt = save_delt*4.0/test
         write(*,*)"test = ",test,"delt = ",delt," nstep =",nstep
! GKS has a storage limit of nstep<=6000 so a restart is needed to go further
         iblock=nstep/6000
         write(*,*)"iblock =",iblock
         if(iblock.gt.0)then
           nstep_tot=nstep
           nstep = 6000
           write(*,*)"nstep = ",nstep
           call gks
           icontinue=1
           do i=1,iblock
             if((i+1)*6000.gt.nstep_tot)nstep=nstep_tot-i*6000
             write(*,*)"nstep = ",nstep
             call gks
           enddo
         else
           call gks
         endif
         icontinue=save_icontinue
         nstep=save_nstep
         delt = save_delt
      endif
!
      if(igks_model.eq.1)call gks_tglf
!
      END SUBROUTINE gstotal
!
      SUBROUTINE gks_tglf
!
      USE gks_var
      USE tglf_pkg
      USE tglf_tg
      IMPLICIT NONE
!
      INCLUDE 'data_exp.m'
      INCLUDE 'input.m'
      INCLUDE 'gks_out.m'
      LOGICAL :: firstcall=.TRUE.
      LOGICAL :: new_gridpoint=.FALSE.
      INTEGER :: n
      REAL :: aky1_save=0.0
      REAL :: rmin_save=-1.0
!
      REAL :: particle_flux_tg(2,2),energy_flux_tg(2,2)
!
      if(firstcall)then
        call tglf_startup
        igeo_tg = igeo
!
        if(igks.le.2)iflux_tg=.TRUE.
        if(cbetae.gt.0.0)then
          use_bper_tg = .TRUE.
          use_mhd_rule_tg=.TRUE.
          if(i_bpar.eq.1)then
            use_bpar_tg = .TRUE.
            use_mhd_rule_tg=.FALSE.
          endif
        endif
!
        zs_tg(1)=-1.0
        zs_tg(2)=1.0
        zs_tg(3)=z2
        mass_tg(1)=amass3
        mass_tg(2)=1.0
        mass_tg(3)=amass2
        ns_tg=2
        if(ncspec2.eq.1)ns_tg=3
        CALL put_species(ns_tg,zs_tg,mass_tg) 
!
        CALL put_signs(sign_Bt_tg,sign_It_tg)
!
        nbasis_min_tg=nbasis_min 
        nbasis_max_tg=nbasis_max 
        ibranch_tg = ibranch
        if(ibranch.eq.-1)nmodes_tg=4
!        write(*,*)"switches",iflux_tg,use_bper_tg,use_bpar_tg,
!     >   ibranch_tg,nmodes_tg,nb_max_tg,nb_min_tg,nxgrid_tg,nky_tg
!   
      CALL put_switches(iflux_tg,use_bper_tg,use_bpar_tg,
     >  use_mhd_rule_tg,use_bisection_tg,ibranch_tg,
     >  nmodes_tg,nbasis_max_tg,nbasis_min_tg,nxgrid_tg,nky_tg)
!
        if(z3.gt.0)adiabatic_elec_tg=.TRUE.
!        write(*,*)"parameters",adiabatic_elec_tg,alpha_p_tg,alpha_e_tg, 
!     >   theta_trap_tg,xnu_fac_tg,debye_factor_tg
!
       CALL put_model_parameters(adiabatic_elec_tg,alpha_p_tg,
     > alpha_e_tg,alpha_kx0_tg,alpha_kx1_tg,alpha_quench_tg,
     > xnu_factor_tg,debye_factor_tg,etg_factor_tg,sat_rule_tg,
     > kygrid_model_tg,xnu_model_tg,vpar_model_tg,
     > vpar_shear_model_tg)
      endif
      firstcall=.FALSE.
!
!      find_width_tg=.FALSE.
      find_width_tg=.TRUE.
      if(aky1.ne.aky1_save)find_width_tg=.TRUE.
      new_gridpoint=.FALSE.
      if(igeo.eq.0.and.rmin_save.ne.eps/epsa)then
        new_gridpoint=.TRUE.
        find_width_tg=.TRUE.
      endif
      if(igeo.eq.1.and.rmin_save.ne.rmin_loc)then
        new_gridpoint=.TRUE.
        find_width_tg=.TRUE.
      endif 
      aky1_save=aky1
      if(igeo.eq.0)then
        rmin_save=eps/epsa
      else
        rmin_save=rmin_loc
      endif
!
      ky_tg = aky1/SQRT(2.0/temp3)
!
      CALL put_kys(ky_tg)
!
!      write(*,*)"find_width_tg=",find_width_tg
!      write(*,*)"width_max_tg=",width_max_tg
!      write(*,*)"width_min=",width_min_tg,"nwidth=",nwidth_tg
!      write(*,*)"width",find_width_tg,width_max_tg,
!     >  width_min_tg,nwidth_tg,find_width_tg
      if(find_width_tg)CALL put_gaussian_width(width_max_tg,
     >  width_min_tg,nwidth_tg,find_width_tg)
!
      rlns_tg(1) = fprim1*fprim3
      rlns_tg(2) = fprim1
      rlns_tg(3) = fprim1*fprim2
      rlts_tg(1) = tprim3
      rlts_tg(2) = tprim1
      rlts_tg(3) = tprim1
      vexb_shear_tg= egamma
      vpar_shear_tg(1)=sign_Bt_tg*uprim3/SQRT(2.0/temp3)
      vpar_shear_tg(2)=sign_Bt_tg*uprim1/SQRT(2.0/temp3)
      vpar_shear_tg(3)=sign_Bt_tg*uprim2/SQRT(2.0/temp3)
!      write(*,*)"gradients",rlns_tg(1),rlns_tg(2),rlns_tg(3)
!      write(*,*)rlts_tg(1),rlts_tg(2),rlts_tg(3)
!      write(*,*)"vexb_shear_tg=",vexb_shear_tg
!      write(*,*)vpar_shear_tg(1),vpar_shear_tg(2)
!
      CALL put_gradients(rlns_tg,rlts_tg,vpar_shear_tg,
     > vexb_shear_tg)
!
      if(igeo_tg.eq.0)then
!
        taus_tg(1) = 1.0
        taus_tg(2) = 1.0/temp3
        taus_tg(3) = 1.0/temp3
        as_tg(1) = 1.0
        as_tg(2) = an1
        as_tg(3) = an2
        vpar_tg(1)=sign_Bt_tg*mach3
        vpar_tg(2)=sign_Bt_tg*mach3
        vpar_tg(3)=sign_Bt_tg*mach2
        if(ncspec2.eq.0)as_tg(3)=0.0  ! treat impurity as dilution only
        betae_tg = beta*temp3   
        xnue_tg = cnewk3*SQRT(2.0/temp3)*vnewk3
        zeff_tg = zeff
        debye_tg = debyelorhos
!        write(*,*)"averages",betae_tg,xnue_tg,zeff_tg,debye_tg
!        write(*,*)taus_tg(1),taus_tg(2),taus_tg(3)
!        write(*,*)as_tg(1),as_tg(2),as_tg(3)
!
      CALL put_averages(taus_tg,as_tg,vpar_tg,betae_tg,xnue_tg,
     >  zeff_tg,debye_tg)
!
        rmin_tg = eps/epsa
        rmaj_tg = 1.0/epsa
        q_tg = 2.0*epsa/pk
        shat_tg = shat
        alpha_tg = shift
        if(igyro_fix.eq.1)b_model_tg=1
!        write(*,*)"geometry",rmin_tg,rmaj_tg,q_tg,shat_tg,alpha_tg
!
        if(new_gridpoint)CALL put_s_alpha_geometry(rmin_tg,rmaj_tg,q_tg,
     > shat_tg,alpha_tg,xwell_tg,theta0_tg,b_model_tg,ft_model_tg)
!        write(*,*)"debug",new_gridpoint,xwell_tg,theta0_tg,
!     >   b_model_tg,ft_model_tg
!
      elseif(igeo_tg.eq.1)then
!
        taus_tg(1) = 1.0
        taus_tg(2) = tiote_loc
        taus_tg(3) = tiote_loc
        as_tg(1) = 1.0
        as_tg(2) = nione_loc
        as_tg(3) = (1.0 - (nione_loc + fastionfrac_loc))/z2
        if(ncspec2.eq.0)as_tg(3)=0.0
        betae_tg = beta/tiote_loc
        xnue_tg = xnu_loc
        zeff_tg = zeff
        debye_tg = debyelorhos
!        write(*,*)"averages",betae_tg,xnue_tg,zeff_tg,debye_tg
!        write(*,*)taus_tg(1),taus_tg(2),taus_tg(3)
!        write(*,*)as_tg(1),as_tg(2),as_tg(3)
!
      CALL put_averages(taus_tg,as_tg,vpar_tg,betae_tg,
     >  xnue_tg,zeff_tg,debye_tg)
!
        rmin_tg = rmin_loc
        rmaj_tg = rmaj0_loc
        zmaj_tg=0.0
        q_tg = q_loc
        q_prime_tg = q_prime_loc
        p_prime_tg = p_prime_loc 
        drmajdx_tg = shift_loc
        drmindx_tg=1.0
        dzmajdx_tg=0.0
        kappa_tg = kappa_loc
        s_kappa_tg = s_kappa_loc
        delta_tg = delta_loc
        s_delta_tg = s_delta_loc*SQRT(1.0-delta_loc**2)  ! gyro conventions now used
        zeta_tg=0.0
        s_zeta_tg=0.0
!        write(*,*)"inout rmin_tg=",rmin_tg
!        write(*,*)"inout rmaj_tg=",rmaj_tg
!        write(*,*)"inout q_tg=",q_tg
!        write(*,*)"inout q_prime_tg=",q_prime_tg
!        write(*,*)"inout p_prime_tg=",p_prime_tg
!        write(*,*)"inout shift_tg=",shift_tg
!        write(*,*)"inout kappa_tg=",kappa_tg
!        write(*,*)"inout s_kappa_tg=",s_kappa_tg
!        write(*,*)"inout delta_tg=",delta_tg
!        write(*,*)"inout s_delta_tg=",s_delta_tg
!
        if(new_gridpoint)then
          CALL put_Miller_geometry(rmin_tg,rmaj_tg,zmaj_tg,drmindx_tg,
     >    drmajdx_tg,dzmajdx_tg,kappa_tg,s_kappa_tg,delta_tg,s_delta_tg,
     >    zeta_tg,s_zeta_tg,q_tg,q_prime_tg,p_prime_tg)
        endif
      else
        write(*,*)"igeo_tg invalid",igeo_tg
        stop
      endif
!
      if(save_tglf)then
        CALL write_tglf_overwrite
        save_tglf = .FALSE.
      endif
!
      if(overwrite_tglf)CALL read_tglf_overwrite
!
      CALL tglf
!
!      write(*,*)"quench",alpha_quench_tg,vexb_shear_tg
      kys(1)=ky_tg
      agammas(2)=get_growthrate(2)
      afreqs(2)=get_frequency(2)
      ne_te_phase_k(2)=get_ne_te_phase(2)
      agammas(3)=get_growthrate(1)
      afreqs(3)=get_frequency(1)
      ne_te_phase_k(3)=get_ne_te_phase(1)
      n=2
      if(agammas(3).gt.agammas(2))n=3
      agammas(1)=agammas(n)
      afreqs(1)=afreqs(n)
      dgammas(1)=1.0e-10
      ne_te_phase_k(1)=ne_te_phase_k(n)
!      write(*,*)"width=",get_gaussian_width()
      if(iflux_tg)then
        particle_flux_tg(1,1) = get_QL_particle_flux(1,1,1)
        particle_flux_tg(1,2) = get_QL_particle_flux(1,2,1)
        particle_flux_tg(2,1) = get_QL_particle_flux(2,1,1)
        particle_flux_tg(2,2) = get_QL_particle_flux(2,2,1)
        energy_flux_tg(1,1) = get_QL_energy_flux(1,1,1)
        energy_flux_tg(1,2) = get_QL_energy_flux(1,2,1)
        energy_flux_tg(2,1) = get_QL_energy_flux(2,1,1)
        energy_flux_tg(2,2) = get_QL_energy_flux(2,2,1)
!
! set fluxes for the most unstable mode in GKS units
!
        n=1
        if(agammas(2).gt.agammas(3))n=2
        peflxa=particle_flux_tg(n,1)*SQRT(2.0/temp3)/temp3
        eeflxa=energy_flux_tg(n,1)*SQRT(2.0/temp3)/temp3
        eiflxa=energy_flux_tg(n,2)*SQRT(2.0/temp3)/temp3
        if(ncspec2.ne.0)then
          eiflxa=eiflxa+get_QL_energy_flux(n,3,1)*SQRT(2.0/temp3)/temp3
        endif
! save the fluctuation intensities summed over nmodes_tg
        phi_bar_k = 0.0
        ne_bar_k = 0.0
        te_bar_k = 0.0
        ti_bar_k = 0.0
        do n=1,nmodes_tg
          phi_bar_k = phi_bar_k + get_phi_bar(n)
          ne_bar_k = ne_bar_k + get_N_bar(n,1)
          te_bar_k = te_bar_k + get_T_bar(n,1)
          ti_bar_k = ti_bar_k + get_T_bar(n,2)
        enddo
      endif
      mhd_DR_k = 0.0
      if(igeo_tg.eq.1)then
        mhd_DR_k = get_DR()
      endif
!   
!      if(iflux_tg.ne.0)then
!        write(*,*) 'particle_flux_tg(1,1) = ',particle_flux_tg(1,1)
!        write(*,*) 'particle_flux_tg(1,2) = ',particle_flux_tg(1,2)
!        write(*,*) 'particle_flux_tg(2,1) = ',particle_flux_tg(2,1)
!        write(*,*) 'particle_flux_tg(2,2) = ',particle_flux_tg(2,2)
!        write(*,*) 'energy_flux_tg(1,1) = ',energy_flux_tg(1,1)
!        write(*,*) 'energy_flux_tg(1,2) = ',energy_flux_tg(1,2)
!        write(*,*) 'energy_flux_tg(2,1) = ',energy_flux_tg(2,1)
!        write(*,*) 'energy_flux_tg(2,2) = ',energy_flux_tg(2,2)
!      endif           
 900  CONTINUE
 20   FORMAT (a)
!
!  
      END SUBROUTINE gks_tglf  
!
       SUBROUTINE write_tglf_overwrite
!
       USE tglf_pkg
       USE tglf_tg
       IMPLICIT NONE
!
       INCLUDE 'input.m'
!
      
!
! Open the tglf_overwrite file
!
       OPEN (unit=3,file='tglf_overwrite',status='unknown')
!
       WRITE(3,*)"&tglfin"
!  averages
       WRITE(3,*)"taus_tg(1)=",taus_tg(1)
       WRITE(3,*)"taus_tg(2)=",taus_tg(2)
       WRITE(3,*)"taus_tg(3)=",taus_tg(3)
       WRITE(3,*)"as_tg(1)=",as_tg(1)
       WRITE(3,*)"as_tg(2)=",as_tg(2)
       WRITE(3,*)"as_tg(3)=",as_tg(3)
       WRITE(3,*)"vpar_tg(1)=",vpar_tg(1)
       WRITE(3,*)"vpar_tg(2)=",vpar_tg(2)
       WRITE(3,*)"vpar_tg(3)=",vpar_tg(3)
       WRITE(3,*)"betae_tg=",betae_tg
       WRITE(3,*)"xnue_tg=",xnue_tg
       WRITE(3,*)"zeff_tg=",zeff_tg
       WRITE(3,*)"debye_tg=",debye_tg
!  gradients
       WRITE(3,*)"rlns_tg(1)=",rlns_tg(1)
       WRITE(3,*)"rlns_tg(2)=",rlns_tg(2)
       WRITE(3,*)"rlns_tg(3)=",rlns_tg(3)
       WRITE(3,*)"rlts_tg(1)=",rlts_tg(1)
       WRITE(3,*)"rlts_tg(2)=",rlts_tg(2)
       WRITE(3,*)"rlts_tg(3)=",rlts_tg(3)
       WRITE(3,*)"vpar_shear_tg(1)=",vpar_shear_tg(1)
       WRITE(3,*)"vpar_shear_tg(2)=",vpar_shear_tg(2)
       WRITE(3,*)"vpar_shear_tg(3)=",vpar_shear_tg(3)
       WRITE(3,*)"vexb_shear_tg=",vexb_shear_tg
! species
       WRITE(3,*)"mass_tg(1)=",mass_tg(1)
       WRITE(3,*)"mass_tg(2)=",mass_tg(2)
       WRITE(3,*)"mass_tg(3)=",mass_tg(3)
       WRITE(3,*)"zs_tg(1)=",zs_tg(1)
       WRITE(3,*)"zs_tg(2)=",zs_tg(2)
       WRITE(3,*)"zs_tg(3)=",zs_tg(3)
! common geometry
       WRITE(3,*)"rmin_tg=",rmin_tg
       WRITE(3,*)"rmaj_tg=",rmaj_tg
       WRITE(3,*)"q_tg=",q_tg
       if(igeo_tg.eq.0)then
! Shifted cicle inputs
       WRITE(3,*)"theta0_tg=",theta0_tg
       WRITE(3,*)"shat_tg=",shat_tg
       WRITE(3,*)"alpha_tg=",alpha_tg
       WRITE(3,*)"xwell_tg=",xwell_tg
       elseif(igeo_tg.eq.1)then
! Miller inputs
       WRITE(3,*)"drmajdx_tg=",drmajdx_tg
       WRITE(3,*)"kappa_tg=",kappa_tg
       WRITE(3,*)"s_kappa_tg=",s_kappa_tg
       WRITE(3,*)"delta_tg=",delta_tg
       WRITE(3,*)"s_delta_tg=",s_delta_tg
       WRITE(3,*)"zeta_tg=",zeta_tg
       WRITE(3,*)"s_zeta_tg=",s_zeta_tg
       WRITE(3,*)"q_prime_tg=",q_prime_tg
       WRITE(3,*)"p_prime_tg=",p_prime_tg
       endif
!
      WRITE(3,*)"/"
!
       CLOSE(3)
!
       RETURN
       END SUBROUTINE write_tglf_overwrite
!
      SUBROUTINE read_tglf_overwrite
!
      USE tglf_pkg
      USE tglf_tg
      IMPLICIT NONE
!
      INCLUDE 'input.m'
!
! read the namelist 
!
        OPEN (unit=3,file='tglf_overwrite',status='old')
        READ(3,nml=tglfin)
        CLOSE(3)
!
        CALL put_species(ns_tg,zs_tg,mass_tg) 
        CALL put_gradients(rlns_tg,rlts_tg,vpar_shear_tg,
     > vexb_shear_tg)
        CALL put_averages(taus_tg,as_tg,vpar_tg,betae_tg,xnue_tg,
     >  zeff_tg,debye_tg)
!
      if(igeo_tg.eq.0)then
        CALL put_s_alpha_geometry(rmin_tg,rmaj_tg,q_tg,
     >  shat_tg,alpha_tg,xwell_tg,theta0_tg,b_model_tg,ft_model_tg)
      elseif(igeo_tg.eq.1)then
        CALL put_Miller_geometry(rmin_tg,rmaj_tg,zmaj_tg,drmindx_tg,
     >  drmajdx_tg,dzmajdx_tg,kappa_tg,s_kappa_tg,delta_tg,s_delta_tg,
     >  zeta_tg,s_zeta_tg,q_tg,q_prime_tg,p_prime_tg)
      endif
!
      RETURN
      END SUBROUTINE read_tglf_overwrite
!__________________________________________
!
      SUBROUTINE tglf_startup
!
      USE tglf_pkg
      USE tglf_tg
      IMPLICIT NONE
!
      INCLUDE 'input.m'    
!
!  1.81  version (APS07)
!
       if(tglf_defaults.eq.1)then
         xnu_model_tg=1
         alpha_p_tg=0.0
         alpha_quench_tg=0.0
         use_bpar_tg=.FALSE.
         use_mhd_rule_tg=.TRUE.
         filter_tg=2.0
       endif
!
! 1.82 version (TGLF09)
!
       if(tglf_defaults.eq.2)then
         xnu_model_tg=2
         alpha_quench_tg=0.0
         use_bpar_tg=.FALSE.
         use_mhd_rule_tg=.TRUE.
         filter_tg=2.0
       endif
!
! 1.93 version (TGLF10)
!
       if(tglf_defaults.eq.3)then
         xnu_model_tg=2
         alpha_quench_tg=1.0
         use_bpar_tg=.FALSE.
         use_mhd_rule_tg=.TRUE.
         filter_tg=2.0
       endif
!
      if(tglf_defaults.eq.0)then
!
! read the namelist 
!
        OPEN (unit=3,file='tglfin',status='old')
        READ(3,nml=tglfin)
        CLOSE(3)
      endif
!
      sign_It_tg=1.0
      sign_Bt_tg=bt_exp/ABS(bt_exp)
      CALL put_signs(sign_Bt_tg,sign_It_tg)
!
      CALL put_rare_switches(theta_trapped_tg,park_tg,ghat_tg
     > ,gchat_tg,wd_zero_tg,Linsker_factor_tg,gradB_factor_tg,
     > filter_tg,damp_psi_tg,damp_sig_tg)      
!
      RETURN
      END SUBROUTINE tglf_startup
!

