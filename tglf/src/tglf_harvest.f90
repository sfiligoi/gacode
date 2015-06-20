      SUBROUTINE tglf_harvest

      USE tglf_species
      USE tglf_coeff
      USE tglf_global
      USE tglf_pkg
      USE tglf_tg

      INTEGER :: j,ierr

      CHARACTER(LEN=65507) :: harvest_sendline
      CHARACTER(LEN=3) :: NUM
      CHARACTER NUL
      PARAMETER(NUL = CHAR(0))

      ierr=set_harvest_verbose(1)
      ierr=set_harvest_table('TGLF_harvest?'//NUL)
      ierr=set_harvest_host('127.0.0.1'//NUL)
      ierr=set_harvest_port(32000)
      ierr=init_harvest(harvest_sendline,65507)

      WRITE(*,*)'===HARVEST starts==='

      ierr=set_harvest_payload_dbl(harvest_sendline,'ky'//NUL,ky_in)
      ierr=set_harvest_payload_int(harvest_sendline,'sign_Bt'//NUL,INT(sign_Bt_in))
      ierr=set_harvest_payload_int(harvest_sendline,'sign_Ip'//NUL,INT(sign_It_in))
      ierr=set_harvest_payload_dbl(harvest_sendline,'vexb_shear'//NUL,vexb_shear_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'vexb'//NUL,vexb_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'betae'//NUL,betae_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'xnue'//NUL,xnue_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'zeff'//NUL,zeff_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'debye'//NUL,debye_in)
      ierr=set_harvest_payload_bol(harvest_sendline,'iflux'//NUL,iflux_in)

      ierr=set_harvest_payload_bol(harvest_sendline,'use_bper'//NUL,use_bper_in)
      ierr=set_harvest_payload_bol(harvest_sendline,'use_bpar'//NUL,use_bpar_in)
      ierr=set_harvest_payload_bol(harvest_sendline,'use_mhd_rule'//NUL,use_mhd_rule_in)
      ierr=set_harvest_payload_bol(harvest_sendline,'use_bisection'//NUL,use_bisection_in)
      ierr=set_harvest_payload_int(harvest_sendline,'ibranch'//NUL,ibranch_in)
      ierr=set_harvest_payload_int(harvest_sendline,'nmodes'//NUL,nmodes_in)
      ierr=set_harvest_payload_int(harvest_sendline,'nbasis_max'//NUL,nbasis_max_in)
      ierr=set_harvest_payload_int(harvest_sendline,'nbasis_min'//NUL,nbasis_min_in)
      ierr=set_harvest_payload_int(harvest_sendline,'nxgrid'//NUL,nxgrid_in)
      ierr=set_harvest_payload_int(harvest_sendline,'nky'//NUL,nky_in)

      ierr=set_harvest_payload_dbl(harvest_sendline,'theta_trapped'//NUL,theta_trapped_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'park'//NUL,park_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'ghat'//NUL,ghat_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'gchat'//NUL,gchat_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'wd_zero'//NUL,wd_zero_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'Linsker_factor'//NUL,Linsker_factor_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'gradB_factor'//NUL,gradB_factor_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'filter'//NUL,filter_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'damp_psi'//NUL,damp_psi_in)

      ierr=set_harvest_payload_bol(harvest_sendline,'adiabatic_elec'//NUL,adiabatic_elec_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'alpha_mach'//NUL,alpha_mach_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'alpha_p'//NUL,alpha_p_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'alpha_e'//NUL,alpha_e_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'alpha_quench'//NUL,alpha_quench_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'xnu_factor'//NUL,xnu_factor_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'debye_factor'//NUL,debye_factor_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'etg_factor'//NUL,etg_factor_in)
      ierr=set_harvest_payload_int(harvest_sendline,'sat_rule'//NUL,sat_rule_in)
      ierr=set_harvest_payload_int(harvest_sendline,'xnu_model'//NUL,xnu_model_in)
      ierr=set_harvest_payload_int(harvest_sendline,'kygrid_model'//NUL,kygrid_model_in)
      ierr=set_harvest_payload_int(harvest_sendline,'vpar_model'//NUL,vpar_model_in)
      ierr=set_harvest_payload_int(harvest_sendline,'vpar_shear_model'//NUL,vpar_shear_model_in)

      ierr=set_harvest_payload_dbl(harvest_sendline,'rmin_sa'//NUL,rmin_sa)
      ierr=set_harvest_payload_dbl(harvest_sendline,'rmaj_sa'//NUL,rmaj_sa)
      ierr=set_harvest_payload_dbl(harvest_sendline,'q_sa'//NUL,q_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'shat_sa'//NUL,shat_sa)
      ierr=set_harvest_payload_dbl(harvest_sendline,'alpha_sa'//NUL,alpha_sa)
      ierr=set_harvest_payload_dbl(harvest_sendline,'xwell_sa'//NUL,xwell_sa)
      ierr=set_harvest_payload_dbl(harvest_sendline,'theta0_sa'//NUL,theta0_sa)
      ierr=set_harvest_payload_int(harvest_sendline,'b_model_sa'//NUL,b_model_sa)
      ierr=set_harvest_payload_int(harvest_sendline,'ft_model_sa'//NUL,ft_model_sa)

      ierr=set_harvest_payload_dbl(harvest_sendline,'rmin_loc'//NUL,rmin_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'rmaj_loc'//NUL,rmaj_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'zmaj_loc'//NUL,zmaj_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'q_loc'//NUL,q_in)
      ierr=set_harvest_payload_dbl(harvest_sendline,'p_prime_loc'//NUL,p_prime_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'q_prime_loc'//NUL,q_prime_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'drmindx_loc'//NUL,drmindx_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'drmajdx_loc'//NUL,drmajdx_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'dzmajdx_loc'//NUL,dzmajdx_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'kappa_loc'//NUL,kappa_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'s_kappa_loc'//NUL,s_kappa_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'delta_loc'//NUL,delta_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'s_delta_loc'//NUL,s_delta_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'zeta_loc'//NUL,zeta_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'s_zeta_loc'//NUL,s_zeta_loc)
      ierr=set_harvest_payload_dbl(harvest_sendline,'kx0_loc'//NUL,kx0_loc)

!      write(*,*)'q_fourier=',q_fourier_in
!      write(*,*)'q_fourier=',q_in
!      write(*,*)'p_prime_fourier=',p_prime_fourier_in
!      write(*,*)'q_prime_fourier=',q_prime_fourier_in
!      write(*,*)'nfourier=',nfourier_in
!      write(*,*)'fourier=',fourier_in

      write(*,*)"---------------------------"

      ierr=set_harvest_payload_int(harvest_sendline,'ns'//NUL,ns_in)

      do j=1,ns_in
         if (j < 10) then
            write (NUM, "(A1,I01)") "_", j
         else
            write (NUM, "(A1,I02)") "_", j
         endif

         ierr=set_harvest_payload_dbl(harvest_sendline,'Z'//NUM//NUL,zs_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'A'//NUM//NUL,mass_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'rlns'//NUM//NUL,rlns_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'rlts'//NUM//NUL,rlts_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'vpar_shear'//NUM//NUL,vpar_shear_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'vns_shear'//NUM//NUL,vns_shear_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'vts_shear'//NUM//NUL,vts_shear_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'taus'//NUM//NUL,taus_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'as'//NUM//NUL,as_in(j))
         ierr=set_harvest_payload_dbl(harvest_sendline,'vpar'//NUM//NUL,vpar_in(j))

!        electrostatic
         ierr=set_harvest_payload_dbl(harvest_sendline,'particle_flux'//NUM//NUL,get_particle_flux(j,1))
         ierr=set_harvest_payload_dbl(harvest_sendline,'energy_flux'//NUM//NUL,get_energy_flux(j,1))
         ierr=set_harvest_payload_dbl(harvest_sendline,'stress_par'//NUM//NUL,get_stress_par(j,1))
         ierr=set_harvest_payload_dbl(harvest_sendline,'stress_tor'//NUM//NUL,get_stress_tor(j,1))
         ierr=set_harvest_payload_dbl(harvest_sendline,'exchange'//NUM//NUL,get_exchange(j,1))
         write(*,*)"---------------------------"
      enddo

      ierr=set_harvest_payload_dbl(harvest_sendline,'gmax'//NUL,get_growthrate(1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'fmax'//NUL,get_frequency(1))

      WRITE(*,*)'===HARVEST sends==='

      ierr=harvest_send(harvest_sendline)

      WRITE(*,*)'===HARVEST ends==='

!!$!
!!$!
!!$!     MODULE tglf_species 
!!$! 
!!$! species parameters
!!$!
!!$      IF(ALLOCATED(ei_exch))DEALLOCATE(ei_exch)
!!$      IF(ALLOCATED(resist))DEALLOCATE(resist)
!!$      IF(ALLOCATED(zs))DEALLOCATE(zs)
!!$      IF(ALLOCATED(mass))DEALLOCATE(mass)
!!$      IF(ALLOCATED(vs))DEALLOCATE(vs)
!!$      IF(ALLOCATED(rlts))DEALLOCATE(rlts)
!!$      IF(ALLOCATED(rlns))DEALLOCATE(rlns)
!!$      IF(ALLOCATED(vpar_shear_s))DEALLOCATE(vpar_shear_s)
!!$      IF(ALLOCATED(as))DEALLOCATE(as)
!!$      IF(ALLOCATED(taus))DEALLOCATE(taus)
!!$      IF(ALLOCATED(vpar_s))DEALLOCATE(vpar_s)
!!$
!!$!      write(*,*)"deallocated species"
!!$!
!!$!
!!$!-------------------------------------------------
!!$!
!!$!      MODULE tglf_coeff
!!$!
!!$! store the hermite basis matrix coefficients
!!$!
!!$! ave_h
!!$      IF(ALLOCATED(ave_hn))DEALLOCATE(ave_hn)
!!$      IF(ALLOCATED(ave_hp1))DEALLOCATE(ave_hp1)
!!$      IF(ALLOCATED(ave_hp3))DEALLOCATE(ave_hp3)
!!$      IF(ALLOCATED(ave_hr11))DEALLOCATE(ave_hr11)
!!$      IF(ALLOCATED(ave_hr13))DEALLOCATE(ave_hr13)
!!$      IF(ALLOCATED(ave_hr33))DEALLOCATE(ave_hr33)
!!$      IF(ALLOCATED(ave_hw113))DEALLOCATE(ave_hw113)
!!$      IF(ALLOCATED(ave_hw133))DEALLOCATE(ave_hw133)
!!$      IF(ALLOCATED(ave_hw333))DEALLOCATE(ave_hw333)
!!$      IF(ALLOCATED(ave_ht1))DEALLOCATE(ave_ht1)
!!$      IF(ALLOCATED(ave_ht3))DEALLOCATE(ave_ht3)
!!$      IF(ALLOCATED(ave_hu1))DEALLOCATE(ave_hu1)
!!$      IF(ALLOCATED(ave_hu3))DEALLOCATE(ave_hu3)
!!$      IF(ALLOCATED(ave_hu33))DEALLOCATE(ave_hu33)
!!$      IF(ALLOCATED(ave_hu3ht1))DEALLOCATE(ave_hu3ht1)
!!$      IF(ALLOCATED(ave_hu3ht3))DEALLOCATE(ave_hu3ht3)
!!$      IF(ALLOCATED(ave_hu33ht1))DEALLOCATE(ave_hu33ht1)
!!$      IF(ALLOCATED(ave_hu33ht3))DEALLOCATE(ave_hu33ht3)
!!$      IF(ALLOCATED(ave_hninv))DEALLOCATE(ave_hninv)
!!$      IF(ALLOCATED(ave_hp1inv))DEALLOCATE(ave_hp1inv)
!!$      IF(ALLOCATED(ave_hp3inv))DEALLOCATE(ave_hp3inv)
!!$      IF(ALLOCATED(ave_gradhp1))DEALLOCATE(ave_gradhp1)
!!$      IF(ALLOCATED(ave_gradhr11))DEALLOCATE(ave_gradhr11)
!!$      IF(ALLOCATED(ave_gradhr13))DEALLOCATE(ave_gradhr13)
!!$      IF(ALLOCATED(ave_gradhp1p1))DEALLOCATE(ave_gradhp1p1)
!!$      IF(ALLOCATED(ave_gradhr11p1))DEALLOCATE(ave_gradhr11p1)
!!$      IF(ALLOCATED(ave_gradhr13p1))DEALLOCATE(ave_gradhr13p1)
!!$      IF(ALLOCATED(ave_gradhu1))DEALLOCATE(ave_gradhu1)
!!$      IF(ALLOCATED(ave_gradhu3))DEALLOCATE(ave_gradhu3)
!!$      IF(ALLOCATED(ave_hnp0))DEALLOCATE(ave_hnp0)
!!$      IF(ALLOCATED(ave_hp1p0))DEALLOCATE(ave_hp1p0)
!!$      IF(ALLOCATED(ave_hp3p0))DEALLOCATE(ave_hp3p0)
!!$      IF(ALLOCATED(ave_hr11p0))DEALLOCATE(ave_hr11p0)
!!$      IF(ALLOCATED(ave_hr13p0))DEALLOCATE(ave_hr13p0)
!!$      IF(ALLOCATED(ave_hr33p0))DEALLOCATE(ave_hr33p0)
!!$      IF(ALLOCATED(ave_hnb0))DEALLOCATE(ave_hnb0)
!!$      IF(ALLOCATED(ave_hp1b0))DEALLOCATE(ave_hp1b0)
!!$      IF(ALLOCATED(ave_hp3b0))DEALLOCATE(ave_hp3b0)
!!$      IF(ALLOCATED(ave_hr11b0))DEALLOCATE(ave_hr11b0)
!!$      IF(ALLOCATED(ave_hr13b0))DEALLOCATE(ave_hr13b0)
!!$      IF(ALLOCATED(ave_hr33b0))DEALLOCATE(ave_hr33b0)
!!$      IF(ALLOCATED(ave_hw113b0))DEALLOCATE(ave_hw113b0)
!!$      IF(ALLOCATED(ave_hw133b0))DEALLOCATE(ave_hw133b0)
!!$      IF(ALLOCATED(ave_hw333b0))DEALLOCATE(ave_hw333b0)
!!$      IF(ALLOCATED(ave_hnbp))DEALLOCATE(ave_hnbp)
!!$      IF(ALLOCATED(ave_hp1bp))DEALLOCATE(ave_hp1bp)
!!$      IF(ALLOCATED(ave_hp3bp))DEALLOCATE(ave_hp3bp)
!!$      IF(ALLOCATED(ave_hr11bp))DEALLOCATE(ave_hr11bp)
!!$      IF(ALLOCATED(ave_hr13bp))DEALLOCATE(ave_hr13bp)
!!$      IF(ALLOCATED(ave_hr33bp))DEALLOCATE(ave_hr33bp)
!!$      IF(ALLOCATED(ave_hw113bp))DEALLOCATE(ave_hw113bp)
!!$      IF(ALLOCATED(ave_hw133bp))DEALLOCATE(ave_hw133bp)
!!$      IF(ALLOCATED(ave_hw333bp))DEALLOCATE(ave_hw333bp)
!!$      IF(ALLOCATED(ave_hnp0b0))DEALLOCATE(ave_hnp0b0)
!!$      IF(ALLOCATED(ave_gradhp1p0))DEALLOCATE(ave_gradhp1p0)
!!$      IF(ALLOCATED(ave_gradhr11p0))DEALLOCATE(ave_gradhr11p0)
!!$      IF(ALLOCATED(ave_gradhr13p0))DEALLOCATE(ave_gradhr13p0)
!!$      IF(ALLOCATED(ave_wdhu3ht1))DEALLOCATE(ave_wdhu3ht1)
!!$      IF(ALLOCATED(ave_wdhu3ht3))DEALLOCATE(ave_wdhu3ht3)
!!$      IF(ALLOCATED(ave_wdhu33ht1))DEALLOCATE(ave_wdhu33ht1)
!!$      IF(ALLOCATED(ave_wdhu33ht3))DEALLOCATE(ave_wdhu33ht3)
!!$      IF(ALLOCATED(ave_modwdhu3))DEALLOCATE(ave_modwdhu3)
!!$      IF(ALLOCATED(ave_modwdhu33))DEALLOCATE(ave_modwdhu33)
!!$      IF(ALLOCATED(ave_modwdhu3ht1))DEALLOCATE(ave_modwdhu3ht1)
!!$      IF(ALLOCATED(ave_modwdhu3ht3))DEALLOCATE(ave_modwdhu3ht3)
!!$      IF(ALLOCATED(ave_modwdhu33ht1))DEALLOCATE(ave_modwdhu33ht1)
!!$      IF(ALLOCATED(ave_modwdhu33ht3))DEALLOCATE(ave_modwdhu33ht3)
!!$      IF(ALLOCATED(ave_c_tor_par_hp1p0))DEALLOCATE(ave_c_tor_par_hp1p0)
!!$      IF(ALLOCATED(ave_c_tor_par_hr11p0))DEALLOCATE(ave_c_tor_par_hr11p0)
!!$      IF(ALLOCATED(ave_c_tor_par_hr13p0))DEALLOCATE(ave_c_tor_par_hr13p0)
!!$!      write(*,*)"deallocated ave_h"
!!$! ave_g
!!$      IF(ALLOCATED(ave_gn))DEALLOCATE(ave_gn)
!!$      IF(ALLOCATED(ave_gp1))DEALLOCATE(ave_gp1)
!!$      IF(ALLOCATED(ave_gp3))DEALLOCATE(ave_gp3)
!!$      IF(ALLOCATED(ave_gr11))DEALLOCATE(ave_gr11)
!!$      IF(ALLOCATED(ave_gr13))DEALLOCATE(ave_gr13)
!!$      IF(ALLOCATED(ave_gr33))DEALLOCATE(ave_gr33)
!!$      IF(ALLOCATED(ave_gw113))DEALLOCATE(ave_gw113)
!!$      IF(ALLOCATED(ave_gw133))DEALLOCATE(ave_gw133)
!!$      IF(ALLOCATED(ave_gw333))DEALLOCATE(ave_gw333)
!!$      IF(ALLOCATED(ave_gt1))DEALLOCATE(ave_gt1)
!!$      IF(ALLOCATED(ave_gt3))DEALLOCATE(ave_gt3)
!!$      IF(ALLOCATED(ave_gu1))DEALLOCATE(ave_gu1)
!!$      IF(ALLOCATED(ave_gu3))DEALLOCATE(ave_gu3)
!!$      IF(ALLOCATED(ave_gu33))DEALLOCATE(ave_gu33)
!!$      IF(ALLOCATED(ave_gu3gt1))DEALLOCATE(ave_gu3gt1)
!!$      IF(ALLOCATED(ave_gu3gt3))DEALLOCATE(ave_gu3gt3)
!!$      IF(ALLOCATED(ave_gu33gt1))DEALLOCATE(ave_gu33gt1)
!!$      IF(ALLOCATED(ave_gu33gt3))DEALLOCATE(ave_gu33gt3)
!!$      IF(ALLOCATED(ave_gninv))DEALLOCATE(ave_gninv)
!!$      IF(ALLOCATED(ave_gp1inv))DEALLOCATE(ave_gp1inv)
!!$      IF(ALLOCATED(ave_gp3inv))DEALLOCATE(ave_gp3inv)
!!$      IF(ALLOCATED(ave_gradgp1))DEALLOCATE(ave_gradgp1)
!!$      IF(ALLOCATED(ave_gradgr11))DEALLOCATE(ave_gradgr11)
!!$      IF(ALLOCATED(ave_gradgr13))DEALLOCATE(ave_gradgr13)
!!$      IF(ALLOCATED(ave_gradgp1p1))DEALLOCATE(ave_gradgp1p1)
!!$      IF(ALLOCATED(ave_gradgr11p1))DEALLOCATE(ave_gradgr11p1)
!!$      IF(ALLOCATED(ave_gradgr13p1))DEALLOCATE(ave_gradgr13p1)
!!$      IF(ALLOCATED(ave_gradgu1))DEALLOCATE(ave_gradgu1)
!!$      IF(ALLOCATED(ave_gradgu3))DEALLOCATE(ave_gradgu3)
!!$      IF(ALLOCATED(ave_gnp0))DEALLOCATE(ave_gnp0)
!!$      IF(ALLOCATED(ave_gp1p0))DEALLOCATE(ave_gp1p0)
!!$      IF(ALLOCATED(ave_gp3p0))DEALLOCATE(ave_gp3p0)
!!$      IF(ALLOCATED(ave_gr11p0))DEALLOCATE(ave_gr11p0)
!!$      IF(ALLOCATED(ave_gr13p0))DEALLOCATE(ave_gr13p0)
!!$      IF(ALLOCATED(ave_gr33p0))DEALLOCATE(ave_gr33p0)
!!$      IF(ALLOCATED(ave_gnb0))DEALLOCATE(ave_gnb0)
!!$      IF(ALLOCATED(ave_gp1b0))DEALLOCATE(ave_gp1b0)
!!$      IF(ALLOCATED(ave_gp3b0))DEALLOCATE(ave_gp3b0)
!!$      IF(ALLOCATED(ave_gr11b0))DEALLOCATE(ave_gr11b0)
!!$      IF(ALLOCATED(ave_gr13b0))DEALLOCATE(ave_gr13b0)
!!$      IF(ALLOCATED(ave_gr33b0))DEALLOCATE(ave_gr33b0)
!!$      IF(ALLOCATED(ave_gw113b0))DEALLOCATE(ave_gw113b0)
!!$      IF(ALLOCATED(ave_gw133b0))DEALLOCATE(ave_gw133b0)
!!$      IF(ALLOCATED(ave_gw333b0))DEALLOCATE(ave_gw333b0)
!!$      IF(ALLOCATED(ave_gnbp))DEALLOCATE(ave_gnbp)
!!$      IF(ALLOCATED(ave_gp1bp))DEALLOCATE(ave_gp1bp)
!!$      IF(ALLOCATED(ave_gp3bp))DEALLOCATE(ave_gp3bp)
!!$      IF(ALLOCATED(ave_gr11bp))DEALLOCATE(ave_gr11bp)
!!$      IF(ALLOCATED(ave_gr13bp))DEALLOCATE(ave_gr13bp)
!!$      IF(ALLOCATED(ave_gr33bp))DEALLOCATE(ave_gr33bp)
!!$      IF(ALLOCATED(ave_gw113bp))DEALLOCATE(ave_gw113bp)
!!$      IF(ALLOCATED(ave_gw133bp))DEALLOCATE(ave_gw133bp)
!!$      IF(ALLOCATED(ave_gw333bp))DEALLOCATE(ave_gw333bp)
!!$      IF(ALLOCATED(ave_gradgp1p0))DEALLOCATE(ave_gradgp1p0)
!!$      IF(ALLOCATED(ave_gradgr11p0))DEALLOCATE(ave_gradgr11p0)
!!$      IF(ALLOCATED(ave_gradgr13p0))DEALLOCATE(ave_gradgr13p0)
!!$      IF(ALLOCATED(ave_wdgu3gt1))DEALLOCATE(ave_wdgu3gt1)
!!$      IF(ALLOCATED(ave_wdgu3gt3))DEALLOCATE(ave_wdgu3gt3)
!!$      IF(ALLOCATED(ave_wdgu33gt1))DEALLOCATE(ave_wdgu33gt1)
!!$      IF(ALLOCATED(ave_wdgu33gt3))DEALLOCATE(ave_wdgu33gt3)
!!$      IF(ALLOCATED(ave_modwdgu3))DEALLOCATE(ave_modwdgu3)
!!$      IF(ALLOCATED(ave_modwdgu33))DEALLOCATE(ave_modwdgu33)
!!$      IF(ALLOCATED(ave_modwdgu3gt1))DEALLOCATE(ave_modwdgu3gt1)
!!$      IF(ALLOCATED(ave_modwdgu3gt3))DEALLOCATE(ave_modwdgu3gt3)
!!$      IF(ALLOCATED(ave_modwdgu33gt1))DEALLOCATE(ave_modwdgu33gt1)
!!$      IF(ALLOCATED(ave_modwdgu33gt3))DEALLOCATE(ave_modwdgu33gt3)
!!$      IF(ALLOCATED(ave_c_tor_par_gp1p0))DEALLOCATE(ave_c_tor_par_gp1p0)
!!$      IF(ALLOCATED(ave_c_tor_par_gr11p0))DEALLOCATE(ave_c_tor_par_gr11p0)
!!$      IF(ALLOCATED(ave_c_tor_par_gr13p0))DEALLOCATE(ave_c_tor_par_gr13p0)
!!$!      write(*,*)"deallocate ave_g"
!!$! ave_wd_h
!!$      IF(ALLOCATED(ave_wdhp1p0))DEALLOCATE(ave_wdhp1p0)
!!$      IF(ALLOCATED(ave_wdhr11p0))DEALLOCATE(ave_wdhr11p0)
!!$      IF(ALLOCATED(ave_wdhr13p0))DEALLOCATE(ave_wdhr13p0)
!!$      IF(ALLOCATED(ave_wdhp1b0))DEALLOCATE(ave_wdhp1b0)
!!$      IF(ALLOCATED(ave_wdhr11b0))DEALLOCATE(ave_wdhr11b0)
!!$      IF(ALLOCATED(ave_wdhr13b0))DEALLOCATE(ave_wdhr13b0)
!!$      IF(ALLOCATED(ave_wdhp1bp))DEALLOCATE(ave_wdhp1bp)
!!$      IF(ALLOCATED(ave_wdhr11bp))DEALLOCATE(ave_wdhr11bp)
!!$      IF(ALLOCATED(ave_wdhr13bp))DEALLOCATE(ave_wdhr13bp)
!!$      IF(ALLOCATED(ave_wdhu1))DEALLOCATE(ave_wdhu1)
!!$      IF(ALLOCATED(ave_wdhu3))DEALLOCATE(ave_wdhu3)
!!$      IF(ALLOCATED(ave_wdhu33))DEALLOCATE(ave_wdhu33)
!!$      IF(ALLOCATED(ave_modwdhu1))DEALLOCATE(ave_modwdhu1)
!!$      IF(ALLOCATED(ave_wdht1))DEALLOCATE(ave_wdht1)
!!$      IF(ALLOCATED(ave_wdht3))DEALLOCATE(ave_wdht3)
!!$      IF(ALLOCATED(ave_modwdht1))DEALLOCATE(ave_modwdht1)
!!$      IF(ALLOCATED(ave_modwdht3))DEALLOCATE(ave_modwdht3)
!!$!      write(*,*)"deallocated wd_h"
!!$! ave_wd_g
!!$      IF(ALLOCATED(ave_wdgp1p0))DEALLOCATE(ave_wdgp1p0)
!!$      IF(ALLOCATED(ave_wdgr11p0))DEALLOCATE(ave_wdgr11p0)
!!$      IF(ALLOCATED(ave_wdgr13p0))DEALLOCATE(ave_wdgr13p0)
!!$      IF(ALLOCATED(ave_wdgp1b0))DEALLOCATE(ave_wdgp1b0)
!!$      IF(ALLOCATED(ave_wdgr11b0))DEALLOCATE(ave_wdgr11b0)
!!$      IF(ALLOCATED(ave_wdgr13b0))DEALLOCATE(ave_wdgr13b0)
!!$      IF(ALLOCATED(ave_wdgp1bp))DEALLOCATE(ave_wdgp1bp)
!!$      IF(ALLOCATED(ave_wdgr11bp))DEALLOCATE(ave_wdgr11bp)
!!$      IF(ALLOCATED(ave_wdgr13bp))DEALLOCATE(ave_wdgr13bp)
!!$      IF(ALLOCATED(ave_wdgu1))DEALLOCATE(ave_wdgu1)
!!$      IF(ALLOCATED(ave_wdgu3))DEALLOCATE(ave_wdgu3)
!!$      IF(ALLOCATED(ave_wdgu33))DEALLOCATE(ave_wdgu33)
!!$      IF(ALLOCATED(ave_modwdgu1))DEALLOCATE(ave_modwdgu1)
!!$      IF(ALLOCATED(ave_wdgt1))DEALLOCATE(ave_wdgt1)
!!$      IF(ALLOCATED(ave_wdgt3))DEALLOCATE(ave_wdgt3)
!!$      IF(ALLOCATED(ave_modwdgt1))DEALLOCATE(ave_modwdgt1)
!!$      IF(ALLOCATED(ave_modwdgt3))DEALLOCATE(ave_modwdgt3)
!!$!      write(*,*)"deallocate wd_g"
!!$! kpar_h
!!$      IF(ALLOCATED(ave_kparhnp0))DEALLOCATE(ave_kparhnp0)
!!$      IF(ALLOCATED(ave_kparhp1p0))DEALLOCATE(ave_kparhp1p0)
!!$      IF(ALLOCATED(ave_kparhp3p0))DEALLOCATE(ave_kparhp3p0)
!!$      IF(ALLOCATED(ave_kparhp1b0))DEALLOCATE(ave_kparhp1b0)
!!$      IF(ALLOCATED(ave_kparhr11b0))DEALLOCATE(ave_kparhr11b0)
!!$      IF(ALLOCATED(ave_kparhr13b0))DEALLOCATE(ave_kparhr13b0)
!!$      IF(ALLOCATED(ave_kparhnbp))DEALLOCATE(ave_kparhnbp)
!!$      IF(ALLOCATED(ave_kparhp3bp))DEALLOCATE(ave_kparhp3bp)
!!$      IF(ALLOCATED(ave_kparhp1bp))DEALLOCATE(ave_kparhp1bp)
!!$      IF(ALLOCATED(ave_kparhr11bp))DEALLOCATE(ave_kparhr11bp)
!!$      IF(ALLOCATED(ave_kparhr13bp))DEALLOCATE(ave_kparhr13bp)
!!$      IF(ALLOCATED(ave_kparhu1))DEALLOCATE(ave_kparhu1)
!!$      IF(ALLOCATED(ave_kparhu3))DEALLOCATE(ave_kparhu3)
!!$      IF(ALLOCATED(ave_kparht1))DEALLOCATE(ave_kparht1)
!!$      IF(ALLOCATED(ave_kparht3))DEALLOCATE(ave_kparht3)
!!$      IF(ALLOCATED(ave_modkparhu1))DEALLOCATE(ave_modkparhu1)
!!$      IF(ALLOCATED(ave_modkparhu3))DEALLOCATE(ave_modkparhu3)
!!$!      write(*,*)"deallocated kpar_h"
!!$! kpar_g
!!$      IF(ALLOCATED(ave_kpargnp0))DEALLOCATE(ave_kpargnp0)
!!$      IF(ALLOCATED(ave_kpargp1p0))DEALLOCATE(ave_kpargp1p0)
!!$      IF(ALLOCATED(ave_kpargp3p0))DEALLOCATE(ave_kpargp3p0)
!!$      IF(ALLOCATED(ave_kpargp1b0))DEALLOCATE(ave_kpargp1b0)
!!$      IF(ALLOCATED(ave_kpargr11b0))DEALLOCATE(ave_kpargr11b0)
!!$      IF(ALLOCATED(ave_kpargr13b0))DEALLOCATE(ave_kpargr13b0)
!!$      IF(ALLOCATED(ave_kpargnbp))DEALLOCATE(ave_kpargnbp)
!!$      IF(ALLOCATED(ave_kpargp3bp))DEALLOCATE(ave_kpargp3bp)
!!$      IF(ALLOCATED(ave_kpargp1bp))DEALLOCATE(ave_kpargp1bp)
!!$      IF(ALLOCATED(ave_kpargr11bp))DEALLOCATE(ave_kpargr11bp)
!!$      IF(ALLOCATED(ave_kpargr13bp))DEALLOCATE(ave_kpargr13bp)
!!$      IF(ALLOCATED(ave_kpargu1))DEALLOCATE(ave_kpargu1)
!!$      IF(ALLOCATED(ave_kpargu3))DEALLOCATE(ave_kpargu3)
!!$      IF(ALLOCATED(ave_kpargt1))DEALLOCATE(ave_kpargt1)
!!$      IF(ALLOCATED(ave_kpargt3))DEALLOCATE(ave_kpargt3)
!!$      IF(ALLOCATED(ave_modkpargu1))DEALLOCATE(ave_modkpargu1)
!!$      IF(ALLOCATED(ave_modkpargu3))DEALLOCATE(ave_modkpargu3)
!!$!      write(*,*)"deallocate kpar_g"
!!$! gradB_h
!!$      IF(ALLOCATED(ave_gradBhp1))DEALLOCATE(ave_gradBhp1)
!!$      IF(ALLOCATED(ave_gradBhp3))DEALLOCATE(ave_gradBhp3)
!!$      IF(ALLOCATED(ave_gradBhr11))DEALLOCATE(ave_gradBhr11)
!!$      IF(ALLOCATED(ave_gradBhr13))DEALLOCATE(ave_gradBhr13)
!!$      IF(ALLOCATED(ave_gradBhr33))DEALLOCATE(ave_gradBhr33)
!!$      IF(ALLOCATED(ave_gradBhu1))DEALLOCATE(ave_gradBhu1)
!!$      IF(ALLOCATED(ave_gradBhu3))DEALLOCATE(ave_gradBhu3)
!!$      IF(ALLOCATED(ave_gradBhu33))DEALLOCATE(ave_gradBhu33)
!!$!      write(*,*)"deallocate gradB_h"
!!$! gradB_g
!!$      IF(ALLOCATED(ave_gradBgp1))DEALLOCATE(ave_gradBgp1)
!!$      IF(ALLOCATED(ave_gradBgp3))DEALLOCATE(ave_gradBgp3)
!!$      IF(ALLOCATED(ave_gradBgr11))DEALLOCATE(ave_gradBgr11)
!!$      IF(ALLOCATED(ave_gradBgr13))DEALLOCATE(ave_gradBgr13)
!!$      IF(ALLOCATED(ave_gradBgr33))DEALLOCATE(ave_gradBgr33)
!!$      IF(ALLOCATED(ave_gradBgu1))DEALLOCATE(ave_gradBgu1)
!!$      IF(ALLOCATED(ave_gradBgu3))DEALLOCATE(ave_gradBgu3)
!!$      IF(ALLOCATED(ave_gradBgu33))DEALLOCATE(ave_gradBgu33)
!!$!      write(*,*)"deallocated gradB_g"
!!$!  ave_theta
!!$!
!!$      IF(ALLOCATED(ave_kx))DEALLOCATE(ave_kx)
!!$      IF(ALLOCATED(ave_c_tor_par))DEALLOCATE(ave_c_tor_par)
!!$      IF(ALLOCATED(ave_c_tor_per))DEALLOCATE(ave_c_tor_per)
!!$      IF(ALLOCATED(ave_c_par_par))DEALLOCATE(ave_c_par_par)
!!$      IF(ALLOCATED(ave_wdh))DEALLOCATE(ave_wdh)
!!$      IF(ALLOCATED(ave_wdg))DEALLOCATE(ave_wdg)
!!$      IF(ALLOCATED(ave_modwdh))DEALLOCATE(ave_modwdh)
!!$      IF(ALLOCATED(ave_modwdg))DEALLOCATE(ave_modwdg)
!!$      IF(ALLOCATED(ave_gradB))DEALLOCATE(ave_gradB)
!!$      IF(ALLOCATED(ave_lnB))DEALLOCATE(ave_lnB)
!!$      IF(ALLOCATED(ave_b0))DEALLOCATE(ave_b0)
!!$      IF(ALLOCATED(ave_b0inv))DEALLOCATE(ave_b0inv)
!!$      IF(ALLOCATED(ave_kpar))DEALLOCATE(ave_kpar)
!!$      IF(ALLOCATED(ave_modkpar))DEALLOCATE(ave_modkpar)
!!$      IF(ALLOCATED(ave_p0))DEALLOCATE(ave_p0)
!!$      IF(ALLOCATED(ave_p0inv))DEALLOCATE(ave_p0inv)
!!$      IF(ALLOCATED(ave_bp))DEALLOCATE(ave_bp)
!!$      IF(ALLOCATED(ave_bpinv))DEALLOCATE(ave_bpinv)
!!$      IF(ALLOCATED(ave_kpar_eff))DEALLOCATE(ave_kpar_eff)
!!$      IF(ALLOCATED(ave_modkpar_eff))DEALLOCATE(ave_modkpar_eff)
!!$!      write(*,*)"deallocated ave_theta"
!!$!
      END SUBROUTINE tglf_harvest
!-----------------------------------------------------
!
