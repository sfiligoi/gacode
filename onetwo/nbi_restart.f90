 MODULE nbi_restart
!
!Module sets up the initial thermal ion densities for the
!nubeam restart runs (eg if nubeam_restart =1).
!
!nubeam_restart = +1  
!nubeam_restart = -1  initial guess of beam density 
!                     from stand-alone pre-run of nubeam 
!                     during nubeam_back_delt
!                     with fixed plasma profiles at time0 JMP2005.12.10
!                    
 USE param,     ONLY : kk,kj
 USE fusion,    ONLY : beam_thermal_dtntot, beam_thermal_ddntot,       &
                       beam_thermal_ddptot, beam_thermal_tt2ntot,      &
                       beam_beam_dtntot, beam_beam_ddntot,             &
                       beam_beam_tt2ntot,beam_beam_ddptot,             &
                       beam_beamddn,beam_beamdtn,beam_beamddp,         &
                       beam_beamtt2n,beam_thermaltth_df,               &
                       beam_thermalddp, beam_thermalddn,               &
                       beam_thermaltt2n,beam_thermaldth_tf 

 USE transp,    ONLY :                                                 &
                    nubeam_restart,nubeam0_dt,nubeam_back_delt,        & !JMP
                    nubeam_calls,qbeame_nub,qbeami_nub,                &
                    qbeame_nubp,qbeami_nubp,nubeam_steps,              &
                    qbth_nubp,qbth_nub,                                &
                    qbeame_intg_nubp,qbeame_intg_nub,                  &
                    qbeami_intg_nubp,qbeami_intg_nub,                  &
                    qbth_intg_nubp,qbth_intg_nub,                      &
                    curb_nub,curb_nubp,curb_intg_nub,                  &
                    curb_intg_nubp,storqueb_nub,                       &
                    storqueb_nubp,storqueb_intg_nub,                   &
                    storqueb_intg_nubp,sprbeame_nub,                   &
                    sprbeame_nubp,sprbeame_intg_nub,                   &
                    sprbeame_intg_nubp,sprbeami_nub,                   &
                    sprbeami_nubp,sprbeami_intg_nub,                   &
                    sprbeami_intg_nubp,beam_data,                      &
                    wbeam_nub,wbeam_nubp,wbeam_intg_nub,               &
                    enbeam_intg_nubp,wbeam_intg_nubp,                  &
                    enbeam_nub,enbeam_nubp,enbeam_intg_nub,            &
                    enbeam_intg_nubp,pwf_tot_intg_nub,                 &
                    pwf_tot_intg_nubp,                                 &
                    beam_thermal_dtntot_nub,beam_thermal_dtntot_nubp,  &
                    beam_thermal_ddntot_nub,beam_thermal_ddntot_nubp,  &
                    beam_thermal_ddptot_nub,beam_thermal_ddptot_nubp,  &
                    beam_thermal_tt2ntot_nub,beam_thermal_tt2ntot_nubp,&
                    beam_beam_dtntot_nub,beam_beam_dtntot_nubp,        &
                    beam_beam_ddntot_nub,beam_beam_ddntot_nubp,        &
                    beam_beam_ddptot_nub,beam_beam_ddptot_nubp,        &
                    beam_beam_tt2ntot_nub,beam_beam_tt2ntot_nubp,      &
                    beam_beamddp_nub,beam_beamddp_nubp,                &
                    beam_beamtt2n_nub,beam_beamtt2n_nubp,              &
                    beam_beamddn_nub,beam_beamddn_nubp,                &
                    beam_beamdtn_nub,beam_beamdtn_nubp,                &
                    beam_thermaltth_df_nub,beam_thermaltth_df_nubp,    &
                    beam_thermalddp_nub,beam_thermalddp_nubp,          &
                    beam_thermalddn_nub,beam_thermalddn_nubp,          &
                    beam_thermaltt2n_nub,beam_thermaltt2n_nubp,        &
                    beam_thermaldth_tf_nub,beam_thermaldth_tf_nubp,    &
                    beam_data,nubeam_evolve,ifix_nubeam_dt !JMP
                    

 USE nub,       ONLY : nbeams
 USE nub2,      ONLY : enbeam,enbs,enbsav,enb,wbeam,                   &
                       wb,enbmin !JMP
 USE numbrs,    ONLY : nj,nk,nion
 USE verbose,   ONLY : tportvb
 USE soln,      ONLY : u,usave,te,ti,ene,enesav,en,rbp
 USE tordlrot,  ONLY : iangrot,angrot, storqueb,storqueb_intg,         &
                       sprbeame,sprbeami,sprbeame_intg,                &
                       sprbeami_intg
 USE sourc,     ONLY : qbeame,qbeami,qbeame_intg,qbeami_intg,          &
                       qbth_intg,qbth,curb_intg,curb,wbeam_intg,       &
                       enbeam_intg,curbi

 USE solcon,    ONLY : time,time0 !JMP
 USE mhdcom,    ONLY : mhdmethd !JMP
 USE ename     ,ONLY : eqdsk_tdem !JMP
 
 IMPLICIT NONE

 CONTAINS

   SUBROUTINE beam_prof_init(time_init)

     REAL *8,intent(in) :: time_init
     INTEGER err_state,err_xplasma,pload_12
     INTEGER j
     INTEGER nubeam_iter,nubeam_iter_tot
     REAL *8 nubeam_last_dt,nubeam0_dt_temp,time0_temp,time0_input
     REAL *8 time_eqdsk !JMP
     INTEGER ifix_nubeam_dt_save !JMP
     
     !get beam density from restart file and use it to
     !get consistent thermal ion densitites

     if (nubeam_restart .eq. 1) then !JMP

       nubeam_evolve = -1 !prevents enbeam =0.0 in sub source 

       !check that the three required restart files are readable:
       !(eq xplasma and nubeam state cdf files and profile file for onetwo restart)
       CALL check_restart(err_state,err_xplasma)
       IF(err_state + err_xplasma .GT. 0)THEN
         PRINT *,'check_restart reports :'
         PRINT *,'err_state,err_xplasma + ', err_state,err_xplasma
         CALL STOP('beam_prof_init',1)
       ENDIF

       CALL read_restart_profs   

       !load 12 variables  
       pload_12 = 1  
       CALL get_12_fiprof(time_init,pload_12)

       enb(:,:,:) = 0.0
       enbsav(:,:,:) = 0.0
       enb(:,1,1) = enbeam(:)   !define only the first componenet
       enbsav(:,1,1) = enbeam(:)
       enbs(1:nj) = enbeam(1:nj)

       !redate copies te,etc into u:
       call redate (u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
       call savesol   (u, usave, nk, nj, kk) ! sets usave =u 
       call update (u, en, te, ti, rbp, nk, nj, kj, kk,iangrot,angrot)
       call copya (ene, enesav, nj)

     elseif (nubeam_restart .eq. -1) then !JMP BLOCK START

       print *,'*** ' 
       print *,'*** calling nubeam for initial guess of beam density'
       print *,'*** ' 

       ifix_nubeam_dt_save = ifix_nubeam_dt
       ifix_nubeam_dt = 1
       
       call set_nubeam_init_profs 

       nubeam_iter_tot = nubeam_back_delt/nubeam0_dt
       nubeam_last_dt = nubeam_back_delt-nubeam_iter_tot*nubeam0_dt

       !hack the starting time 
       !this is consistent with the boundary conditions and 
       !will not cause any side-effects and 
       
       if(eqdsk_tdem .eq. 'tdem' ) then
         time_eqdsk = time0
         call wrt_tdem_eqdsk(time_eqdsk,'g0.restart') 
         print *,'g0.restart written at',time0
       end if  

       time0_temp=time
       time = time0-nubeam_back_delt
       
       do nubeam_iter = 1,nubeam_iter_tot+2 

         print *,'***' 
         print *,'***',nubeam_iter,'th initial nubeam call',nubeam_iter_tot,time0
         print *,'***' 

         if (nubeam_iter .gt. 1) call update_nubeam_profs

         nubeam0_dt_temp = nubeam0_dt
         if (nubeam_iter .eq. nubeam_iter_tot+1) then
           if (nubeam_last_dt .gt. 0.0) then
             nubeam0_dt = nubeam_last_dt
           else
             cycle
           endif
         endif 

         call zen
         call NTCC_driver

         pload_12 = 0  
         call get_12_fiprof(time,pload_12)
         !wbeam,enbeam are set in nubeam ( for time = time+dt )
         wb(:,1,1) = wbeam(:)/0.62415064e16  !only define for (Kj:1,1)
         enb(:,1,1) = enbeam(:)
         do j=1,SIZE( enb,dim=1)
           enb(j,1,1) = MAX(enb(j,1,1),enbmin) !zero doesnt work
         enddo

         time = time+nubeam0_dt
         nubeam0_dt = nubeam0_dt_temp

       end do
       
       call redate (u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
       call savesol (u, usave, nk, nj, kk)

       time = time0_temp
       ifix_nubeam_dt = ifix_nubeam_dt_save

     endif !JMP BLOCK END

   END SUBROUTINE beam_prof_init

   SUBROUTINE set_fiprof
      INTEGER j,ksymp1,ksymp2
      ksymp1 = ((3*nbeams+1)*3*nbeams)/2+1      !for beam_beam
      ksymp2  = 3*nbeams+1                      !for beam_thermal
            qbeame(1:nj)    = qbeame_nubp(1:nj)
            qbeami(1:nj)    = qbeami_nubp(1:nj)
            qbth(1:nj)      = qbth_nubp(1:nj)
            storqueb(1:nj)  = storqueb_nubp(1:nj)
            curbi(1:nj)      = curb_nubp(1:nj)
            storqueb(1:nj)  = storqueb_nubp(1:nj)
            sprbeame(1:nj)  = sprbeame_nubp(1:nj)
            sprbeami(1:nj)  = sprbeami_nubp(1:nj)
            wbeam(1:nj)     = wbeam_nubp(1:nj)
            enbeam(1:nj)    = enbeam_nubp(1:nj)
            beam_beamdtn(1:nj,ksymp1) = beam_beamdtn_nubp(1:nj)
            beam_beamddp(1:nj,ksymp1) = beam_beamddp_nubp(1:nj)
            beam_beamddn(1:nj,ksymp1) = beam_beamddn_nubp(1:nj) 
            beam_beamtt2n(1:nj,ksymp1) = beam_beamtt2n_nubp(1:nj) 
            beam_thermaltth_df(1:nj,ksymp1) = beam_thermaltth_df_nubp(1:nj)
            beam_thermalddp(1:nj,ksymp1) = beam_thermalddp_nubp(1:nj)
            beam_thermalddn(1:nj,ksymp1) = beam_thermalddn_nubp(1:nj)
            beam_thermaltt2n(1:nj,ksymp1) =  beam_thermaltt2n_nubp(1:nj)
            beam_thermaldth_tf(1:nj,ksymp1) = beam_thermaldth_tf_nubp(1:nj)

            qbeame_intg     = qbeame_intg_nubp
            qbeami_intg     = qbeami_intg_nubp
            qbth_intg       = qbth_intg_nubp
            curb_intg       = curb_intg_nubp
            storqueb_intg   = storqueb_intg_nubp
            sprbeame_intg   = sprbeame_intg_nubp
            sprbeami_intg   = sprbeami_intg_nubp
            beam_data%pwf_tot_intg = pwf_tot_intg_nubp
            wbeam_intg      = wbeam_intg_nubp
            enbeam_intg     = enbeam_intg_nubp
            beam_thermal_dtntot  = beam_thermal_dtntot_nubp
            beam_thermal_ddntot  = beam_thermal_ddntot_nubp
            beam_thermal_ddptot  = beam_thermal_ddptot_nubp 
            beam_thermal_tt2ntot = beam_thermal_tt2ntot_nubp 
            beam_beam_dtntot     = beam_beam_dtntot_nubp
            beam_beam_ddntot     = beam_beam_ddntot_nubp
            beam_beam_ddptot     = beam_beam_ddptot_nubp 
            beam_beam_tt2ntot    = beam_beam_tt2ntot_nubp 

   END SUBROUTINE set_fiprof

 END MODULE nbi_restart
