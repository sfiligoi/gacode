     MODULE nub4
     USE param,only : kj,kbs,ke,kb
     implicit none

!           new common block related specifically to 
!           time dependent beam module
!           HSJ  4/22/00
	integer ns,nbe,npulse(kb),n21s,                              &
               n_pulse,kj1,kbs1,ke1,kb1,kt1,                         &
               method,time_dep_beam,                                 &
               beam_thermal_cutoff,beam_init_restart_file,           &
               beam_time_init,beam_pulse_control,kt           
!	parameter (kion = 3,kj=51,kb=2,ke=3,kbs =2,kt = 200)
        parameter (kt = 200)
        logical no_fast_losses
        logical *1 pssv(kt,kj,kbs,ke,kb), pssvoff(kt,kj,kbs,ke,kb)
        character beam_restart_file*256, beam_mode*8
	real *8                                                     &
             tau0_vlj(kt,kj,kbs,ke,kb),time_start,                  &
             tau0_vljoff(kt,kj,kbs,ke,kb),mass_beami(kb),           &
             tau0(kj),enbeam_tot(kj),tauslow,bstime,time0_beam,     &
             wenbeam_part(kj,kbs,ke,kb),wenbeam_tot(kj),            &
             Qfi_part(kj,kbs,ke,kb),Qfe_tot(kj),Qfi_tot(kj),        &
             enbeam_part_nl(kj,kbs,ke,kb),enbeam_tot_nl(kj),        &
             pitch_angle,Rfi_part(kj,kbs,ke,kb),                    &
             Rfe_part(kj,kbs,ke,kb),Rfi_tot(kj),Rfe_tot(kj),        &
             beam_intensity(kt,kj,kbs,ke,kb),pbeamOn(kt,kbs,kb),    &
             pbeamOff(kt,kbs,kb), source2_phase(kb),                &
             enbeam_part(kj,kbs,ke,kb),vcrit(kj,kb),vbeam(ke,kb),   &
             vthi(kj),beam_thermal_speed,therm_frac,                &
             nf_tot_loss(kj),nf_tot_therm(kj),d_nf_tot_dt(kj),      &
             enbeam_tot_prev(kj),nf_tot_source(kj),nf_conf_time,    &
             enbeam_tot_prev_nl(kj),pwf_tot_source(kj)

   

!            n21s =1    !restricts the source index to 1 instead of the
                        !usual 2. This is done because the freya results
                        ! are returned with both sources mixed into a
                        !single array index. (this is OK so long as
                        !the two sources act at the same time, which is
                        !what is currently assumed here )
      data      beam_restart_file /'null'/

     END MODULE nub4
