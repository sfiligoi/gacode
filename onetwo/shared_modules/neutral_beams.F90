 MODULE neutral_beams
    USE nrtype ,                             ONLY : DP, I4B,I2B
    USE ions_gcnmp,                          ONLY : name_size
!#ifndef GCNMP
    ! use of kz here is ok for Onetwo.  In P_Nfreya the kz dimensioned arrays
    ! are replaced with allocatable arrays in zonal_data.f90
    USE nf_param,                            ONLY : kb,nap,kj,kz,ke,kion,kcm,kf, &
                                                    kprim,kimp,kbs
 
!#endif
    IMPLICIT NONE 
    INTEGER(I4B), PARAMETER                   :: maxchr = 256  
    INTEGER(I4B), PARAMETER                   :: shortchr = 16    
    INTEGER(I4B), PARAMETER                   :: nb_ec = 3       ! # energy components per beamline
    INTEGER,PARAMETER,PUBLIC                  :: nbeam_ml= 4     ! # elemens in nameb_ml
    INTEGER,PARAMETER                         :: npitch = 25     ! # allowed pitch angles
                                                                 ! in orbloss , not implemented
    REAL(DP)  nb_net_torque                                      ! total beam torque *sign Ip

             
    INTEGER(I4B)  nbeams                                         ! # active neutral beams
    INTEGER(I4B)  no_injectors                                   ! total,physical plus pseduo
                                                                 ! in namelist this is number
                                                                 ! of physical injectors and is a
                                                                 ! pseudonym for nobeams
    INTEGER(I4B)  no_physical_injectors                          ! code may assign additional copies of injectors
    REAL(DP), SAVE,ALLOCATABLE,DIMENSION(:,:) :: enbeam          ! beam density by species (sumed over 3 energy componenets)
    REAL(DP), SAVE,ALLOCATABLE,DIMENSION(:)   :: wbeam           ! beam ion stored energy density
    REAL(DP), SAVE,ALLOCATABLE,DIMENSION(:)   :: storqueb,enbeam_tot
    REAL(DP), SAVE,DIMENSION(nb_ec)           :: nb_energy
    INTEGER(I4B),PUBLIC                       :: nbion           ! # beam ions
    REAL(DP) fd_beam                                             ! fraction of d in beam

    REAL(DP)   beam_sim_time_start,beam_sim_time_end

    LOGICAL no_fast_losses, split_injectors

    CHARACTER(len = name_size) ,PUBLIC, DIMENSION(:), POINTER :: nameb

!#ifndef GCNMP
    REAL(DP) sfrac1(kb), bvofset(kb),bhofset(kb), fbcur(ke, kb)

    CHARACTER(len = name_size) ,PUBLIC, DIMENSION(nbeam_ml)   :: nameb_ml 
    INTEGER(I4B),DIMENSION(nbeam_ml),PUBLIC                   :: nameb_index
    LOGICAL use_ufile,calc_variance ,nfreya_vb
    LOGICAL iterate_beam,beam_iteration,randomize_seed
    CHARACTER(len = maxchr) write_performance_data  

     INTEGER(I4B), PARAMETER ::  kt = 200       ! # pulses ,verify??

     INTEGER(I4B), PARAMETER :: ntimbplt = 15


     INTEGER(I4B), ALLOCATABLE,DIMENSION(:)    ::  pseudo_injectors,         &
                                                   sorted_physical_injectors,&
                                                   injctr_wk_list,           &
                                                   phys_injctr_wk_list
     INTEGER(I4B), ALLOCATABLE,DIMENSION(:,:)  ::  npart_pseudo,nmbrz

     INTEGER(I4B)                                                           &
                      iexcit, ilorent, mstate, ncont, kdene,                &
                      kdeni, kdenz, ksvi, ksvz, ksve, krad, ngh,            &
                      ngl, iz(kimp), izstrp(kimp),iatype(nap,kb)


     INTEGER(I4B)     ns,nbe,npulse(kb),n21s,                               &
                      n_pulse,kj1,kbs1,ke1,kb1,kt1,                         &
                      method,time_dep_beam,no_birth_points,                 &
                      beam_thermal_cutoff,beam_init_restart_file,           &
                      beam_time_init,beam_pulse_control  

      INTEGER(I4B)       fast_ion_target, neg_ion_source(kb),master_injector_no

      INTEGER(I4B)                                                         &
             imaxbeamconvg, naptr, ne_tk, ibeam, ibion,               &
             inubplt, npart, npskip, nsourc,namelist_nbion,mi,mj,     &
             mim1,mjm1,nbdep,    & ! nbdep is unit for beamdep output file in sub nbdep2d
             npart_all_beamlines   ! total number of pseudo particles all beam lines combined
                                   ! individual beamslines get their fraction from this
                                   ! according to beam power, note that npart_all_beamlines
                                   ! is a pseudonym for npart 


     INTEGER(I4B), DIMENSION(2) :: npat

     INTEGER(I4B)  norb,nrpat,nzpat

     INTEGER(I4B), DIMENSION(5) :: idebug

     INTEGER*4,PUBLIC    ::  itrapfi, itrapech, ibcur, ibcx,       &
                             iborb, ibslow,inubpat, n_izpt
        
     INTEGER(I2B) ncorin,nfreya_plot_unit

     CHARACTER(maxchr) nfreya_plot_file_name

     CHARACTER (shortchr) inj_list ! indicates if master did none,all or 1 injector
 


     LOGICAL *1 pssv(kt,kj,kbs,ke,kb), pssvoff(kt,kj,kbs,ke,kb)

 

     CHARACTER beam_restart_file*256, beam_mode*8

     CHARACTER(len = name_size)  nbshape(kb), nashape(nap,kb)
 
     CHARACTER(len = name_size) ,PUBLIC, DIMENSION(:), POINTER :: namelist_nameb  



      REAL(DP)     drpat,dzpat,elong,bntot,time_now
 
      REAL (DP),PUBLIC,SAVE ::                                        &
             anglev(kb), angleh(kb),                                  &
             aheigh(nap, kb), awidth(nap, kb), bcur(kb),              &
             alen(nap, kb), bleni(kb), blenp(kb),        &
             bptor(kb), bheigh(kb),                      &
             bwidth(kb), bhfoc(kb), bvfoc(kb), bhdiv(kb),             &
             bvdiv(kb), beam_on(kb), beam_end(kb), beam_time(kb),     &
             beam_off(kb),beam_cycles(kb),                            &
             ebkev(kb), ennub(kj),                     &
             hicm(kj,ke,kb,kcm), ds_tk,                               &
             fe_tk, zeta(kj,ke,kb), csgn, hibrs(kj,kb),               &
             ftrapfit(kj,ke,kb), angmpf(kj,ke,kb),  fdbeam,           &
             psif(kf), rowpsi(kj), relnub,                &
             rpivot(kb), zpivot(kb), tenub(kj),                       &
             timbplt(ntimbplt), sb(kj,ke,kb), sbcx(kj,2), sbion(kj),  &
             spb(kj,ke,kb), qb(kj,ke,kb), qbbe(kj,ke,kb),             &
             qbbi(kj,ke,kb), qbf(kj,ke,kb), zzi(kz,kion+1),           &
             zne(kz), zni(kz,kion), zte(kz), zangrot(kz),             &
             qbeami_rot(kj), qbeame_rot(kj), psivol(kz),              &
             freyr(kf), fionx, hdepsmth, sbpure(kj,ke,kb),            &
             atw_beam, rtstcx, relaxden, relaxden_err,de_tk

      REAL (DP),PUBLIC,SAVE ::                                      &
             ecrit(kj), emzrat(kj), enbs(kj),                       &
             enbav0(kj), enbav1(kj),tbeam(kj),taus(kj),             &
             olossc(kj,ke,kb), oloss(kj,ke,kb),                     &
             wbeamperp(kj),wbeampara(kj),                           &
             nb_thresh ,enbmin ,tmin_curray,enbmin_curray

      REAL (DP)                                                     &
             tau0_vlj(kt,kj,kbs,ke,kb),time_start,                  &
             tau0_vljoff(kt,kj,kbs,ke,kb),mass_beami(kb),           &
             tau0(kj),tauslow,bstime,time0_beam,     &
             wenbeam_part(kj,kbs,ke,kb),wenbeam_tot(kj),            &
             Qfi_part(kj,kbs,ke,kb),Qfe_tot(kj),Qfi_tot(kj),        &
             enbeam_part_nl(kj,kbs,ke,kb),enbeam_tot_nl(kj),        &
             pitch_angle,Rfi_part(kj,kbs,ke,kb),                    &
             Rfe_part(kj,kbs,ke,kb),Rfi_tot(kj),Rfe_tot(kj),        &
             beam_intensity(kt,kj,kbs,ke,kb),pbeamOn(kt,kbs,kb),    &
             pbeamOff(kt,kbs,kb), source2_phase(kb),                &
             enbeam_part(kj,kbs,ke,kb),vcrit(kj,kb),                &
             vthi(kj),beam_thermal_speed,therm_frac,                &
             nf_tot_loss(kj),nf_tot_therm(kj),d_nf_tot_dt(kj),      &
             enbeam_tot_prev(kj),nf_tot_source(kj),nf_conf_time,    &
             enbeam_tot_prev_nl(kj),pwf_tot_source(kj),             &
             thetp(kb),     &
             thetpp(kb),costp(kb),sintp(kb),costpp(kb),sintpp(kb)



     REAL(DP),ALLOCATABLE,DIMENSION(:,:)     :: ebeam,pbeam,vbeam,bion,  &
                                                bneut,forb, fb11,        &
                                                fb10, fb01, fb00, fber,  &  
                                                wb11, wb10, wb01, wb00,  &
                                                fap,fwall


     REAL(DP),ALLOCATABLE,DIMENSION(:,:,:)   :: bencap, fbe, fbi,fbth,   &
                                                bke, bki,enb, enbsav,    &
                                                enbav,hdep,hdepz, pb0,   &
                                                ppb, ppbsav, ppbav,qbsav,&
                                                sbsav,spbsav, taupb,     &
                                                tauppb, taueb, wb, wbsav,&
                                                wbav,zetaz,hibr,hibrz,   &
                                                ftrapfi,angmpz

     REAL(DP),ALLOCATABLE,DIMENSION(:,:,:,:) :: hicmz

     REAL(DP),ALLOCATABLE,DIMENSION(:)       :: cangv,cangh,sangv,sangh, &
                                                rtan_ls,rtan_rs
                                                

     REAL(DP),ALLOCATABLE,DIMENSION(:,:,:)   :: vx_izpt,vy_izpt,vz_izpt, &
                                                pitch_a,x_izpt,y_izpt,   &
                                                z_izpt,r_izpt

     INTEGER(I4B),ALLOCATABLE,DIMENSION(:,:) :: nsample_izpt

      ! no_*d keeps track of the number of items (od to 4d) that will be 
      ! passed by mpi in a buffer (mpi_buffer):
      INTEGER(i4B)   no_0di,no_3dr_p,no_2di,no_2dr,no_3dr,no_4dr,no_2di2

!            n21s =1    !restricts the source index to 1 instead of the
                        !usual 2. This is done because the freya results
                        ! are returned with both sources mixed into a
                        !single array index. (this is OK so long as
                        !the two sources act at the same time, which is
                        !what is currently assumed here )

      DATA      beam_restart_file /'null'/

      DATA      nbdep /0_I2B/

      DATA enbmin / 1.e3_DP /                     !min fast ion density



  DATA nameb_ml /'h', 'd', 't', 'dt'/
  DATA nameb_index /1,2,3,4/

!#endif

    CONTAINS 
      SUBROUTINE neut_beam_allocate
      !-------------------------------------------------------------------------
      ! -- allocate some  arrays associated with neutral beam statefile output:
      !-------------------------------------------------------------------------
        USE nrtype,                    ONLY : DP,I4B,I2B

        USE Plasma_properties ,        ONLY : neut_beam

        USE common_constants,          ONLY : convert,zeroc,izero

 !       USE grid_class,                ONLY : nj

        USE nf_param,                  ONLY : kcm

        IMPLICIT NONE

        INTEGER(I4B),PARAMETER :: ke_bm     = 3         ! # beam energies
  


      IF(neut_beam%nbeams .GT. 0)THEN

  
         IF(.NOT. ASSOCIATED(neut_beam%sb) ) THEN
            ALLOCATE(neut_beam%sb(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%sb(:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%qb))THEN
            ALLOCATE(neut_beam%qb(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%qb(:,:,:) =zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%angmpf))THEN
            ALLOCATE(neut_beam%angmpf(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%angmpf(:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%pb0))THEN
            ALLOCATE(neut_beam%pb0(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%pb0(:,:,:)= zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%spb))THEN
            ALLOCATE(neut_beam%spb(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%spb(:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%spbr))THEN
            ALLOCATE(neut_beam%spbr(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%spbr(:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%ebeam))THEN
            ALLOCATE(neut_beam%ebeam(ke_bm,neut_beam%nbeams))
            neut_beam%ebeam(:,:)= zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%pbeam))THEN
            ALLOCATE(neut_beam%pbeam(ke_bm,neut_beam%nbeams))
            neut_beam%pbeam(:,:)= zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%fber))THEN
            ALLOCATE(neut_beam%fber(ke_bm,neut_beam%nbeams))
            neut_beam%fber(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%fb00))THEN
            ALLOCATE(neut_beam%fb00(ke_bm,neut_beam%nbeams))
            neut_beam%fb00(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%fb01))THEN
            ALLOCATE(neut_beam%fb01(ke_bm,neut_beam%nbeams))
            neut_beam%fb01(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%fb10))THEN
            ALLOCATE(neut_beam%fb10(ke_bm,neut_beam%nbeams))
            neut_beam%fb10(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%fb11))THEN
            ALLOCATE(neut_beam%fb11(ke_bm,neut_beam%nbeams))
            neut_beam%fb11(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%wb00))THEN
            ALLOCATE(neut_beam%wb00(ke_bm,neut_beam%nbeams))
            neut_beam%wb00(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%wb01))THEN
            ALLOCATE(neut_beam%wb01(ke_bm,neut_beam%nbeams))
            neut_beam%wb01(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%wb10))THEN
            ALLOCATE(neut_beam%wb10(ke_bm,neut_beam%nbeams))
            neut_beam%wb10(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%wb11))THEN
            ALLOCATE(neut_beam%wb11(ke_bm,neut_beam%nbeams))
            neut_beam%wb11(:,:) = zeroc
         ENDIF


         IF(.NOT. ASSOCIATED(neut_beam%hicm))THEN
            ALLOCATE(neut_beam%hicm(neut_beam%nj_beam,ke_bm,neut_beam%nbeams,kcm))
            neut_beam%hicm(:,:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%fwall))THEN
            ALLOCATE(neut_beam%fwall(ke_bm,neut_beam%nbeams))
            neut_beam%fwall(:,:) = zeroc
         ENDIF



         IF(.NOT. ASSOCIATED(neut_beam%fap))THEN
            ALLOCATE(neut_beam%fap(ke_bm,neut_beam%nbeams))
            neut_beam%fap(:,:) = zeroc
         ENDIF
   
         IF(.NOT. ASSOCIATED(neut_beam%forb))THEN
            ALLOCATE(neut_beam%forb(ke_bm,neut_beam%nbeams))
            neut_beam%forb(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%bion))THEN
            ALLOCATE(neut_beam%bion(ke_bm,neut_beam%nbeams))
            neut_beam%bion(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%bneut))THEN
            ALLOCATE(neut_beam%bneut(ke_bm,neut_beam%nbeams))
            neut_beam%bneut(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%hibr))THEN
            ALLOCATE(neut_beam%hibr(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%hibr(:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%hdep))THEN
            ALLOCATE(neut_beam%hdep(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%hdep(:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%zeta))THEN
            ALLOCATE(neut_beam%zeta(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))
            neut_beam%zeta(:,:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%ebeam))THEN
            ALLOCATE(neut_beam%ebeam(ke_bm,neut_beam%nbeams))
            neut_beam%ebeam(:,:) = zeroc
         ENDIF

         IF(.NOT. ASSOCIATED(neut_beam%prompt_pwr_in_plasma))THEN
            ALLOCATE(neut_beam%prompt_pwr_in_plasma(ke_bm,neut_beam%nbeams))
            neut_beam%prompt_pwr_in_plasma(:,:) = zeroc
         ENDIF


         IF(.NOT. ASSOCIATED(neut_beam%fbcur))THEN
            ALLOCATE(neut_beam%fbcur(ke_bm,neut_beam%nbeams))
            neut_beam%fbcur(:,:) = zeroc
         ENDIF
          
         IF(.NOT. ASSOCIATED(neut_beam%rhog_beam))THEN
            ALLOCATE(neut_beam%rhog_beam(neut_beam%nj_beam))
            neut_beam%rhog_beam(:) = zeroc
         ENDIF
          
      ENDIF


      RETURN


      END SUBROUTINE neut_beam_allocate
 END MODULE neutral_beams
