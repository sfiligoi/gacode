     MODULE P_Nfreya_12_interface


       USE nrtype,                ONLY : DP,I4B

       USE common_constants,      ONLY : zeroc,izero 

       IMPLICIT NONE
       PUBLIC 
       LOGICAL use_P_Nfreya,P_Nfreya_read,sent_ufile

       CHARACTER(LEN=256) P_Nfreya_run_directives,P_Nfreya_run_directives_base

       DATA P_Nfreya_run_directives       /'P_nfreya_run_dir_namelist'/
       DATA P_Nfreya_run_directives_base  /'P_nfreya_run_dir_namelist'/
       DATA  use_P_Nfreya /.FALSE./

     CONTAINS


       SUBROUTINE map_nfreya_data(nf_zones,nubeam_namelist,nubeam_ufile,    &
                                  adas_xsct,nameb_nb)
!---------------------------------------------------------------------
! -- given Nfreya related  data in inone copy it  into the slots
! -- required to write the namelist
!---------------------------------------------------------------------

         USE param,                 ONLY : ke
 
         USE fusion,                ONLY : fdbeam

         USE nub,                   ONLY : anglev,angleh,sfrac1,mf,               &
                                           npart,nsourc,nbeams,ranseed,           &
                                           ebkev,fbcur,bptor,bvofset,bheigh,      &
                                           bhofset,alen,bleni,blenp,bwidth,       &
                                           bhfoc,bvfoc,bhdiv,bvdiv,bcur,nashape,  &
                                           nbshape,rpivot,zpivot,aheigh,awidth,   &
                                           naptr,nbeams,fe_tk,fionx

         USE nub2,                  ONLY : iborb

         USE nub3,                  ONLY : iexcit,kdene,ncont,mstate,kdeni,       &
                                           ngl,iz,izstrp,ngh,ilorent,ksve,ksvi,   &
                                           ksvz,krad,kdenz

         USE ions,                  ONLY : nameb

         USE transp,                ONLY : use_ufile,BEAM_DATA_UFILE,             &
                                           BEAM_DATA_NAMELIST,nubeam_back_delt,   &
                                           nubeam_namelist_cpy

         USE neutral_beams,         ONLY :                                          &
                                           anglev_nb => anglev, angleh_nb => angleh,    &
                                           use_ufile_nb => use_ufile,                   &
                                           sfrac1_nb => sfrac1,npart_nb => npart,       &
                                           beam_sim_time_start,beam_sim_time_end,       &
                                           calc_variance,write_performance_data,        &
                                           split_injectors, nfreya_plot_file_name,      &
                                           no_injectors,iterate_beam,iborb_nb =>iborb,  &
                                           randomize_seed ,fap_nb => fap,               &
                                           ebkev_nb => ebkev,fbcur_nb => fbcur,         &
                                           bptor_nb => bptor,bvofset_nb => bvofset,     &
                                           bheigh_nb => bheigh,bhofset_nb => bhofset,   &
                                           alen_nb => alen ,bleni_nb => bleni,          &
                                           blenp_nb => blenp,bwidth_nb => bwidth,       &
                                           bhfoc_nb => bhfoc, bvfoc_nb  => bvfoc,       &
                                           bhdiv_nb => bhdiv, bvdiv_nb => bvdiv,        &
                                           bcur_nb => bcur,nashape_nb => nashape,       &
                                           nbshape_nb => nbshape,rpivot_nb => rpivot,   &
                                           zpivot_nb => zpivot,aheigh_nb => aheigh,     &
                                           awidth_nb => awidth, nsourc_nb =>nsourc,     &
                                           naptr_nb => naptr,nbeams_nb => nbeams,       &
                                           iexcit_nb => iexcit,ilorent_nb => ilorent,   &
                                           mstate_nb => mstate,kdene_nb => kdene,       &
                                           ncont_nb => ncont,kdeni_nb => kdeni,         &
                                           ngl_nb => ngl,ngh_nb => ngh,                 &
                                           fe_tk_nb => fe_tk,krad_nb => krad,           &
                                           ksve_nb => ksve, ksvi_nb => ksvi,            &
                                           ksvz_nb => ksvz,kdenz_nb => kdenz,           &
                                           izstrp_nb => izstrp
                                           
                                           


         USE io_gcnmp,                  ONLY : switch_statefile_output,                &
                                               statefile_output_name


         USE solcon,                    ONLY : time,time0



         IMPLICIT NONE

         REAL(DP) nubeam_dt
         INTEGER(I4B) sz_12,sz_nf,sz,sap,sap_12,sap_nf,iap,jb,ks,kee
         INTEGER(I4B) nf_zones,dclnf,dclnf2
         CHARACTER*(*) nubeam_namelist,nubeam_ufile,adas_xsct,nameb_nb
 
!-------------------------------------------------------------------------------------
! -- items that should  be obtained from inone:
!-------------------------------------------------------------------------------------
         write_performance_data    = "'performance_data.txt'"    ! should set in inone
         split_injectors           = .FALSE.    ! should set in inone
         nfreya_plot_file_name     = "'P_NF_bpltfil'" 
         iterate_beam              = .FALSE. ! this is should be left off for P_Nfreya
         adas_xsct                 = "'/usc-data/p2/linux/onetwo/'" 



      
         ! Assume that if ranseed is left at its default value of 7**7
         ! then we do not want to randomize the seed:
           randomize_seed            = .TRUE.
           IF(ranseed == 7**7) randomize_seed = .FALSE.


           dclnf = LEN(nbshape_nb(1))   ! declared length of nbshape in neutral_beams, may be different
           dclnf2 = LEN(nashape_nb(1,1))
           sz = LEN_TRIM(nameb)

         nameb_nb                  ='"'//nameb(1:sz)//'"'

         iborb_nb                  = iborb
         IF(iborb_nb .GT. 1)iborb  = 1 ! iborb =2 not implemented in P_Nfreya
         no_injectors              = nbeams
         nubeam_ufile              = '"'//BEAM_DATA_UFILE(1:LEN_TRIM(BEAM_DATA_UFILE))//'"'
         !nubeam_namelist           = '"'//BEAM_DATA_NAMELIST(1:LEN_TRIM(BEAM_DATA_NAMELIST))//'"'
         nubeam_namelist           = '"'//TRIM(nubeam_namelist_cpy)//'"'



         use_ufile_nb              = USE_UFILE 
         nf_zones                  = mf
         npart_nb                  = npart
         calc_variance             = .TRUE. 

         statefile_output_name(:)  = ' '

         switch_statefile_output   = 0 ! For P_Nfreya we dont allow toggling
                                       ! from .nc to .txt and vice versa
         


         sz_12   = SIZE(anglev)
         sz_nf   = SIZE(anglev_nb)
         sz     = MIN(sz_12,sz_nf)
         sap_12 = SIZE(alen,1)
         sap_nf = SIZE(alen_nb,1)
         sap    = MIN(sap_12,sap_nf)


         iexcit_nb                = iexcit
         ilorent_nb               = ilorent
         mstate_nb                = mstate
         nbeams_nb                = nbeams
         nsourc_nb                = nsourc
         naptr_nb                 = naptr
         fe_tk_nb                 = fe_tk
         ngl_nb                   = ngl
         ngh_nb                   = ngh
         kdene_nb                 = kdene
         kdeni_nb                 = kdeni
         ksvz_nb                  = ksvz
         ksvi_nb                  = ksvi
         ksve_nb                  = ksve
         krad_nb                  = krad
         kdenz_nb                 = kdenz
         ncont_nb                 = ncont

         anglev_nb(1:sz)          = anglev(1:sz)
         angleh_nb(1:sz)          = angleh(1:sz)
         sfrac1_nb(1:sz)          = sfrac1(1:sz)
         ebkev_nb(1:sz)           = ebkev(1:sz)
         fbcur_nb(1:ke,1:sz)      = fbcur(1:ke,1:sz)
         bptor_nb(1:sz)           = bptor(1:sz)
         bvofset_nb(1:sz)         = bvofset(1:sz)
         bheigh_nb(1:sz)          = bheigh(1:sz) 
         bhofset_nb(1:sz)         = bhofset(1:sz)

        
         bleni_nb(1:sz)           = bleni(1:sz)
         blenp_nb(1:sz)           = blenp(1:sz)
         bwidth_nb(1:sz)          = bwidth(1:sz)
         bvdiv_nb(1:sz)           = bvdiv(1:sz)
         bhdiv_nb(1:sz)           = bhdiv(1:sz)
         bvfoc_nb(1:sz)           = bvfoc(1:sz)
         bhfoc_nb(1:sz)           = bhfoc(1:sz)
         bcur_nb(1:sz)            = bcur(1:sz)

         zpivot_nb(1:sz)          = zpivot(1:sz)
         rpivot_nb(1:sz)          = rpivot(1:sz)


         awidth_nb(1:sap,1:sz)    = awidth(1:sap,1:sz)
         aheigh_nb(1:sap,1:sz)    = aheigh(1:sap,1:sz) 
         alen_nb(1:sap,1:sz)      = alen(1:sap,1:sz)
        nbshape_nb(:)(1:dclnf)    = ' '
        nashape_nb(:,:)(1:dclnf2) = ' '
         DO jb =1,sz
           nbshape_nb(jb)(1:LEN_TRIM(nbshape(jb))) =               &
                             nbshape(jb)(1:LEN_TRIM(nbshape(jb)))
           DO iap = 1,sap
              nashape_nb(iap,jb)(1:LEN_TRIM(nashape(iap,jb))) =      &
                          nashape(iap,jb)(1:LEN_TRIM(nashape(iap,jb)))
           ENDDO
         ENDDO


!---------------------------------------------------------------
! we need to fill in slots not occupied. Otherwise the namelist
! output will not be readable
!----------------------------------------------------------------
         IF(sap_nf .GT. sap_12 )THEN
            DO iap = sap+1,sap_nf
              awidth_nb(iap ,1:sz) =  awidth_nb(1,1:sz)
              aheigh_nb(iap ,1:sz) =  aheigh_nb(1,1:sz)
              nashape_nb(iap,1:sz) =  nashape_nb(1,1:sz)
              alen_nb(iap,1:sz)    =  alen_nb(1,1:sz)
            ENDDO
         ENDIF
         IF(sz_nf .GT. sz_12)THEN
            DO jb = sz+1,sz_nf
              DO iap = 1,sap_nf
                 awidth_nb(iap ,jb) =  awidth_nb(iap,1)
                 aheigh_nb(iap ,jb) =  aheigh_nb(iap,1)
                 nashape_nb(iap,jb) =  nashape_nb(iap,1)
                 alen_nb(iap,jb)    =  alen_nb(iap,1)
              ENDDO
              zpivot_nb(jb)         = zpivot(1)
              rpivot_nb(jb)         = rpivot(1)
              nbshape_nb(jb)        = nbshape_nb(1)
              bcur_nb(jb)           = bcur(1)
              bhfoc_nb(jb)          = bhfoc_nb(1)
              bvfoc_nb(jb)          = bvfoc_nb(1)
              bhdiv_nb(jb)          = bhdiv_nb(1)
              bvdiv_nb(jb)          = bvdiv_nb(1)
              bwidth_nb(jb)         = bwidth_nb(1)
              blenp_nb(jb)          = blenp_nb(1)
              bleni_nb(jb)          = bleni_nb(1)
              bhofset_nb(jb)        = bhofset_nb(1)
              bheigh_nb(jb)         = bheigh_nb(1)
              bvofset_nb(jb)        = bvofset_nb(1)
              bptor_nb(jb)          = bptor_nb(1)
              fbcur_nb(1:ke,jb)     = fbcur(1:ke,1)
              ebkev_nb(jb)          = ebkev_nb(1)
              sfrac1_nb(jb)         = sfrac1_nb(1)
              angleh_nb(jb)         = angleh_nb(1)
              anglev_nb(jb)         = anglev_nb(1)
            ENDDO
         ENDIF
         DO jb =1,sz_nf
            DO iap = 1,sap_nf
               nashape_nb(iap,jb)  = '"'//nashape_nb(iap,jb)(1:LEN_TRIM(nashape_nb(iap,jb)))//'"'
            ENDDO
            !The quotes are now taken care of in sub add_quotes_to_namelist
            !nbshape_nb(jb)         = '"'//nbshape_nb(jb)(1:LEN_TRIM(nbshape_nb(jb)))//'"'
         ENDDO


         RETURN
       END SUBROUTINE map_nfreya_data



       SUBROUTINE P_Nf_events
!---------------------------------------------------------------------
! -- Set up a time array that holds the beam switching on and off
! -- event times. We do not need to know if the switching time is 
! -- a beam on or a beam off time. All we need is that a change occured.
! -- INPUT
!        beam_data%tbona
!        beam_data%tboffa
!
!
! -- OUTPUT
!        beam_data%beam_sw sorted in ascending time
!      
!---------------------------------------------------------------------


         USE transp,                        ONLY : beam_data

         IMPLICIT NONE
         INTEGER(I4b) j,jb,jb1,jb2,km,kmm,i,k
         INTEGER(I4B) INDEX(2*beam_data%nbeam)
         REAL(DP) beam_sw(2*beam_data%nbeam)



         jb1 = beam_data%nbeam
     
         jb2 = 2*jb1

         beam_sw(1:jb1)     = beam_data%tbonac(1:jb1)
         beam_sw(jb1+1:jb2) = beam_data%tboffac(1:jb1)

         DO j=1,jb2
            INDEX(j) = j
         ENDDO

         !CALL sort_unique(beam_sw,jb2)
         CALL bsort(beam_sw,jb2,index,.TRUE.)
         IF(jb1 .GT.1)THEN
            km=jb2
            DO WHILE( 1 .GT. 0 ) !lf95 doesnt like DO WHILE(1) 
               kmm=0
               DO j=2,km
                  IF(beam_sw(j) .EQ. beam_sw(j-1))THEN
                     DO i=j+1,jb2
                        beam_sw(i-1) = beam_sw(i)
                        INDEX(i-1)   = INDEX(i)
                     ENDDO 
                     km = km-1
                     kmm=1
                     EXIT
                  ENDIF
               ENDDO
               IF(kmm .EQ. 0) EXIT
            ENDDO
         ENDIF


         IF(ASSOCIATED(beam_data%beam_sw))DEALLOCATE(beam_data%beam_sw)
         ALLOCATE(beam_data%beam_sw(jb2))
         beam_data%beam_sw(1:jb2) = beam_sw(1:jb2)

         RETURN

       END SUBROUTINE P_Nf_events




       SUBROUTINE process_beamlets(beamlet_ct,beamlet_nosw,time,run_P_Nfreya, &
                                   recall_P_Nfreya)
!-----------------------------------------------------------------------------------
! -- Returns 
! -- beam_data%beamlet_active(jb) =  1  if beamlet jb is currently on
! -- beam_data%beamlet_active(jb) =  0  if beamlet jb hs not yet been turned on 
! -- beam_data%beamlet_active(jb) = -1  if beamlet was on but is now off
! -- beam_data%beamlet_cp(:) holds state of beam_data%beamlet_active on exit from this routine
! -- jb = 1,..beam_data%nbeam
! -- beamlet_ct     total  number of active and possibly now inactive 
! --                because it was on but is now off turned off
! --                (=0 if no active beamlets are found up to ths time)
! -- recall_P_Nfreya = T or F 
! -- run_P_Nfreya 1 if p_Nfreya should be run
!                -1 otherwise
! -- beamlet_nosw 
!---------------------------------------------------------------HSJ---2/14/2012-----
       USE transp,                          ONLY : beam_data

       USE nub2

       USE events,                          ONLY : isave_event,       &
                                                   index_p_nf_start,  &
                                                   index_p_nf_stop

       USE solcon,                          ONLY : time_tol

       USE tordlrot,                        ONLY : spbrsav,pprbsav
       IMPLICIT NONE
       LOGICAL recall_P_Nfreya
       INTEGER(i4b) beamlet_ct,run_P_Nfreya,beamlet_nosw,jb,icall
       REAL(DP) time,dton,dtoff

      DATA icall /izero/


             run_P_Nfreya = -1
             recall_P_Nfreya = .FALSE.
             beamlet_ct = izero        ! remains at zero if no beamlets are or were on

             beam_data%beamlet_active(:) = 0
             DO   jb=1,beam_data%nbeam
                ! beamlet is curently on:
                dton  = ABS(time - beam_data%tbonac(jb) )
                dtoff = ABS(time - beam_data%tboffac(jb)) 
                IF(time .GE. beam_data%tbonac(jb)-time_tol .AND. time .LT. & 
                                    beam_data%tboffac(jb)-time_tol)THEN
                   beam_data%beamlet_active(jb) = 1
                   run_P_Nfreya = 1
                   beamlet_ct  = beamlet_ct +1
                ENDIF

                ! beamlet was on now off ==> zero source
                IF(time .GE. beam_data%tboffac(jb)-time_tol                &
                             .AND. beam_data%beamlet_active(jb) == 1)THEN
                   beam_data%beamlet_active(jb) = -1
                   sbsav(1:kj,1:ke,jb)   = zeroc
                   qbsav(1:kj,1:ke,jb)   = zeroc
                   spbsav(1:kj,1:ke,jb)  = zeroc
                   spbrsav(1:kj,1:ke,jb) = zeroc
                   pb0(1:kj,1:ke,jb)     = zeroc
                   beamlet_ct            = beamlet_ct +1
                ENDIF

                ! beamlet not yet on:
                IF(time .LT. beam_data%tbonac(jb)-time_tol)THEN
                   ! zero the sources
                   beam_data%beamlet_active(jb)    = izero
                   sbsav(1:kj,1:ke,jb)   = zeroc
                   qbsav(1:kj,1:ke,jb)   = zeroc
                   spbsav(1:kj,1:ke,jb)  = zeroc
                   spbrsav(1:kj,1:ke,jb) = zeroc
                   pb0(1:kj,1:ke,jb)     = zeroc
                   ! zero profiles derived from sources:
                   enbsav(1:kj,1:ke,jb)  = zeroc
                   wbsav(1:kj,1:ke,jb)   = zeroc
                   ppbsav(1:kj,1:ke,jb)  = zeroc
                   pprbsav(1:kj,1:ke,jb) = zeroc
                ENDIF
             ENDDO

!             beamlet_nosw = izero

!             IF(index_p_nf_start .le. isave_event .and.         &
!                              isave_event .le. index_p_nf_stop)THEN
!                   run_P_Nfreya = 1
!                   IF(beamlet_ct .GT. izero) recall_P_Nfreya = .TRUE.
                   ! get beamlet no causing this trigger:
                   ! make sure this beamlet is switched on or off as
                   ! the case may be. account for other beamlets that
                   ! might have the same switching time as well:
!                   beamlet_nosw = (isave_event - index_p_nf_start)/2+1
 
!              ENDIF


       DO jb =1,beam_data%nbeam
          IF(beam_data%beamlet_cp(jb) .NE.                  &
                     beam_data%beamlet_active(jb) .AND.     &
                     beamlet_ct .GT. izero)THEN     
               recall_P_Nfreya = .TRUE.
               run_P_Nfreya    = 1
          ENDIF
       ENDDO

       beam_data%beamlet_cp(:) = beam_data%beamlet_active(:)

!  print *,'beamlet_nosw,time =', beamlet_nosw,time  ! 888899
!  print *,'isave_event,index_p_nf_s,e =',isave_event,index_p_nf_start,index_p_nf_stop
!  print *,'recall_P_Nfreya,beamlet_ct =',recall_P_Nfreya,beamlet_ct


      icall = icall +1
!  write(977,5)icall
!  write(977,1)time,run_P_Nfreya,recall_P_Nfreya
1 FORMAT(2x,'time,run_P_Nfreya,recall_P_Nfreya =',1pe12.6,x,i2,x,l8)
!  write(977,2)beamlet_ct,beamlet_nosw
!  print *,'beamlet_nosw ******************************=',beamlet_nosw
2  FORMAT(2x,'beamlet_ct,beamlet_nosw =',i5,2x,i5)
!  write(977,3)(beam_data%beamlet_active(jb),jb=1,beam_data%nbeam)
3 FORMAT('bm_actv =',10(i2,2x))
!  write(977,4)isave_event,index_p_nf_start,index_p_nf_stop
!  write(977,6)sbsav(1,1:3,1)
4 FORMAT('is ,start,stop =',3(2x,i5))
5 FORMAT('call no ',i5)
6 FORMAT('sbsav bm1 =',3(2x,1pe12.4))
          RETURN


      END SUBROUTINE process_beamlets



      SUBROUTINE bsort (a,n,index,ascend)
!-------------------------------------------------------------------
! -- bubble sort array a in ascending ascend = true
! -- or descneding (ascend = false) order 
! -- change array index in same way as array a is changed
!------------------------------------------------------------HSJ-------
      USE nrtype,                          ONLY : DP,I4B



      IMPLICIT NONE


      INTEGER(I4B) i,j,jm,itemp,n
      LOGICAL ascend,switch
      REAL(DP), DIMENSION(n) :: a
      INTEGER(I4B), DIMENSION(n) :: index

      REAL(DP) temp

      IF(n .GT. 1)THEN ! return with a(1),index(1) =1 otherwise
         IF(.NOT. ascend)THEN
            jm=n-1
            DO i=1,n-1
               switch = .FALSE.
               ! each time the j do loop is completed
               ! the curently smalles  element will be
               ! in the correct position:
               DO  j=1,jm
                  IF(a(j) .GT. a(j+1))CYCLE
                  temp=a(j)
                  a(j)=a(j+1)
                  a(j+1)=temp
                  itemp=INDEX(j)
                  INDEX(j)=INDEX(j+1)
                  INDEX(j+1)=itemp
                  switch = .TRUE.
               ENDDO
               IF(.NOT. switch ) EXIT
               jm = jm-1
            ENDDO
         ELSE
            jm=n-1
            DO i=1,n-1
               switch = .FALSE.
               ! each time the j do loop is completed
               ! the curentl largest element will be
               ! in the correct position:
               DO  j=1,jm
                  IF(a(j) .LT. a(j+1))CYCLE
                  temp=a(j)
                  a(j)=a(j+1)
                  a(j+1)=temp
                  itemp=INDEX(j)
                  INDEX(j)= INDEX(j+1)
                  INDEX(j+1)= itemp
                  switch = .TRUE.
               ENDDO
               IF(.NOT. switch ) EXIT
               jm = jm-1
            ENDDO
         ENDIF
       ENDIF

        RETURN

      END SUBROUTINE bsort


     END MODULE P_Nfreya_12_interface
