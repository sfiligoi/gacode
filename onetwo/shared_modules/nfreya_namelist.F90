   MODULE Nfreya_namelist

    USE nrtype,              ONLY : DP,I4B,I2B

    USE nf_param,            ONLY : kb

    USE error_handler,       ONLY : iomaxerr,lerrno,terminate

    USE Plasma_properties ,  ONLY : mhd_dat,pellet

    USE solcon_gcnmp,        ONLY : time0,time,time_max,dt,dtmin,                     &
                                    dtmax,use_glf23, dtstart,max_steps,tGCNMf,tGCNMs, &
                                    fdigits,gradtol,steptol,bandwidth,jacobian_type,  &
                                    fvectol,gradtol,ssqrmin,maxfev,steady_state,      &
                                    conv_skip,wrt_nwt,outputvb,theta,                 &
                                    nsteps_plot,plot_time_incr,save_incr_plot_file,   &
                                    nsteps_restart,create_restart_file,prt_turbm,     &
                                    restart_time_incr,save_incr_restart_file

    USE common_constants,      ONLY : izero,zeroc

    USE grid_class,            ONLY : nj,nj_out,allow_regrid,reset_cap_parms,       &
                                      set_cap1,nj_max,nj_start,use_compact_schemes, &
                                      frac_core_points

    USE io_gcnmp,              ONLY : namelist_filename,nlog,ncrt,                  &
                                      rpc,iterdb_12_input_file,maxchr,              &
                                      switch_statefile_output,                      &
                                      statefile_output_name,maxchr

    USE ions_gcnmp,            ONLY : nprim,nimp,nion,name_size,ni_sc

    USE curden_terms,          ONLY : jboot,rbsaxis,irfc,ibcur,currf,curbeam

    USE bc_values_gcnmp,       ONLY : u_vloop_bc,mult_den,mult_flux,gammat_bc,    &
                                      axis_bc_type

    USE neutral_beams,         ONLY :                                          &
                                       timbplt, beam_on, beam_time,            &
                                       statefile_nameb => nameb,iborb,         &
                                       namelist_nameb,namelist_nbion,          &
                                       anglev, angleh,nashape, aheigh,relnub,  &
                                       awidth, bcur, bptor, blenp, nbshape,    &
                                       bleni,bheigh, bwidth, bhfoc, bvfoc,     &
                                       bhdiv, bvdiv, ebkev, fbcur, nbeams,     &
                                       naptr, alen, bvofset, bhofset, nsourc,  &
                                       sfrac1, npart, npskip, rpivot,          &
                                       zpivot, randomize_seed, fionx,fd_beam,  &
                                       ilorent,iexcit,ngl, kdeni, kdenz, ksvi, &
                                       ksvz, ksve, krad, ngh,ncont,kdene,ne_tk,&
                                       fe_tk,ds_tk,fdbeam,izstrp,              &
                                       nameb_ml,nbeam_ml,nameb_index,          &
                                       iterate_beam,ilorent,mstate,            &
                                       npart_all_beamlines,no_injectors,       &
                                       beam_sim_time_start,beam_sim_time_end,  &
                                       nfreya_plot_file_name,time_dep_beam,    &
                                       split_injectors,write_performance_data, &
                                       nfreya_vb,calc_variance,use_ufile


    USE nub,                   ONLY:   hdepsmth,fidiff_on,bfr_neutrlz

#ifndef ONETWO
    USE P_nfreya_interface,    ONLY:   beam_data_namelist,beam_data_ufile

    USE source_terms_gcnmp,    ONLY:   beam_th_mult, stfuse_mult, sbfuse_mult,qrad_mult

    USE xsct,                  ONLY:   adas_xsct_path

    USE zonal_data,            ONLY:   mf,mfm1

    USE MPI_data,              ONLY:   myid,master,parallel_model
#else 
    USE string_util,           ONLY:  to_upper_case1,to_lower_case1
#endif

    IMPLICIT NONE 

    INTEGER(I2B)  ionml,ionml_temp
    INTEGER(I2B), EXTERNAL :: get_next_io_unit 
    INTEGER(I4B)  j,ienct, nf_zones
    INTEGER, PARAMETER :: nion_max = 10
 
    REAL(DP) irfc_mult,ibcur_mult,nictr,niavg,nistd,nirho
 

    REAL(DP) pellet_rmin,pellet_dep_width,pellet_np,pellet_freq

    REAL(DP) P_Nfreya_dt

!    CHARACTER(LEN = LEN(namelist_filename)+2) name_temp ! gfortran doesnt like this
    CHARACTER(LEN = maxchr +2) name_temp 
    CHARACTER(LEN =132) jbootstrap

! local definitions to satisfy Onetwo and P_Nfreya 
    CHARACTER(LEN = 256) adas_xsct 
    CHARACTER(LEN = 256) nubeam_namelist
    CHARACTER(LEN = 256) nubeam_ufile

   
    CHARACTER(len = name_size),PUBLIC ::  nameb  !nameb from namelist (not statefile)
    LOGICAL pellet_inject,nameb_set



    NAMELIST /run_data/                                                       &
       timbplt, beam_on, beam_time, nameb, relnub,  anglev, angleh,           &
       nashape, aheigh, awidth, bcur, bptor, blenp, nbshape, bleni,           &
       bheigh, bwidth, bhfoc, bvfoc, bhdiv, bvdiv, ebkev, fbcur,              &
       nbeams, naptr, alen, bvofset, bhofset, nsourc, sfrac1,nf_zones,        &
       npart, npskip, rpivot, zpivot, randomize_seed, fionx,  fdbeam,iexcit,  &
       mstate,fe_tk,izstrp,adas_xsct,iborb,iterate_beam,ilorent,ncont,        &
       ngl, kdeni, kdenz, ksvi,ksvz, ksve, krad, ngh,kdene,                   &
       npart_all_beamlines,no_injectors,beam_sim_time_start,                  &
       beam_sim_time_end,nfreya_plot_file_name,split_injectors,               &
       write_performance_data,nfreya_vb,calc_variance,                        &
       switch_statefile_output,statefile_output_name,nubeam_namelist,         &
       nubeam_ufile,use_ufile,time0,time_max,hdepsmth,fidiff_on,              &
       bfr_neutrlz,ne_tk


     CONTAINS 

        SUBROUTINE read_run_directives_namelist
! ------------------------------------------------------------------------------
! SET SOME RUN DIRECTIVE DEFAULTS AND read namelist "run_data" 
! NOTE: This routine is called only by the master process
! ---------------------------------------------------------------------HSJ-05/19/05

    INTEGER(I4B)  glf_flx_ctr




#ifndef ONETWO
         ! utils.f90:
    INTERFACE
      SUBROUTINE to_upper_case(string)
         USE nrtype,            ONLY : I4B
         IMPLICIT NONE
         INTEGER(I4B) l
         CHARACTER*(*), INTENT (INOUT) :: string
      END SUBROUTINE to_upper_case


      SUBROUTINE to_lower_case(string) 
           USE nrtype,            ONLY : I4B
           IMPLICIT NONE
           INTEGER(I4B) l
           CHARACTER*(*), INTENT (INOUT) :: string
      END SUBROUTINE to_lower_case
    END INTERFACE
#else
      !string_util.F90  interface declaration wont compile
!    INTERFACE
!      SUBROUTINE  to_lower_case1(string)
!            CHARACTER*(*), INTENT (INOUT) :: string
!            INTEGER l
!      END SUBROUTINE  to_lower_case1
!
!      SUBROUTINE  to_upper_case1(string)
!            CHARACTER*(*), INTENT (INOUT) :: string
!!            INTEGER l
!      END SUBROUTINE  to_upper_case1
!    END INTERFACE 
#endif

 
    !set local namelist defaults
        nictr = -1._DP ;  niavg = -1._DP  ;   nistd = -1._DP ; nirho =-1.
        nameb = namelist_nameb(1)





 
    adas_xsct       = ''       ! local variable default to nothing
    nubeam_namelist = ''
    nubeam_ufile    = ''
    nf_zones        = 0

    ionml      = get_next_io_unit()
    ionml_temp = get_next_io_unit()



    name_temp = ADJUSTL(namelist_filename(1:LEN_TRIM(namelist_filename)))

#ifndef ONETWO
        CALL to_upper_case(name_temp)
#else


        CALL to_upper_case1(name_temp)
#endif
    IF(name_temp(1:LEN_TRIM(name_temp)) .EQ.                &
         namelist_filename(1:LEN_TRIM(name_temp)))          &
         name_temp = name_temp(1:LEN_TRIM(name_temp))//'_1' 
    OPEN (unit = ionml, file = namelist_filename, status = 'OLD',err = 1)
    OPEN (unit = ionml_temp,file =name_temp,status ='UNKNOWN',err =1)
    CALL strip_comments (ionml, ionml_temp)
    CLOSE(unit = ionml,status = 'KEEP')
    REWIND(UNIT=ionml_temp)
    READ(ionml_temp,nml=run_data,END = 2) !  better to let it
                                          !  fail than redirect with ERR = 
    CLOSE(unit = ionml_temp,status = 'DELETE')

#ifndef ONETWO             ! avoids having to bring some modules into Onetwo
    ! use default values if not set in namelist:
    IF(LEN_TRIM(adas_xsct) .GT. 0)                       & 
         adas_xsct_path(1:LEN_TRIM(adas_xsct))           &
                        = adas_xsct(1:LEN_TRIM(adas_xsct))   

    IF(LEN_TRIM(nubeam_namelist) .GT. 0)THEN 
         beam_data_namelist(1:LEN(beam_data_namelist))=''
         beam_data_namelist(1:LEN_TRIM(nubeam_namelist)) &
                        = nubeam_namelist(1:LEN_TRIM(nubeam_namelist))
    ENDIF
    IF(LEN_TRIM(nubeam_ufile) .GT. 0)THEN
       beam_data_ufile(1:LEN(beam_data_ufile))=''
         beam_data_ufile(1:LEN_TRIM(nubeam_ufile))    &
                        = nubeam_ufile(1:LEN_TRIM(nubeam_ufile))
    ENDIF


    IF(nf_zones .GT. 0)mf = nf_zones
    mfm1      = mf - 1

    IF(ne_tk .LE. izero)ne_tk = 20 ! same as default value    
    IF(fe_tk .le. 1.1) fe_tk = 1.4
#endif

    IF(npart .GT. izero)THEN
       npart_all_beamlines = npart
    ELSEIF(npart_all_beamlines .GT. izero)THEN
       npart = npart_all_beamlines
    ELSE
       lerrno = 229 + iomaxerr
       CALL terminate(lerrno,nlog)
    ENDIF

    IF(nbeams .GT. izero)THEN
       no_injectors  = nbeams
    ELSEIF(no_injectors .GT. izero)THEN
       nbeams  = no_injectors
    ELSE
       lerrno = 230 + iomaxerr
       CALL terminate(lerrno,nlog)
    ENDIF

    IF(kb .LT. nbeams)THEN
       lerrno = iomaxerr + 232_I4B
       call terminate(lerrno,nlog)
    ENDIF

    dt = dtstart

    ni_sc(1)  = nictr ; ni_sc(2) = niavg ; ni_sc(3)= nistd ; ni_sc(4) = nirho



!    IF(time0 .LT. zeroc)THEN
!        time0 = time          !set initial time to iterdb file time
!        WRITE(nlog,FMT='(" start time =" ,f12.6," taken from iterdb file")')time0
!    ELSE
!       WRITE(nlog,FMT='(" start time =" ,f12.6," taken from namelist file")')time0
!    ENDIF
!    WRITE(nlog,FMT ='(" max number of time steps allowed = ",i12,/)')max_steps
!    tGCNMs = time0

!    IF(time_max .LT. zeroc .AND. steady_state == 1 )THEN
!       time_max = tGCNMf !set final time to iterdb file time
!       WRITE(nlog,FMT='(" end  time = ",f12.6," taken from iterdb file")')time_max
!    ELSEIF(time_max .LT. zeroc .AND. steady_state == 0 )THEN
!       WRITE(nlog,FMT='(" run into steady state")')
!    ELSE
!      IF(steady_state == 1)THEN
!          WRITE(nlog,FMT='(" end time = ",f12.6," taken from namelist file")')time_max
!      ELSE
!          WRITE(nlog,FMT='(" run into steady state")')
!      ENDIF
!    ENDIF



!    IF(time_max .LT. time0) lerrno = iomaxerr + 1_I4B
!    IF(lerrno .GT. 0) CALL terminate(lerrno,nlog)
!    tGCNMf = time_max


    IF(LEN_TRIM(statefile_output_name) .GT. maxchr)lerrno = iomaxerr + 245 
    IF(lerrno .GT. 0) CALL terminate(lerrno,nlog)
 
    IF(set_cap1 .GT. 0 ) CALL reset_cap_parms(set_cap1)

    IF(ABS(irfc_mult -FLOAT( irfc)) .GT. 1.e-5 .AND. irfc .NE. 0)THEN
       IF(ABS(irfc_mult) .GT. 1.e-8)THEN
          currf(:) = currf(:)*irfc_mult
       ELSE
          currf(:) = 0.0_DP
          irfc = 0
       ENDIF

    ENDIF




!---------------------------------------------------------------------------
! -- Sort out beam species issues
!---------------------------------------------------------------------------
     IF(nameb =='none')THEN
        WRITE(ncrt,FMT ='("ERROR: sub read_namelist, nameb not correct")')
        lerrno = 214 + iomaxerr
        CALL terminate(lerrno,nlog)
     ELSE
       ! use lower case in rest of code:
#ifndef ONETWO
        CALL to_lower_case(nameb)
#else
        CALL to_lower_case1(nameb)
#endif
        nameb_set = .FALSE.
        namelist_nbion = 1
        namelist_nameb(1) = nameb
        DO j=1,nbeam_ml 
           IF(nameb == nameb_ml(j))THEN ! ALLOWED are  'h', 'd', 't', 'dt' 
               nameb_set = .TRUE.
               nameb_index(j) = j ! use this info afte reading primary 
                                  ! ion names in statefile
           ENDIF
        ENDDO
        IF(.NOT. nameb_set)THEN
           WRITE(ncrt,FMT ='("ERROR: sub read_namelist, nameb not correct")')
           lerrno = 214 + iomaxerr
           CALL terminate(lerrno,nlog)
        ENDIF
     ENDIF


     IF((beam_sim_time_start .LT. -Huge(1._DP)/2._DP .OR. beam_sim_time_end .LT. -Huge(1._DP)/2._DP) &
           .OR. (beam_sim_time_end .LE. beam_sim_time_start)) THEN
           WRITE(ncrt,FMT ='("ERROR: sub read_namelist, beam_sim_time_start and/or end not set")')
           lerrno = 231 + iomaxerr
           CALL terminate(lerrno,nlog)

     ENDIF

     IF(hdepsmth .LT. 0.0)fidiff_on = .FALSE. ! master overide of fas ion smoothing with diffsion  


#ifndef ONETWO
        CALL to_upper_case(write_performance_data)
#else
        CALL to_upper_case1(write_performance_data)
#endif


  RETURN




1   lerrno = 4_I4B
    CALL terminate(lerrno,nlog)
  RETURN
2   lerrno = iomaxerr + 2_I4B
    CALL terminate(lerrno,nlog)
  RETURN
3   lerrno = iomaxerr + 3_I4B
    CALL terminate(lerrno,nlog)


  RETURN


  END SUBROUTINE read_run_directives_namelist


#ifdef ONETWO
  SUBROUTINE  write_P_Nfreya_run_directives(time)
!--------------------------------------------------------------
! -- this subroutine writes the namelist run directives file
! -- that will be read by the stand  alone P_Nfreya module.
! -- A problem here is that we need to write character data
! -- using the namelist mechanism. This is problematic because
! -- we have to put leading and trailing quotes around the 
! -- character variables before writting the namelist.
! -- This is an issue here beacause  all backward compatible code
! -- since day 1 of Onetwo assumes character variable size of 8. But at least one
! -- quantity (nbshape) is too large to allow for quotes in addition
! -- to the source type(eg "rect-lps" is 10 characters). 
! -- The path of least resistance appears to be
! -- first writting the namelist and then rewritting it again with the
! -- quotes included. (see sub add_quotes_to_namelist) . In this way
! -- we retain  the 8 character limit on nbshape,etc..
!
! -- IMPORTANT NOTE:
! -- Ifort does not dump the namelist character variables nashape,nbshape
! -- correctly. Hence this wont work with ifort ??????!!!!!
!---------------------------------------------HSJ-4/15/11------

     USE nrtype,                           ONLY : DP,I4B,I2B

     USE P_Nfreya_12_interface,            ONLY : P_Nfreya_run_directives, &
                                                  map_nfreya_data

     USE io,                               ONLY : versid


     IMPLICIT NONE
     INTEGER(I2B) io_unit
     REAL(DP) time
     CHARACTER(LEN=96) string
 
     CALL map_nfreya_data(nf_zones,nubeam_namelist,nubeam_ufile,         &
                          adas_xsct,nameb)        ! P_Nfreya_12_interface.f90

     CALL P_Nfreya_output_name(time)


    IF(npart .GT. izero)THEN
       npart_all_beamlines = npart
    ELSEIF(npart_all_beamlines .GT. izero)THEN
       npart = npart_all_beamlines
    ELSE
       lerrno = 229 + iomaxerr
       CALL terminate(lerrno,nlog)
    ENDIF

     io_unit  = get_next_io_unit()



     OPEN (unit = io_unit,file = P_Nfreya_run_directives,status ='UNKNOWN',err =1)
     string = 'Run directives file created by '//versid
     WRITE(io_unit,FMT ='(a)')string
     string = '-----------------------------------------------------'
     WRITE(io_unit,FMT ='(a)')string

     WRITE(io_unit,NML = run_data)
   
 
     CALL add_quotes_to_namelist(io_unit)


     RETURN
1    lerrno = 315 + iomaxerr
     CALL terminate(lerrno,nlog)
  END SUBROUTINE  write_P_Nfreya_run_directives




  SUBROUTINE  add_quotes_to_namelist(io_in)
!-----------------------------------------------------------------------
! -- see explanation in write_P_Nfreya_run_directives
! -- for existence of this routine
!---------------------------------------------------------------HSJ-----

     USE nrtype,                           ONLY : DP,I4B,I2B

     USE error_handler,                    ONLY : iomaxerr,lerrno,terminate

     IMPLICIT NONE 
     INTEGER(I2B) io_in,io_out,get_next_io_unit
     INTEGER(I4B) ll,lr,j,lj
     CHARACTER(LEN=256) line,prt1,prt2,filename,filenamec
     LOGICAL open

     INQUIRE(unit = io_in, OPENED = open,NAME  = filename)
    
     IF(open)THEN
        REWIND(io_in)
        io_out = get_next_io_unit ()
        filenamec = TRIM(filename)//'c'
        OPEN (unit = io_out,file = filenamec ,status ='UNKNOWN',err =2)

     ELSE
1       lerrno = 356  + iomaxerr
        CALL terminate(lerrno,nlog)
     ENDIF 
     DO WHILE(1)
        line(1:LEN(line)) =' '
        READ(io_in,FMT='(a)',END = 10)line
        IF(INDEX(line,'NBSHAPE') .NE. 0 .OR. INDEX(line,'nbshape') .NE. 0)THEN
         ! first part is nbshape = *****
           ll = INDEX(line,"=")
           IF(ll == 0)THEN
              lerrno = 356 + iomaxerr
              CALL terminate(lerrno,nlog)
           ENDIF
 
           lj = -1
           prt1(1:ll) = line(1:ll)
           !search for first non blank character after = sign:
           DO j = ll+1,LEN_TRIM(line)
              IF(line(j:j) == ' ')CYCLE
              lj = j ! first non blank position after "="
              EXIT
           ENDDO
           IF(lj == -1)THEN
              lerrno = 356 + iomaxerr
              CALL terminate(lerrno,nlog)
           ENDIF
           prt2(:) = ' '
!           prt2(1:LEN_TRIM(line)-lj-1) = line(lj:LEN_TRIM(line)-1)
           prt2(1:LEN_TRIM(line)-lj) = line(lj:LEN_TRIM(line))


           line(1:LEN(line)) =' '
           line = prt1(1:ll)//' '//'"'//TRIM(prt2)//'"'//','
           WRITE(io_out,FMT='(a)')TRIM(line)

         ! succeeding lines are of the form    rect-lps  (for example)
         ! process all array elements which appaer in the namelist
         ! until we encounter the next variable
           lj = -1
           DO WHILE(1)     ! read lines until next variable in namelist is found
              prt1(:) = ' '
              READ(io_in,FMT = '(a)')prt1 
              IF(INDEX(TRIM(prt1),'=') .Eq. 0)THEN
                 DO j=1,LEN_TRIM(prt1)
                    IF(prt1(j:j) .NE. ' ')THEN
                       lj = j
                       EXIT
                    ENDIF
                 ENDDO
                 lr = INDEX(prt1,',')
                 line = '  '//'"'//prt1(lj:lr-1)//'"'//','
                 WRITE(io_out,FMT='(a)')TRIM(line)
              ELSE
                 EXIT ! leave do while
              ENDIF
           ENDDO
              line = prt1      ! new variable line,will be written  below
        ENDIF
        !write remainder of namelist :
        WRITE(io_out,FMT='(a)')TRIM(line)
     END DO

 10  CLOSE(UNIT=io_in, STATUS = 'DELETE') ! this is filename
     CLOSE(UNIT=io_out,STATUS = 'KEEP')    ! this is filenamec
     ! change filenamec to filename:
      call system("mv " // trim(filenamec) // " " // trim(filename))

     RETURN
2       lerrno = 356  + iomaxerr
        CALL terminate(lerrno,nlog)

  END SUBROUTINE  add_quotes_to_namelist





  SUBROUTINE P_Nfreya_output_name(outpt_time)
!---------------------------------------------------------------------------
! -- Generate name of statefile tah P_Nfreya will write as its output file
!---------------------------------------------------------------------------
        CHARACTER*14 record
        REAL(DP) outpt_time
        WRITE(record,FMT='(1pe14.6)')outpt_time
        statefile_output_name(1:LEN(statefile_output_name))=' '
        statefile_output_name = '"P_Nfreya_output_'//TRIM(ADJUSTL(record(1:LEN_TRIM(record))))//'.nc"'

     RETURN
  END SUBROUTINE P_Nfreya_output_name
  
#endif
  END MODULE Nfreya_namelist
