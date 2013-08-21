
   MODULE gcnmp_interface
! ----------------------------------------------------------------------
! --- 
! --- Create namelist input file for use in GCNMP code
! ---
! ----------------------------------------------------------HSJ-06/29/07



    USE nrtype,                ONLY : I4B,DP

 
    USE gcnmp_input,           ONLY : gcnmp_nprocs,                     &
                                      write_iterdb_txt,                 &
                                      iterdbf,                          &
                                      swioin => switch_iterdb_output,   &
                                      eq_split_gcnmp,itran_gcnmp,       &
                                      save_incr_restart_file_gcnmp,     &
                                      gcnmp_macro_dt,dtmin_gcnmp,       &
                                      gcnmp_nml_filename,gcnmp_host,    &
                                      gcnmp_remote_dir,                 &
                                      gcnmp_iterdb_filename
    USE io_gcnmp,              ONLY : switch_iterdb_output,rpc,         &
                                      iterdb_12_input_file

    USE iterdbmd,              ONLY : irwflag_12 => irwflag,            &
                                      iterdbfilename      

    USE iterdbmd_gcnmp,        ONLY : irwflag_gcnmp => irwflag,&
                                      iterdb_file_name


    USE solcon_gcnmp,          ONLY : std_file_size,use_glf23,          &
                                      single_density_simulation

    USE GCNMP_namelist 

    USE  set_12_gcnmp_vars,    ONLY : set_onetwo_vars, set_gcnmp_vars


    USE ions_gcnmp,            ONLY :  name_size

    IMPLICIT NONE 



     Contains


       SUBROUTINE write_gcnmp_namelist

         USE bc_values_gcnmp,    ONLY : u_vloop_bc_12 => u_vloop_bc



         USE numbrs,             ONLY : nprim,nimp

         USE solcon,             ONLY : time_12  => time,                   &
                                        dtmax_12 => dtmax,                  &
                                        dt_12    => dt,                     &
                                        nmax,timmax,                        &
                                        theta_12 => theta,                  &
                                        steady_state_12 => steady_state

         USE soln2d,             ONLY : rbsaxis_12 => rbsaxis

         USE glf23,              ONLY : include_glf ,                       &
                                        x_alpha_glf_12 => x_alpha_glf,      &
                                        exbmult_glf_12 => exbmult_glf,      &
                                        limp_glf_12 => limp_glf,            &
                                        jroot_iglf,                         &
                                        ibtflag_glf_12 => ibtflag_glf,      &
                                        irotstab_12 => irotstab,            &
                                        jeigen_glf_12 => jeigen_iglf,       &
                                        i_delay_12 => i_delay,              &
                                        iglf_12_chie    => iglf_chie,       &
                                        iglf_12_chii    => iglf_chii,       &
                                        iglf_12_chiv    => iglf_chiv,       &
                                        iglf_12_d       => iglf_d           

         USE io,                ONLY : nout,ncrt


         USE nonlin,             ONLY : switch_method_12  => switch_method, &
                                        tot_iters_max_12  => tot_iters_max, &
                                        non_lin_method_12 => non_lin_method,&
                                        fdigits_12  => fdigits,             &
                                        gradtol_12  => gradtol,             &
                                        steptol_12  => steptol,             &
                                        bandwidth_12 => bandwidth,          &
                                        jacobian_type_12 => jacobian_type,  &
                                        fvectol_12 => fvectol,              &
                                        ssqrmin_12 => ssqrmin,              &
                                        maxfev_12  => maxfev,               &
                                        iters_freeze_12 => iters_freeze,    &
                                        freeze_type_12 => freeze_type
                                        

        USE sourc,               ONLY : wdelt_12 => wdelt
                                        
        USE aid_newton,          ONLY : jac_skip_12 => jac_skip    
        
        USE verbose,             ONLY : newtonvb_12 => newtonvb

        USE tfact,               ONLY : wneo,set_chie_chii_12 => set_chie_chii, &
                                        jboot

        USE pelcom,              ONLY : nampel


        LOGICAL exists,opened
        INTEGER(I4B) iounit,newtonvb
        REAL(DP) pos_def

! -----------------------------------------------------------------------
! ---   Load gcnmp variables from 12 variables:
! -----------------------------------------------------------------------

!       set default values:

        pos_def = 10._DP*EPSILON(0.0_DP)
        itte = 0  ; itti = 0 ;itxj = 0 
        itw  = 0 ; itenpd(:) = 0 ; itenid(:) =0
        steady_state = 0 
        pellet_name  = "'"//TRIM(nampel)//"'"


!       load transport flags
        itenpd(1:nprim) = itran_gcnmp(1:nprim)
        itenid(1:nimp) = itran_gcnmp(nprim+1:nprim+nimp)
        itte = itran_gcnmp(nprim+nimp+1)
        itti = itran_gcnmp(nprim+nimp+2)
        itxj = itran_gcnmp(nprim+nimp+3)
        itw  = izero
        IF(include_glf > 0) itw  = itran_gcnmp(nprim+nimp+4)

!       load gcnmp flags:
        IF(include_glf > 0)THEN
           !currently no distinction is made on primary and impurity
           !transport models. Either all are glf23 or none are.
           !Note that this does not say which ions are transported and which ones
           !are run in simulation. (That is done with itenpd,itenid).
           !This simply says that if an ion is transported then it will
           !use glf23 for the diffusivity:
           IF(ABS(iglf_12_d) .GT. pos_def)use_glf23(1) =1
           !TE:
           IF(ABS(iglf_12_chie) .GT. pos_def)use_glf23(2) =1
           !TI:
           IF(ABS(iglf_12_chii) .GT. pos_def)use_glf23(3) =1
           !use_glf23(4) is current density slot that remains unused
           use_glf23(4) = izero
           !toroidal rotation:
           IF(ABS(iglf_12_chiv) .GT. pos_def)use_glf23(5) =1
        ENDIF

        time0 = time_12                 !time_12  is current time in 12
        time_max = MIN(timmax,time_12 + gcnmp_macro_dt)
        IF(ABS(time_max-timmax) .LT. pos_def)time_max = timmax
        dtmax = dtmax_12 
        dtmin = dtmin_gcnmp ! new dtmin in onetwo not in modules 
        dtmin = MIN(dtmin,time_max-time0)
        dtstart = dt_12
        dtstart = MIN(dtstart,dtmin)
        IF(SIZE(wneo,1) == SIZE(neocl_mult,1) .AND. &
                                  SIZE(wneo,2) == SIZE(neocl_mult,2))THEN
           neocl_mult(:,:)          = wneo(:,:)
        ELSE
           WRITE(nout,1)
           WRITE(ncrt,1)
1          FORMAT(2x,'Error detected in sub write_gcnmp_namelist',/, &
                  2x,'Incompatibility in defined # transport variables')
           CALL STOP('error in onetwo/Gcnmp compatitibility settings',1)
        ENDIF

        max_steps                = nmax 
        freeze_type              = freeze_type_12
        switch_method            = switch_method_12               
        tot_iters_max            = tot_iters_max_12
        non_lin_method           = non_lin_method_12
        fdigits                  = fdigitS_12
        gradtol                  = gradtol_12     
        steptol                  = steptol_12
        bandwidth                = bandwidth_12
        jacobian_type            = jacobian_type_12
        fvectol                  = fvectol_12       
        ssqrmin                  = ssqrmin_12
        maxfev                   = maxfev_12
        IF( steady_state_12 > 0.1) steady_state = 1             
        conv_skip                = 0_I4B
        wrt_nwt                  = 0_I4B
        newtonvb                 = newtonvb_12
        iters_freeze             = iters_freeze_12
        random_pert              = 3_I4B
        jac_skip                 = jac_skip_12
        asym_time                = time_max + 1. !turns off this option
        asym_dt_max              = dtmax         !not used (see  asym_time) 
        create_plot_file         = 0_I4B
        use_constd(:)            = 0_I4B
        constd_values(:)         = 0.0_DP
        wdelt                    = wdelt_12
        consts_values(:)         = 0.0_DP
        use_consts(:)            = 0_I4B
        use_avg_glf_chi          = 0_I4B
        use_avg_chi              = 0_I4B
        rbsaxis                  = rbsaxis_12
        vloopvb                  = 0_I4B
        u_vloop_bc               = u_vloop_bc_12
        irfc_mult                = 1.0_DP
        ibcur_mult               = 1.0_DP
        set_cap1                 = 0_I4B
        use_constconv(:)         = 0_I4B
        constconv_values(:)      = 0.0_DP
        set_chie_chii            = set_chie_chii_12
        jbootstrap               = "'Sauter'"
        IF(jboot == 110)jbootstrap = "'Sauter_zeff'" !old way uses zeff

        use_mask                 = 0_I4B
        exbmult_glf              = exbmult_glf_12
        x_alpha_glf              = x_alpha_glf_12
        theta                    = theta_12
        jroot_glf                = jroot_iglf
        limp_glf                 = limp_glf_12
        ibtflag_glf              = ibtflag_glf_12
        irotstab                 = irotstab_12
        dv_method                = 0_I4B
        spline_gsc               = 0_I4B          ! spline/dont spline profiles requires dv_method =1
        dv_delt                  = 0.01_DP        ! for dv_method only
        itte_dv                  = 0_I4B
        itti_dv                  = 0_I4B
        itangrot_dv              = 0_I4B
        itenp_dv                 = 0_I4B
        itene_dv                 = 0_I4B
        jeigen_glf               = jeigen_glf_12
        i_delay                  = i_delay_12
        nsteps_plot              = 0_I4B
        plot_time_incr           = 0.0_DP
        save_incr_plot_file      = 0_I4B
        nsteps_restart           = 0_I4B
        create_restart_file      = 0_I4B
        restart_time_incr        = 0.0_DP    
        parallel_model           = 0_I4B
        max_bwcycle              = 3_I4B
        max_bandwidth            = 5_I4B
        min_bandwidth            = 0_I4B
        cycle_bandwidth          = .FALSE.
        mul_flux(:)              = 0.0_DP
        mul_den(:)               = 1.0_DP
        gamma_bc(:)              = 0.0_DP 
        switch_iterdb_output     = swioin
        axis_bc_type             = 0_I4B
        cgd                      = 0.0_DP
        cgexp                    = 1.0_DP
        crit_grad                = 0.0_DP
        t_delay_glf23            = 0.0_DP


        !when called from onetwo we always use rpc=1 to indicate that
        !the iterdb output file created by gcnmp, which must be read
        !by onetwo once gcnmp has finished, is to be named
        !iterdb_12_input_file
        rpc                      = 1
        IF(write_iterdb_txt .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in text mode and the gcnmp
           !output file is in netcdf form
             iterdb_12_input_file     = "'iterdb_12_input.nc'" ! also set below
        ELSEIF(.NOT. write_iterdb_txt .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in netcdf mode and the gcnmp
           !output file is in text form
           iterdb_12_input_file     = "'iterdb_12_input.txt'" ! also set below
        ELSEIF(write_iterdb_txt .AND. switch_iterdb_output ==0)THEN
           !here the gcnmp input and output files are in text form
           iterdb_12_input_file     = "'iterdb_12_input.txt'" ! also set below
        ELSEIF(.NOT. write_iterdb_txt .AND. switch_iterdb_output == 0)THEN
           !here the gcnmp input and output files are in netcdf  form
           iterdb_12_input_file     = "'iterdb_12_input.nc'" ! also set below
        ENDIF



        !NOTE: There are some important input variables introduced in gcnmp
        !      that have no counterpart in Onetwo:
        !      eq_split and save_incr_restart_file are gcnmp code variables
        !      they are set to eq_split_gcnmp and save_incr_restart_file_gcnmp
        !      because we use the _gcnmp extension in inone to clearly indicate
        !      that these variables are for gcnmp only and will not be used in
        !      Onetwo otherwise.
        eq_split                 = eq_split_gcnmp        
        save_incr_restart_file   = save_incr_restart_file_gcnmp 




! ---------------------------------------------------------------------------------
! --- Finally create the namelist file:
! ---------------------------------------------------------------------------------
        INQUIRE(FILE = gcnmp_nml_filename, NUMBER=iounit,EXIST = exists,OPENED=opened )
        IF(exists)THEN
           IF(opened)CLOSE(UNIT=iounit)
           CALL DESTROY (gcnmp_nml_filename )
        ENDIF
        iounit = 200
        CALL  getioun(iounit,iounit)
        OPEN(UNIT=iounit,  file=gcnmp_nml_filename,  status='new')
        WRITE(iounit,NML = run_data)
        
        CALL giveupus(iounit)
        CLOSE(UNIT=iounit)
        !the value of iterdb_12_input_file was set to "'iterdb_12_input_file'"
        !above so that teh namelist output would be correct. We now need
        !to strip away the double quotes so iterdb_12_input_file can be used 
        !elsewhere:
        IF(write_iterdb_txt .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in text mode and the gcnmp
           !output file is in netcdf form
           iterdb_12_input_file     = 'iterdb_12_input.nc'
        ELSEIF(.NOT. write_iterdb_txt .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in netcdf mode and the gcnmp
           !output file is in text form
           iterdb_12_input_file     = 'iterdb_12_input.txt'
        ELSEIF(write_iterdb_txt .AND. switch_iterdb_output ==0)THEN
           !here the gcnmp input and output files are in text form
           iterdb_12_input_file     = 'iterdb_12_input.txt'
        ELSEIF(.NOT. write_iterdb_txt .AND. switch_iterdb_output == 0)THEN
           !here the gcnmp input and output files are in netcdf  form
           iterdb_12_input_file     = 'iterdb_12_input.nc'
        ENDIF
        RETURN

       END        SUBROUTINE write_gcnmP_namelist







      SUBROUTINE gcnmp_driver
      USE io,                    ONLY : nout,ncrt

!-----------------------------------------------------------------------------     
! create the namelist that gcnmp will read:
!------------------------------------------------------------------------------
        CALL write_gcnmp_namelist 






!------------------------------------------------------------------------------
! create the startup  statefile for gcnmp:
!-------------------------------------------------------------------------------
        IF(write_iterdb_txt)THEN
           irwflag_gcnmp = 0 
           !NOTE: iterdbfilename = 'iterdb',set in cray101.f,  if gcnmp is not used
           !check for extension of .txt in gcnmp_iterdb_filename:
            IF(INDEX(gcnmp_iterdb_filename,'.txt') == 0)THEN
              iterdb_file_name = TRIM(ADJUSTL(gcnmp_iterdb_filename))//'.txt'
              gcnmp_iterdb_filename = TRIM(ADJUSTL(iterdb_file_name))
           ENDIF
           iterdbfilename = TRIM(ADJUSTL(gcnmp_iterdb_filename))
           ! print *,'gcnmp_iterdb_filename =',gcnmp_iterdb_filename
           CALL set_gcnmp_vars
           CALL iter_dbase_txt       ! writes text file iterdbfilename
           iterdbf = iterdbfilename(1:LEN_TRIM(iterdbfilename))
        ELSE
           !check for extension of .nc in gcnmp_iterdb_filename
           IF(INDEX(gcnmp_iterdb_filename,'.nc') == 0 )THEN
              iterdb_file_name = TRIM(ADJUSTL(gcnmp_iterdb_filename))//'.nc'
              gcnmp_iterdb_filename = TRIM(ADJUSTL(iterdb_file_name))
           ENDIF
           irwflag_gcnmp =0
           CALL set_gcnmp_vars  ! netcdf I/O uses gcnmp variables  
           iterdb_file_name = TRIM(ADJUSTL(gcnmp_iterdb_filename))
           CALL iter_dbase_nc          ! writes netcdf  file iterdb_file_name
           iterdbf = iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
        ENDIF





!-----------------------------------------------------------------------
!ship statefile  and namelist file to remote destination and run code :
!------------------------------------------------------------------------
        
        CALL gcnmp_rpc                 ! remote proceedure call to gcnmp




        print *,'------------Returned control to Onetwo----------------'
!--------------------------------------------------------------------------
        !gcnmp has run and a new statefile exist in cwd. read it back
        !into Onetwo:
!--------------------------------------------------------------------------
        IF(write_iterdb_txt .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in text mode and the gcnmp
           !output file is in netcdf form
           irwflag_gcnmp =1
           iterdb_file_name = TRIM(ADJUSTL(iterdb_12_input_file))
           CALL iter_dbase_nc             ! reads gcnmp output in netcdf form
           CALL set_onetwo_vars           ! load 12 vars with units conversion 
        ELSEIF(.NOT. write_iterdb_txt .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in netcdf mode and the gcnmp
           !output file is in text form
           irwflag_gcnmp = 1
           iterdb_file_name = TRIM(ADJUSTL(iterdb_12_input_file))
           CALL iter_dbase_txt          !reads  gcnmp output in text form
           !this option needs  CALL set_onetwo_vars to handle npsi difference
           !in onetwo versus gcnmp
           CALL set_onetwo_vars           ! load 12 vars with units conversion 
           !CALL STOP(' GCNMP interface not coded ',1)
        ELSEIF(write_iterdb_txt .AND. switch_iterdb_output ==0)THEN
           !here the gcnmp input and output files are in text form
           irwflag_gcnmp = 1
           iterdb_file_name = TRIM(ADJUSTL(iterdb_12_input_file))
           CALL iter_dbase_txt            !reads  gcnmp output in text form
           CALL set_onetwo_vars           ! load 12 vars with units conversion 
           !CALL STOP(' GCNMP interface not coded ',2)
        ELSEIF(.NOT. write_iterdb_txt .AND. switch_iterdb_output == 0)THEN
           !here the gcnmp input and output files are in netcdf  form
           irwflag_gcnmp= 1
           iterdb_file_name = TRIM(ADJUSTL(iterdb_12_input_file))
           CALL iter_dbase_nc             ! reads gcnmp output in netcdf form
           CALL set_onetwo_vars 
        ELSE
           write(ncrt,10)write_iterdb_txt,switch_iterdb_output
           write(nout,10)write_iterdb_txt,switch_iterdb_output
           CALL STOP(' Error in reading  gcnmp output',1)
10         FORMAT(' Sub gcnmp_driver flags = ',L,2x,I5)
        ENDIF

      RETURN

      END SUBROUTINE gcnmp_driver







      SUBROUTINE gcnmp_rpc

        USE string_util,                   ONLY : to_upper_case
        USE ext_prog_info,                 ONLY : gcnmp_hosts,gcnmp_paths
        USE io,                            ONLY : nout,ncrt
        !rpc is set up in Python (could  do in Fortran with INTEL RPC LIB)
        INTEGER(I4B) ierr,ISHELL,j,js,GETCWD
        INTEGER hostnm,system
        CHARACTER(LEN = 2*std_file_size)  command
        CHARACTER*8 npc
        CHARACTER(len=LEN(gcnmp_host)) lname
        CHARACTER(LEN = std_file_size) remote_dir,local_dir,gcnmp_dir
        CHARACTER(len = 256) loc_host,gcnmp_loc

   !------------------------------------------------------------------------
        !gcnmp_nprocs  = number  of processors to run gcnmp with
        !gcnmp_host    = name of remote or local machine which has
        !                server listening for gcnmp run commands
        !gcnmp_remote_dir    is directory on server in which gcnmp will be run
        !gcnmp_dir     is directory on server where executable gcnmp 
        !              will be found (must be member of gcnmp_paths)
        !local_dir     is directory on local machine to which output 
        !              from gcnmp will be shipped
        !iterdbf       is iterdb file that gcnmp will use for startup
        !run_dirf      is name of run directives file (namelist file) 
        !              that gcnmp will read
    !------------------------------------------------------------------------

        lname = ADJUSTL(to_upper_case(gcnmp_host))
        ierr  = GETCWD(local_dir) ! sets  local dir to cwd
        ierr  = hostnm(loc_host)
        !hostnm has filled the tail end of loc_host with
        !disk hash which does not equate to blanks.
        !hence we have to do the following:
        command= TRIM(ADJUSTL(loc_host)) 
        loc_host(1:LEN_TRIM(command))=command(1:LEN_TRIM(command))
        js=INDEX(loc_host,'.com')
        js =js+4
        loc_host(js:LEN(loc_host))=' '
        remote_dir ='usr/tmp'
        DO  j=1,SIZE(gcnmp_hosts)
           IF(INDEX(gcnmp_hosts(j),TRIM(lname)) >0) THEN
              gcnmp_dir  = TRIM(gcnmp_paths(j))
              EXIT
           ENDIF
           IF( j == SIZE(gcnmp_hosts))THEN
               WRITE(nout,1)TRIM(lname)
               WRITE(ncrt,1)TRIM(lname)
               CALL STOP('gcnmp_rpc, invalid host name',1)
 1             FORMAT(2x,a,'ERROR,host does not run gcnmp')
           ENDIF
        ENDDO



        WRITE(npc,FMT='(i5)')gcnmp_nprocs

!         print *,'remote_dir =',TRIM(ADJUSTL(gcnmp_remote_dir)) 
!         print *,'local_dir =',TRIM(ADJUSTL(local_dir))
!         print *,'gcnmp-dir =',TRIM(ADJUSTL(gcnmp_dir)) 
!         print *,' host = ',TRIM(ADJUSTL(gcnmp_host))
        gcnmp_loc = '/p/linux/onetwo/pgf90/gcnmp_client.py'
        command = TRIM(ADJUSTL(gcnmp_loc))     &
             //' '// TRIM(ADJUSTL(npc))                                      &
             //' '// TRIM(ADJUSTL(iterdbf))                                  &
             //' '// TRIM(ADJUSTL(gcnmp_nml_filename))                       &
             //' '// TRIM(ADJUSTL(gcnmp_remote_dir))                         &
             //' '// TRIM(ADJUSTL(local_dir))                                &
             //' '// TRIM(ADJUSTL(gcnmp_dir))                                &
             //' '// TRIM(ADJUSTL(gcnmp_host))                               &
             //' '// TRIM(ADJUSTL(loc_host(1:js)))                           &
             //' '// TRIM(ADJUSTL(iterdb_12_input_file))
        command =  TRIM(ADJUSTL(command))
        print *,command

!---------------------------------------------------------------------------
!   execute command:
!----------------------------------------------------------------------------
           ierr = SYSTEM(command)


        IF (ierr .NE.   0 )THEN
           write(ncrt,10)ierr
           write(nout,10)ierr
           CALL STOP(' Error in spawning gcnmp',1)
10         FORMAT(' Sub gcnmp_rpc returned ISHELL error = ',i5)
        ENDIF
      RETURN
      END SUBROUTINE gcnmp_rpc  


   END MODULE gcnmp_interface
 
