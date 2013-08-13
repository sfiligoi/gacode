
  MODULE P_Nfreya_rpc_interface

      USE nrtype,                                 ONLY : DP,I4B

      USE error_handler,                          ONLY : terminate,lerrno,iomaxerr

      USE io_gcnmp,                               ONLY : nlog

      USE string_util,                            ONLY : to_upper_case

      USE io_gcnmp,                               ONLY : switch_statefile_output,       &       
                                                          statefile_output_name
      USE Nfreya_namelist,                        ONLY : nubeam_namelist,nubeam_ufile

      USE common_constants,                       ONLY : zeroc,izero

      INTEGER(I4B), PARAMETER            :: std_file_size = 128

      CHARACTER(LEN = std_file_size)     :: P_Nfreya_remote_host,                &
                                            P_Nfreya_inpt_statefile,             &
                                            P_Nfreya_nml_filename,               &
                                            P_Nfreya_remote_dir,                 &
                                            P_Nfreya_outpt_statefile,            &
                                            P_Nfreya_rpc_client,                 &
                                            P_Nfreya_rpc_server,                 &
                                            MPI_hostfile
      CHARACTER(Len=8)                   :: remote_dir_dispatch,wipe_copy ! valid ARE 'wipe', 'save'
      REAL(DP)      P_Nfreya_call_time
      INTEGER(I4B)  P_Nfreya_nprocs, P_Nfreya_port_no,P_Nfreya_call_no

      DATA P_Nfreya_call_no,P_Nfreya_call_time /0, 0.0_DP/

  CONTAINS 
    
    SUBROUTINE setup_run_P_Nfreya(time,task)
      !----------------------------------------------------------------------
      ! -- Gateway to remote proceedure method of running P_Nfreya on a
      ! -- parallel machine
      !---------------------------------------------------------HSJ-4/22/2011-

      USE nrtype,                                  ONLY :  DP,I4B

      USE Nfreya_namelist,                         ONLY :  write_P_Nfreya_run_directives, &
                                                           P_Nfreya_output_name

      USE P_Nfreya_12_interface,                   ONLY :  P_Nfreya_run_directives,      &
                                                           P_Nfreya_run_directives_base, &
                                                           sent_ufile

      USE iterdbmd_gcnmp,                           ONLY : iterdb_file_name

      USE file_proc,                                ONLY : new_name
   
      USE iterdbmd,                                 ONLY : iterdb

      USE gcnmp_input,                              ONLY : write_iterdb_txt

      USE solcon_gcnmp,                             ONLY : time_max

      USE solcon,                                   ONLY : timmax

      IMPLICIT NONE

      REAL(DP) time
      INTEGER(I4B) task,fl
      CHARACTER*14 stime
      LOGICAL exists

      ! -----------------------------------------------------------------
      ! write the namelist file that describes
      ! the run to be done by P_Nfreya. This file does not change during the 
      ! course of a run and hence only needs to be written once
      ! ufile and nubeam type namelist for beam geometry specification
      ! is known to exist at this point and is readable because we checked
      ! it in sub init.
      ! ------------------------------------------------------------------
!      IF(task == 0)THEN
         WRITE(stime,FMT='(1pe14.6)')time 

         P_Nfreya_run_directives = P_Nfreya_run_directives_base(1:LEN_TRIM(P_Nfreya_run_directives_base)) &
              //'_'//ADJUSTL(stime(1:LEN_TRIM(stime))) ! file name for namelist run directives

         P_Nfreya_nml_filename = P_Nfreya_run_directives
         fl = LEN_TRIM(P_Nfreya_run_directives)

         IF(fl .LE. LEN(P_Nfreya_nml_filename))THEN 
            P_Nfreya_nml_filename(1:fl) = P_Nfreya_run_directives(1:fl)
         ELSE
            lerrno = 352 + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF
         time_max = 2.*timmax  ! HSJ 8/11/11
         CALL write_P_Nfreya_run_directives(time)               ! nfreya_namelist.f90

         ! --------------------------------------------------------------------
         ! write a state file at this time for use in P_Nfreya
         ! to simplify we only write netcdf state file here
         ! --------------------------------------------------------------------
         iterdb = 1 ! take control away from user, ensure that iterdb is written
         WRITE_ITERDB_TXT = .FALSE.

         CALL Statefile_proc                         ! Statefile_proc.f90

         fl = LEN_TRIM(new_name)
         P_Nfreya_inpt_statefile(1:LEN(P_Nfreya_inpt_statefile))=' '

         IF(fl .LE. LEN(P_Nfreya_inpt_statefile))THEN 

            P_Nfreya_inpt_statefile(1:fl) = new_name(1:fl)
         ELSE
            lerrno = 352 + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF
!         CALL P_Nfreya_output_name(time)      ! nfreya_namelist.f90

         P_Nfreya_outpt_statefile(1:LEN(P_Nfreya_outpt_statefile)) = ' '
         fl = LEN_TRIM(statefile_output_name)
         ! P_Nfreya_outpt_statefile is filename without the quotes required for namelist write
         ! in statefile_output_name(1:1) and  statefile_output_name(fl:fl) 
         P_Nfreya_outpt_statefile(1:fl) = statefile_output_name(2:fl-1)

      ! --------------------------------------------------------------------
      ! run   P_Nfreya remotely and wait for return.
      !---------------------------------------------------------------------

         CALL P_Nfreya_rpc_driver

      !-----------------------------------------------------------------------
      ! -- Check if P_Nfreya produced the output statefile
      !-----------------------------------------------------------------------
      INQUIRE(FILE = P_Nfreya_outpt_statefile ,EXIST = exists)
      IF(exists)THEN 
         PRINT *,'******** P_Nfreya output state file found by Onetwo:'
         PRINT *, ADJUSTL(P_Nfreya_outpt_statefile(1:LEN_TRIM(P_Nfreya_outpt_statefile)))
      ELSE
         PRINT *,'ERROR, expected statefile',P_Nfreya_outpt_statefile(1:LEN_TRIM(P_Nfreya_outpt_statefile))
         PRINT *,"       Was not returned, can't continue"
         lerrno = 350 + iomaxerr
         CALL terminate(lerrno,nlog)
      ENDIF

      !-----------------------------------------------------------------------
      ! -- read state file  Monte Carlo data that  P_Nfreya produced:
      !-----------------------------------------------------------------------
 
       CALL read_P_Nfreya_output      ! P_freya_rpc_interface.f90


      RETURN
    END SUBROUTINE setup_run_P_Nfreya





    SUBROUTINE P_Nfreya_rpc_driver
 !------------------------------------------------------------------------
 ! -- read configuration file P_Nfreya_rpc_config
 ! -- do remote proceedure call (CALL P_Nfreya_rpc)
 !--------------------------------------------------------------HSJ--------


      IMPLICIT NONE
      INTEGER, PARAMETER :: arglist = 1
      INTEGER, PARAMETER :: iot  = 10
      INTEGER, PARAMETER :: iot2 = 11
      INTEGER iloc,j,k,value,ierr ,system
      CHARACTER(len=std_file_size) carg(0:arglist)
      CHARACTER(len=std_file_size) inpt_copy,command,command1
      CHARACTER(len=LEN(P_Nfreya_outpt_statefile)) outpt_copy
      CHARACTER(len = 1) nmbr
      CHARACTER(len = 4) ext
      CHARACTER(len = 8) copy_nml
 
      LOGICAL inpt_txt,outpt_txt,exists,file_end


      NAMELIST /CONFIG/                             &
           P_Nfreya_nprocs,                         &
           P_Nfreya_port_no,                        &
           P_Nfreya_remote_dir,                     &
           P_Nfreya_remote_host,                    &
           P_Nfreya_rpc_client,                     &
           P_Nfreya_rpc_server,                     &
           MPI_hostfile,                            &
           remote_dir_dispatch


         remote_dir_dispatch = 'save'  ! directory on remote machine that P_nfreya will use
                                   ! will not be wiped clean before P_nfreya is run.
      carg(0:arglist)(:) = ' '
      ! read info to be passed to remote machine from configuration file.
      ! Note that some information is obtained dynamically:
      !       local working dir 
      !       remote location of P_Nfreya
      !       local host name
      ! NOTE THAT SERVER, P_NFREYA_SERVER.PY, MUST BE RUNNING ON REMOTE MACHINE
      ! and  listening on P_Nfreya_port_no (defined in P_Nfreya_rpc_config)
      OPEN (unit=iot,file='P_Nfreya_rpc_config', status='OLD',ERR = 20)
      READ(iot,nml=CONFIG,END = 2)
      CLOSE(UNIT=iot)


      !make a copy of the namelist file:
      !copy_nml = 'copy_nml'
      !WRITE(command,FMT='(a)')ADJUSTL(TRIM(P_Nfreya_nml_filename))
      !command1 ='cp '//command(1:LEN_TRIM(command))//' '//copy_nml
      !ierr = SYSTEM(command1)




      !get extension on input statefile:
      inpt_copy = to_upper_case(P_Nfreya_inpt_statefile)
      inpt_txt = .FALSE.
      IF(INDEX(inpt_copy,'.TXT') .GT. 0 )inpt_txt = .TRUE.
      ! Now see if the output file is in text or netcdf form
      outpt_copy = to_upper_case(P_Nfreya_outpt_statefile)
      outpt_txt = .FALSE.
      IF(INDEX(outpt_copy,'.TXT') .GT. 0 )outpt_txt = .TRUE.


      !switch extension name if necessary
      value = 0 ! do not use this option, input,output is of *.nc form for P_Nfreya
      IF(value == 1)THEN        ! a value of 1 means toggle extension
         iloc = INDEX(P_Nfreya_outpt_statefile,'.')
         IF(inpt_txt == .TRUE.)THEN ! input is .txt so switch output  to .nc
            ext = '.nc '
         ELSE                       ! input is  .nc so switch output to .txt
            ext = '.txt'
         ENDIF
         PRINT *,'inpt_txt =',inpt_txt
         P_Nfreya_outpt_statefile(1:iloc-1+4) = P_Nfreya_outpt_statefile(1:iloc-1)//ext
      ENDIF

      !-------------------------------------------------------------
      ! -- Do the remote proceedure call
      !-------------------------------------------------------------
           CALL P_Nfreya_rpc 


      RETURN

2        PRINT *,'unexpected end in config file namelist read'
         lerrno = 351 + iomaxerr
         CALL terminate(lerrno,nlog)
10    PRINT *,'Error, namelist file not found ',P_Nfreya_nml_filename
         lerrno = 352 + iomaxerr
         CALL terminate(lerrno,nlog)
20    PRINT *,'Error, input configuration file not found:',TRIM(ADJUSTL(carg(1)))
         lerrno = 353 + iomaxerr
         CALL terminate(lerrno,nlog)

    END SUBROUTINE P_Nfreya_rpc_driver



    SUBROUTINE P_Nfreya_rpc
      !---------------------------------------------------------------------------------------------
      !     NOTE: rpc ==> remote procedure call
      !     P_Nfreya_nprocs         = number  of processors to run P_Nfreya with
      !     P_Nfreya_remote_host    = name of remote or local machine which has
      !                               server listening for P_Nfreya run commands
      !     P_Nfreya_remote_dir       is directory on server in which P_Nfreya will be run
      !     P_Nfreya_dir              is directory on server where executable P_Nfreya
      !                               will be found (must be member of P_Nfreya_paths)
      !     local_dir                 is directory on local machine to which output 
      !                               from P_Nfreya will be shipped
      !     run_dirf                  is name of run directives file (namelist file) 
      !                               that P_Nfreya will read
      !     P_Nfreya_outpt_statefile  Name of state file  that P_Nfreya will produce on
      !                               the server side (and to be shipped back to the client)
      !     remote_dir_dispatch           remove (='wipe'  (or not = 'save')any files in 
      !                               remote directory before transferring files and 
      !                               running P_Nfreya
      !---------------------------------------------------------------------------------------------
      USE  ext_prog_info,                        ONLY : P_Nfreya_hosts,P_Nfreya_paths

      USE  transp,                               ONLY : nubeam_namelist_cpy

      USE P_nfreya_12_interface,                 ONLY : sent_ufile


      !rpc is set up in Python (could  do in Fortran with INTEL RPC LIB)
      INTEGER ierr,ISHELL,j,js,GETCWD,not,ncrt,ks,ke
      INTEGER hostnm,system
      CHARACTER(LEN = 6*std_file_size)  command,command1
      CHARACTER*8 npc,port
      CHARACTER(len=LEN(P_Nfreya_remote_host)) lname
      CHARACTER(LEN = std_file_size) local_dir,P_Nfreya_dir
      CHARACTER(len = 256) loc_host
      CHARACTER(len = LEN(nubeam_ufile))    nubeam_ufile_cpy
      CHARACTER*1 ufs

      ncrt = 6
      nout = 8


      lname = ADJUSTL(to_upper_case(P_Nfreya_remote_host))
      ierr  = GETCWD(local_dir)          ! sets  local dir to cwd
      ierr  = hostnm(loc_host)           ! 

      !hostnm has filled the tail end of loc_host with
      !disk hash which does not equate to blanks.
      !hence we have to do the following:

      command= TRIM(ADJUSTL(loc_host)) 

      loc_host(1:LEN_TRIM(command))=command(1:LEN_TRIM(command))
 
      js=INDEX(loc_host,'.com')

      IF(js .NE. 0)THEN ! ifort does not return fully qualified name ??
         js =js+4
      ELSE
         js =LEN_TRIM(loc_host)
      ENDIF
      loc_host(js:LEN(loc_host))=''

      !----------------------------------------------------------------
      ! NOTE: Fully qualified name of P_Nfreya executable is set in 
      ! ext_prog_info.f90
      ! For example running P_Nfreya on lohan1 through a remote
      ! procedure call from Onetwo assigns this executable:
      ! '/p/linux/onetwo/source/onetwo/nfreya_module/P_Nfreya_pgf_l1'
      !----------------------------------------------------------------
      DO  j=1,SIZE(P_Nfreya_hosts)
         IF(INDEX(P_Nfreya_hosts(j),TRIM(lname)) >0) THEN
            P_Nfreya_dir  = TRIM(P_Nfreya_paths(j))
            EXIT
         ENDIF
         IF( j == SIZE(P_Nfreya_hosts))THEN
            !WRITE(nout,1)TRIM(lname)
            WRITE(ncrt,1)TRIM(lname)
            CALL STOP('P_Nfreya_rpc, invalid host name',1)
1           FORMAT(' ERROR,host ',a,' does not run P_Nfreya')
         ENDIF
      ENDDO
       ufs = 'F'
       if(sent_ufile) ufs = 'T'


      WRITE(npc,FMT='(i5)')P_Nfreya_nprocs
      WRITE(npc,FMT='(i5,a)')P_Nfreya_nprocs,ufs
      WRITE(port,FMT='(i8)')P_Nfreya_port_no
      PRINT *,' npc,port  =',npc,port 
      PRINT *,' remote_dir =',TRIM(ADJUSTL(P_Nfreya_remote_dir)) 
      PRINT *,' local_dir =',TRIM(ADJUSTL(local_dir))
      PRINT *,' P_Nfreya_dir =',TRIM(ADJUSTL(P_Nfreya_dir)) 
      PRINT *,' host = ',TRIM(ADJUSTL(P_Nfreya_remote_host))
      PRINT *,' P_Nfreya output state file =',P_Nfreya_outpt_statefile
      PRINT *,' remote_dir_dispatch =',TRIM(ADJUSTL(remote_dir_dispatch))


      ierr = 0
      command(1:len(command))=''
      wipe_copy =  to_upper_case( remote_dir_dispatch )
      IF(wipe_copy  ==  'WIPE')THEN

         command = TRIM(ADJUSTL(P_Nfreya_rpc_client))           &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(npc))                           &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(port))                          &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(P_Nfreya_remote_dir))           &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(P_Nfreya_remote_host))          &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(MPI_hostfile))                  &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(P_Nfreya_inpt_statefile))       &  !  the rest are determined
           //' '// TRIM(ADJUSTL(P_Nfreya_outpt_statefile))      &  !  in Onetwo
           //' '// TRIM(ADJUSTL(P_Nfreya_nml_filename))         &
           //' '// TRIM(ADJUSTL(nubeam_namelist_cpy))           &
           //' '// TRIM(ADJUSTL('dummy'))                       &  !ufile copy not required here
           //' '// TRIM(ADJUSTL(local_dir))                     &
           //' '// TRIM(ADJUSTL(P_Nfreya_dir))                  &
           //' '// TRIM(ADJUSTL(wipe_copy(1:LEN_TRIM(wipe_copy)))) &  ! from P_Nfreya_rpc_config
           //' '// TRIM(ADJUSTL(loc_host(1:js))) 
!           //' '// ufs
! note command doesnt work unless loc_host is last item               
         PRINT *," Clearing remote directory before use with command"
         PRINT *, command
!----------------------------------------------------------------------------------------------
! the CALL system should execute this:
!/bin/sh -c 
!----------------------------------------------------------------------------------------------
         ierr = SYSTEM(command1)
!         CALL SYSTEM(command,ierr)
      ENDIF


      IF (ierr .NE.   0 )THEN
         WRITE(ncrt,FMT ='("Sub P_Nfreya_rpc failed with command",a)')command
         CALL STOP(' Error in rpc call to  P_Nfreya',-1)
      ENDIF


      ! setup for next command:
      command(1:len(command))=' '
      !remove  quotes from nubeam_namelist, nubeam_ufile before passing command:
      !nubeam_namelist_cpy(:) = ' '
      nubeam_ufile_cpy(:) = ' '
      ks  = izero
      ke  = izero
      ks  = INDEX(nubeam_namelist,'"',.FALSE.)
      ke  = INDEX(nubeam_namelist,'"',.TRUE.)
      IF(ks .GE. 1)THEN
         IF(ke .GT. ks+1)THEN
            ks = ks+1
            ke = ke-1
            IF(nubeam_namelist(ks:ks+1) == './')ks = ks+2
            !nubeam_namelist_cpy(1:ke-ks +1) = nubeam_namelist(ks:ke) 
         ELSE
            lerrno = 354 + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF
      ENDIF

      ks  = izero
      ke  = izero
      ks  = INDEX(nubeam_ufile,'"',.FALSE.)
      ke  = INDEX(nubeam_ufile,'"',.TRUE.)
      IF(ks .GE. 1)THEN
         IF(ke .GT. ks+1)THEN
            ks = ks+1
            ke = ke-1
            nubeam_ufile_cpy(1:ke-ks +1) = nubeam_ufile(ks:ke) 
         ELSE
            lerrno = 354 + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF
      ENDIF
  

      !---------------------------------------------------------------------------------
      ! create a command  to run the P_nfreya client code with arguments .
      ! The argument list consists of itmes known only in Onetwo plus
      ! items needed for the rpc paradigm to work. These additional items
      ! are not necessary in Onetwo but are brought in by reading the
      ! the  P_Nfreya_rpc_config file. The config file could be read
      ! in the client code instead but the client name,P_Nfreya_rpc_client,
      ! must be known in Onetwo. It is expected that this client is
      ! different for different users so it is not included in ext_prog_info.f90.
      ! Consequently we read it from P_Nfreya_rpc_config. Having had to read this file
      ! to get one piece of information, I opted to actually get all
      ! necessary information at this read and then pass it to the client.
      ! Hence the client does not need to read P_Nfreya_rpc_config again.
      !---------------------------------------------------------HSJ---------------------
      wipe_copy ='SAVE'                                            ! must be upper case
      command = TRIM(ADJUSTL(P_Nfreya_rpc_client))              &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(npc))                           &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(port))                          &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(P_Nfreya_remote_dir))           &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(P_Nfreya_remote_host))          &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(MPI_hostfile))                  &  !  from P_Nfreya_rpc_config 
           //' '// TRIM(ADJUSTL(P_Nfreya_inpt_statefile))       &  !  the rest are determined
           //' '// TRIM(ADJUSTL(P_Nfreya_outpt_statefile))      &  !  in Onetwo
           //' '// TRIM(ADJUSTL(P_Nfreya_nml_filename))         &
           //' '// TRIM(ADJUSTL(nubeam_namelist_cpy))           &
           //' '// TRIM(ADJUSTL(nubeam_ufile_cpy))              &
           //' '// TRIM(ADJUSTL(local_dir))                     &
           //' '// TRIM(ADJUSTL(P_Nfreya_dir))                  &
           //' '// TRIM(ADJUSTL(wipe_copy(1:LEN_TRIM(wipe_copy)))) &  ! from P_Nfreya_rpc_config
           //' '// TRIM(ADJUSTL(loc_host(1:js)))   

      command =  TRIM(ADJUSTL(command))
! note command doesnt work unless loc_host is last item

      PRINT *," Onetwo making remote proceedure call to P_Nfreya with "
      PRINT *, command

      !---------------------------------------------------------------------------
      !   execute command:
      !----------------------------------------------------------------------------

      ierr = SYSTEM(command)
!       CALL SYSTEM(command,ierr)
        sent_ufile = .TRUE.  ! subsequent calls will not send ufile

      IF (ierr .NE.   0 )THEN
         WRITE(ncrt,10)ierr
         !          WRITE(nout,10)ierr
         CALL STOP(' Error in rpc call to  P_Nfreya',1)
10       FORMAT(' Sub P_Nfreya_rpc returned SYSTEM  error = ',i5)
      ENDIF

    END SUBROUTINE P_Nfreya_rpc



    SUBROUTINE read_P_Nfreya_output
!---------------------------------------------------------------------------
!   Read statefile produced by P_Nfreya,P_Nfreya_outpt_statefile
!---------------------------------------------------------------------------

      USE iterdbmd,                         ONLY : statefile_name,statefile_type
      USE gcnmp_input,                      ONLY : write_iterdb_txt,switch_iterdb_output
      statefile_name = TRIM(P_Nfreya_outpt_statefile)

      CALL read_statefile   ! Statefile_proc.f90
      switch_iterdb_output  = 0  ! make sure next file to be written is
                                 ! in netcdf form
    RETURN
    END SUBROUTINE read_P_Nfreya_output



  END MODULE P_Nfreya_rpc_interface
