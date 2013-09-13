 SUBROUTINE get_beam_data(name_beam)
! --------------------------------------------------------------------------
! read beam geometry info from file beam_data_namelist
! read beam data from 2d ufile (beam_data_ufile)
! uses NTCC ufiles modules-libraries
! -------------------------------------------------------------HSJ 11/13/03----
    USE transp,                               ONLY :   use_nubeam,beam_data_ufile,     &
                                                       beam_data_namelist,             &
                                                       nubeam0_dt,nubeam_dt,d_fast_ion,&
                                                       use_ufile,nubeam_namelist_cpy,  & 
                                                       nubeam_version
    USE io,                                   ONLY : ncrt,nitre
    USE fast_ion_diffusion ,                  ONLY : fi_init
    USE beam_structure,                       ONLY : neutral_beam, beam_structure_allocate
    USE string_util,                          ONLY : strip_path
    USE solcon,                               ONLY : time0
    USE numbrs,                               ONLY : nprim,nimp
    USE string_util,                          ONLY : to_upper_case1

    IMPLICIT NONE


    INTEGER,PARAMETER :: NSCMAX= 10         ! maximum number of scalars; 
                                            ! array dimension)
    INTEGER ndifbe
    REAL*4 SCVAL (NSCMAX)                   ! (scalar data array)
    CHARACTER*4 TDEV                        ! (tokamak dev, eg DIII-D,etc
    CHARACTER*10 SHDATE                     ! (shot date)
    CHARACTER*10 SCLAB(3,nscmax)            ! (scalar labels array)
    CHARACTER(len = LEN_TRIM(beam_data_namelist)+3) :: beam_data_namelist_copy
    CHARACTER(len = LEN_TRIM(beam_data_namelist)+3) :: longpath,longfile
    REAL*4,DIMENSION(:),ALLOCATABLE  :: X,y   !arrays to hold x,y)
    REAL*4,DIMENSION(:,:),ALLOCATABLE :: F    !or REAL F(NFMAX), array to hold f
    CHARACTER*10 XLAB(3), YLAB(3), FLAB(3)  ! on the VAX (label arrays 
                                            !  for X, Y, and  F) 
                                            !(NOTE CHARACTER*30 XLAB,
                                            !FLAB will also work
    TYPE(neutral_beam)name_beam
    INTEGER ier,ulun,nlun,iostat,limit,ishot,iflag,iproc,nx,ny,nsc,itry, &
            NXMAX,NYMAX,ierr,ngmax,nrhix,nbeam,nfast
    LOGICAL   nexists, uexists,nopened,read_com,uopened,read_ufile
    CHARACTER *64 path,prefix,extension     !64 is limit in ufiles routines !JMP
    INTEGER j,ntemp,nml_size
 

   INTERFACE local
      SUBROUTINE  set_beam_data(F,x,y,nx,ny)
          REAL*4 , INTENT (IN),DIMENSION(:,:) ::     F
          REAL*4 , INTENT (IN),DIMENSION(:) :: x,y
          INTEGER, INTENT (IN) :: nx,ny
      END  SUBROUTINE  set_beam_data
     SUBROUTINE sort_unique(a,km)
       REAL *8, INTENT(INOUT),DIMENSION(:) :: a
       INTEGER, INTENT(OUT) :: km
     END SUBROUTINE sort_unique
   END INTERFACE
! -----------------------------------------------------------------------
! first we need to know if the beam namelist file and the ufile exist:
! they may be located in directories other than the current working directory
! -----------------------------------------------------------------------
    INQUIRE (file = beam_data_namelist, iostat = iostat,               &
                   exist = nexists, opened = nopened,number=nlun)
    INQUIRE (file = beam_data_ufile, iostat = iostat,               &
                   exist = uexists, opened = uopened,number=ulun)

    IF( use_nubeam)THEN

       IF(.NOT. nexists)THEN
          PRINT *,'Error, nubeam selected but input file'
          PRINT *,beam_data_namelist(1:LEN_TRIM(beam_data_namelist))
          PRINT *,'does not exist'
          CALL STOP(' get_beam_data',1)
       ENDIF

       IF((use_ufile).AND.(.NOT. uexists))THEN 
          PRINT *,'Error, ufile selected but input file'
          PRINT *,beam_data_ufile(1:LEN_TRIM(beam_data_ufile))
          PRINT *,'does not exist'
          CALL STOP(' get_beam_data',1)
       ENDIF

       IF ((nprim.EQ.1) .AND. (nubeam_version .LT. 201107) )THEN
          PRINT *,'Error, nubeam selected with nprim = 1'
          CALL STOP(' get_beam_data',1)
       ENDIF

    ELSE  !dont use nubeam but use ufile for input if present

      IF(.NOT. nexists)THEN 
         ! Dont use nubeam and dont use nubeam namelist because it dosen't exist
         ! inone was read and inputs are allready set from inone
         ! hence there is nothing to be done
          RETURN            
      ENDIF        

    ENDIF



! ----------------------------------------------------------------------
!  beam_data_namelist exists and will be read
! ----------------------------------------------------------------------
    read_com = .FALSE.                      !dont read comment field in ufile
  naml:  IF(nexists)THEN
       nubeam_namelist_cpy(:) = ' '
       nubeam_namelist_cpy(1:LEN_TRIM(beam_data_namelist)) = beam_data_namelist(1:LEN_TRIM(beam_data_namelist))
       !get available io unit no.,nlun :
       nlun = 64                               !may be changed by getioun
       CALL getioun(nlun,nlun)                 ! dummy doesnt change  nlun
       !open the file with the name obtained from inone input:
       OPEN(unit=nlun,file = beam_data_namelist,status ='OLD',IOSTAT = iostat)
 


       IF(iostat .EQ. 0)THEN

          !strip comments from namelist file before processing:
          !----------------------------------------------------------
          longfile(:) = ' '
          nml_size = LEN_TRIM(beam_data_namelist)
          nml_size = nml_size + 3  !for "_12", see below
!          beam_data_namelist_copy(1:nml_size-3) = beam_data_namelist(1:nml_size-3)
          longfile(1:nml_size-3) = beam_data_namelist(1:nml_size-3)

          limit = LEN(beam_data_namelist_copy)
          CALL strip_path(longpath,prefix,extension,longfile,1,ishot,limit,ncrt)!JMP
!         call to_upper_case1(beam_data_namelist_copy)
          call to_upper_case1(longfile)
          longpath = ADJUSTL(longpath)
          longfile = ADJUSTL(longfile)
          beam_data_namelist_copy = longpath(1:LEN_TRIM(longpath))//   &
          longfile(1:LEN_TRIM(longfile))

          if(beam_data_namelist_copy(1:nml_size-3) &
              .eq. beam_data_namelist(1:nml_size-3))then
              beam_data_namelist_copy(1:nml_size) = beam_data_namelist_copy(1:nml_size-3)//'_12'
              
          endif

            ntemp = 61
            call getioun(ntemp,ntemp)
            open   (unit = ntemp, file = beam_data_namelist_copy , status = 'UNKNOWN')
            !copy from unit nlun to unit ntemp, removing comments:
            call strip_comments (ncrt, nlun, ntemp)
            !delete association with original file (connected to nlun)
            call giveupus(nlun) ; close  (unit = nlun)
            !close the stripped file and reopen it under the old name:
            call giveupus(ntemp) ; close  (unit = ntemp)
            !below we use only the stripped file,by changing the name here:

            beam_data_namelist = beam_data_namelist_copy 

            !open the strippped file
            CALL getioun(nlun,nlun)
            OPEN(unit=nlun,file = beam_data_namelist,status ='OLD',IOSTAT = iostat)
           ! get size of some arrays by reading beam_data_namelist
             CALL nblist_array_size(nlun,nfast,nbeam,ndifbe,ierr)
             REWIND(nlun)     !nlun to be read again below,after array allocation
             IF(ierr .GT. 0)THEN
                PRINT *,'Error in reading file :',beam_data_namelist
                PRINT *,'values found :'
                PRINT *,'ngmax =',ngmax
                PRINT *,'nrhix =',nrhix
                PRINT *,'nbeam =',nbeam
                PRINT *,'ndifbe =',ndifbe
                CALL STOP(' get_beam_data',2)
             ENDIF

 
           !allocate data structure:
            ngmax = nprim + 1 !allow for possible 'dt' species
                              !reset ngmax in sub set_nubeam_thermal
            nrhix = nimp + 1  !reset in sub set_nubeam_impure
            CALL beam_structure_allocate(name_beam,nbeam,nrhix,ngmax,nfast)

            CALL fi_init(d_fast_ion,ndifbe)
           !set defaults for DIII-D beams:
            CALL set_beam_defaults(name_beam,ngmax,nrhix)
 
           !set defaults for other parameters required by nubeam:
            CALL set_nubeam_misc_defaults

           !set thermal species (nubeam has to know about these)
            CALL set_nubeam_thermal

           !set impurity species:
            CALL set_nubeam_impure
 
           !finally load name_beam data structure
            IF (nubeam_version .lt. 201107) THEN 
                CALL nblist_read(nlun,name_beam,ierr)
            ELSE 
                CALL read_nubeam_namelist(nlun,name_beam,ierr)
            ENDIF

            IF(ierr .GT. 0)THEN
               write(nitre,7)
 7              Format(2x,' ********************************',/, &
               '  namelist read failed,try parsing namelist'/,   &
                 ' ****************************************' )
               !namelist read from file beam_data_namelist
               !failed. This is most likely due to the fact that
               !beam_data_namelist doesnt contain a valid  namelist
               ! but could also mean that the variables in the
               ! the namelist are bad. We cant protect agains the latter
               !case at all but the former case may be treatable by
               !parsing the file manually. The parser is not very robust
               !however so this may or may not work:
               REWIND(nlun)
               CALL nblist_read_parse(nlun,name_beam,ierr)
               if(nubeam0_dt .eq. 0.0)nubeam0_dt = nubeam_dt
               if(ierr .gt. 0)then
                  print *,'tried to read file ',beam_data_namelist
                  print *,' twice, once as a namelist and '
                  print *,' the second time as a manually parsed'
                  print *,' file. Neither read was successful '
                  call stop('namelist read for nubeam failed',1)
               endif
             ENDIF


 
           IF(nbeam  .EQ. 0 .OR. ngmax .EQ. 0)THEN
                !nubeam does not allow nbeam =0
                PRINT *,'zero beam case not allowed'
                PRINT *,'ngmax =',ngmax
                PRINT *,'nbeam =',nbeam
                CALL STOP(' get_beam_data',3)
           ENDIF
 

           !if name_beam%tbona and name_beam%tboffa are set only for
           !one beam but more than one beam is input then fix it here:
            DO j=2,name_beam%nbeam
               IF(name_beam%tbona(j) .LT. -1.e6*time0 )          &
                  name_beam%tbona(j) = name_beam%tbona(j-1)
               IF(name_beam%tboffa(j) .LT. -1.e6*time0 )         &
                  name_beam%tboffa(j) = name_beam%tboffa(j-1)
            ENDDO
! save these times for use in P_Nfreya
                 name_beam%tbonac(:)  = name_beam%tbona(:)
                 name_beam%tboffac(:) = name_beam%tboffa(:)

           IF(use_ufile) THEN !JMP BLOCK START
              DO j=1,name_beam%nbeam
              	name_beam%tbona(j) = time0
              END DO
           END IF !JMP BLOCK END

           IF(ierr .NE. 0)read_com = .TRUE.
       ELSE
            read_com = .TRUE.
       ENDIF
 
       CALL giveupus(nlun)         !done with this file,give up the unit
       CLOSE(nlun)
 

       CALL ufile_bypass(read_ufile)     !do we need to read the ufile ?
!       PRINT *,'done calling ufile bypass,read_ufile =',read_ufile,use_ufile
       IF(use_ufile) go to 300 !JMP
       
       IF(.NOT. read_ufile) go to 200    !if not finish up by 
                                         !going  to nbdrive_nubeam_checknaml
       IF(.NOT. use_ufile) go to 200 !JMP

    ELSE naml                            !beam_data_namelist doesnt exist
       !file doesnt exist, not fatal, we will get the info
       !from ufile instead (see below)
       read_com = .TRUE.                     !read ufile comment fields
    ENDIF naml


    print *,'name_beam%nlbdat=',name_beam%nlbdat
    

! ------------------------------------------------------------------------
! now process the ufile, read comments in ufile if read_com = T
! ------------------------------------------------------------------------
300  CONTINUE
     IF(.NOT. name_beam%data_allocated)  THEN !JMP
       IF(nbeam .NE. 0)THEN
           nfast = nbeam+4
           CALL beam_structure_allocate(name_beam,nbeam,nrhix,ngmax,nfast)
       ENDIF
    ENDIF

    !HERE WE ASSUME, FOR NOW ANYWAY, THAT NDIFBE AND THE ASSOCAITED
    !INITIALIZATION OF TYPE fidif (MODULE FAST_ION_DIFFUSION)
    ! WAS SET ABOVE EVEN IF UFILE READS ARE USED. THIS IS BECAUSE
    ! I THINK UFILES DO NOT SUPPLY A VALUE FOR NIDFBE (BUT I MAY 
    ! BE WRONG - IF SO WE NEED TO DO A CALL TO FI_INIT HERE)
    

    NXMAX =100; NYMAX = 32
    ALLOCATE( X(NXMAX)) ;ALLOCATE(Y(NYMAX))
    ALLOCATE(F(NXMAX,NYMAX) )
    limit = LEN(path)
    iflag =0                               !use 2d form for F
    itry =0

      
    IF (uexists .AND. uopened) THEN   ! file is open
        CALL giveupus(ulun)           ! give up the unit
        CLOSE(ulun)
    ENDIF
    IF(.NOT. uexists )THEN   !error         
         WRITE(ncrt,FMT='(2x,"Error : beam_data_ufile not found :",a)')beam_data_ufile
              CALL STOP ('subroutine get_beam_data: bad INQUIRE', 1)

    ENDIF


    !beam_data_ufile exists. open and process it.
    !get the path if it is part of the file specification:
     CALL strip_path(path,prefix,extension,beam_data_ufile,0,ishot,limit,ncrt) !JMP
          

    !get available io unit no.,ulun :
    ulun = 64                               !may be changed by getioun
    CALL getioun(ulun,ulun)
!   print *,'beam_data_ufile =',beam_data_ufile
!   print *,'prefix =',prefix
!   !print *,'suffix =',suffix
!   !print *,'disk =',disk
!   !print *,'directory =',directory
!   print *,'extension =',extension
!   print *,'path =',path 


    !call setup routine required for ufiles read:
    !UFSETR is used to set up the naming and location convention
    !you will use.  The disk and directory fields may be blank, in 
    !which case the files will default to the current working directory.  
    !call UFSETR(ulun,prefix,suffix,disk,directory)
    !          ulun - integer
    !          prefix - character, max length 16
    !          suffix - character, max length 16
    !          disk - character, max length 16
    !          directory - character, max length 64
     CALL UFSETR(ulun,prefix,extension,' ',path) !JMP
 
     !set the shot number :
 100 CALL UFOPRD(ulun,ishot,ier)     !Open for READONLY access
     IF(ier .NE. 0)THEN
        WRITE(ncrt,1)ier
 1      FORMAT(2x,'get_beam_data returned ier =',i5,/, &
               2x,'while calling UFOPRD')
        CALL STOP('get_beam_data',1)
     ENDIF
 
 
 
 
 
     !now do the actual reading:
     CALL  UF2DRD(ulun, TDEV, SHDATE, F, NXMAX, X, NXMAX, Y, NYMAX, IPROC, &
     NX, NY, IFLAG, FLAB,XLAB, YLAB, NSCMAX, NSC,SCVAL,SCLAB,IER)
 
 
 
     IF(ier .NE. 0)THEN
        itry = itry + 1
        WRITE(ncrt,2)ier
 2      FORMAT(2x,'get_beam_data returned ier =',i5,/, &
               2x,'while calling UF2DRD')
        WRITE(ncrt,3)nx,nxmax,ny,nymax
 3      FORMAT(2x,'nx, nxmax =',2(2x,i5),/, &
               2x,'ny, nymax =',2(2x,i5))
 
        IF(itry .GT. 1) THEN
            CALL STOP('get_beam_data',2)
        ELSE
           WRITE(ncrt,4)
 4         FORMAT(2x,'Reallocating arrays,calling UF2DRD again')
           IF(nx .NE. nxmax)THEN
                  DEALLOCATE(X)
                  ALLOCATE(X(nx))
                  nxmax = nx
            ENDIF
           IF(ny .NE. nymax)THEN
                  DEALLOCATE(y)
                  ALLOCATE(y(ny))
                  nymax = ny
           ENDIF
           DEALLOCATE(F)
           ALLOCATE(F(nx,ny))
           CALL giveupus(ulun)              !give up the unit
           CALL UFCLOS(ulun)                !Close the file
           go to 100
        ENDIF
     ENDIF
 
     CALL giveupus(ulun)              !give up the unit
     CALL UFCLOS(ulun)                !Close the file


    IF(NSC .NE. 1)THEN
       PRINT *,'UFILE read reports number of scalars =',nsc
       PRINT *,' Expected to read NSC =1 scalar'
       CALL STOP('get_beam_data',5)
    ENDIF
  
    IF(name_beam%nbeam .NE. INT(scval(nsc)))THEN
       PRINT *,'UFILE read reports number of beams =',scval(nsc)
       PRINT *,'namelist read has nbeams =',name_beam%nbeam
       CALL STOP('get_beam_data',6)
    ENDIF

    !process data for use in Onetwo:
     CALL set_beam_data(F,x,y,nx,ny)
     CALL  nbdrive_nubeam_checknaml(ierr)



 200  IF(ALLOCATED(x)) DEALLOCATE(x) ;IF(ALLOCATED(y))DEALLOCATE(y) 
     IF(ALLOCATED(f)) DEALLOCATE(f)


    RETURN
    END SUBROUTINE get_beam_data



   SUBROUTINE  set_beam_data(F,x,y,nx,ny)
!  ------------------------------------------------------------------
!  The 2d UFILE  contains nx time points in x,
!  the channel no in y NOTE: r4 not r8
!  ---------------------------------------------------HSJ-11/18/03---
   USE transp,                       ONLY : beam_data,use_nubeam
   USE P_Nfreya_12_interface,        ONLY : use_P_Nfreya
   IMPLICIT NONE
   REAL*4 , INTENT (IN),DIMENSION(:,:) ::     F
   REAL*4 , INTENT (IN),DIMENSION(:) :: x,y
   INTEGER, INTENT (IN) :: nx,ny
   INTEGER ufx,ufy,lfx,lfy,uxx,lxx,uyx,lyx,k,j,l,m

   ufx = UBOUND(F,1)
   ufy = UBOUND(F,2) 
   lfx = LBOUND(F,1)
   lfy = LBOUND(F,2)
   uxx = UBOUND(x,1)
   lxx = LBOUND(x,1)
   uyx = UBOUND(y,1)
   lyx = LBOUND(y,1)

   IF(ufx .NE. uxx .OR. ufy .NE. uyx .OR. lfx .NE. lxx .OR. lfy .NE. lyx) THEN 
      PRINT *,'UFILE error in dimensions'
      PRINT *,'ufx,ufy,lfx,lfy =',ufx,ufy,lfx,lfy
      PRINT *,'uxx,lxx,uyx,lyx =',  uxx,lxx,uyx,lyx
      CALL STOP('set_beam_data',1)
   ENDIF

   IF( ASSOCIATED(beam_data%beam_times))  &
      DEALLOCATE(beam_data%beam_times)
      ALLOCATE(beam_data%beam_times(nx))
   IF( ASSOCIATED(beam_data%beam_inject))  &
      DEALLOCATE(beam_data%beam_inject)
      ALLOCATE(beam_data%beam_inject(nx,ny))
   IF( ASSOCIATED( beam_data%beam_chan))  &
      DEALLOCATE(beam_data%beam_chan)
      ALLOCATE(beam_data%beam_chan(ny))

   beam_data%beam_inject(:,:) = f(:,:)
   beam_data%beam_times(:) = x(:)
   beam_data%beam_chan(:)  = y(:)
   IF( .NOT. use_nubeam .AND. .NOT. use_P_Nfreya) CALL nfreya_load 

   RETURN
   END SUBROUTINE  set_beam_data
   


 SUBROUTINE nfreya_load
! ------------------------------------------------------------------------
! setup nfreya type input for beams
! this data comes from the ufiles and the nubeam namelist input
! NOT from inone !!


 
!  try to construct some beams suitable for use in Freya
!  by analyzing the beam input from ufiles or nubeam namelist.
!  Beams have many attributes. The primary groups are
!  injection geometry,species,energy and power.
!  One set of choices from each group defines a single beam.
!  We need to map these beams into at most kb "pseudo beams"
!  for Freya. Freya  has pseudo beams in the sense that 
!  each beam consists of two sources. Other than that the attributes of
!  each individual beam have to be the same . 
      
!  info for fast species(injected and fusion products) is here:
!  beam_data%label_thions(:)
!  beam_data%label_fastions(:)
!  beam_data%label_impions(:)
!  beam_data%nfusion_species 
!  beam_data%nbeam_species
!  number of fast species = beam_data%nbeam_species 
!                                   + beam_data%nfusion_specie
!  info for beam geometry:
!  In DIII-D the left source of a beamline is the "tangential"
!  source. It has a tangency radius of 114.6cm
!  The right source of a beamline is the "perpendicular" source
!  it has atangency radius  of 74.94 cm
!  Sfrac1(ib) means the perpendicular source of beam line ib (ib =1,2..kb)
!  The tangential source supplies  1.-sfrac1(ib) fraction of the ions.
!  For Freya we have up to kb beams each with two sources
!  but only 1 beam species is allowed. ( DT mixture is
!  treated as one effective species in Freya.)
!  the sources determine the tangency radius.
! --------------------------------------------------HSJ-12/03/03----------


   USE param,                              ONLY : kb
   USE transp,                             ONLY : beam_data

   USE ions,                               ONLY : nameb
   USE fusion,                             ONLY : fdbeam,fd
   USE nub,                                ONLY : fbcur,ebkev,bptor, &
                                                  sfrac1,nbeams,beamon,btime
   USE nub3,                               ONLY : iexcit,izstrp
   USE tordlrot,                           ONLY : angrcple,nbeamtcx

   USE P_Nfreya_12_interface,              ONLY : use_P_Nfreya

   IMPLICIT NONE  
   INTEGER j,l,k,ks,ufx,ufy,lfx,lfy,uxx,lxx,uyx,             &
           lyx,match,nfbeam,kbspec,ksort

   LOGICAL  set_diiid,first

   REAL *8 rtangcy(beam_data%nbeam)            !temporary array
   REAL *8 sfrac2(beam_data%nbeam)          !temporary array
   INTEGER*4 nubeam_group(beam_data%nbeam)  !temporary array

  INTERFACE local
     SUBROUTINE sort_unique(a,km)
       REAL *8, INTENT(INOUT),DIMENSION(:) :: a
       INTEGER, INTENT(OUT) :: km
     END SUBROUTINE sort_unique
     SUBROUTINE set_freya_beam(nfb,nnub,first,set_diiid,match,sfrac2)
       INTEGER,INTENT(in) :: nfb,nnub,match
       REAL*8,INTENT(OUT):: sfrac2(:)
       LOGICAL set_diiid,first
     END SUBROUTINE set_freya_beam
   END INTERFACE

   sfrac1(:) =0.0
   sfrac2(:) =0.0



      !decide on uniqe beam species
       IF(beam_data%nbeam_species .GT. 2)THEN
          PRINT *,'More than 2 injected beam species not allowed' 
          CALL STOP('too many beam species',1)
       ELSE IF(beam_data%nbeam_species .GT. 1)THEN
          !max of 2 allowed, must be 'd' , 't'
          kbspec=2
          DO j=1,beam_data%nbeam_species
             match=0
             DO l=1,kbspec
!                IF(beam_data%label_injected_ions(j) == 'd'      &
!                 .or.  beam_data%label_injected_ions(j) == 't')THEN
!                   match=1
!                   EXIT
!                 ENDIF
              ENDDO
              IF(match == 0) THEN
                   PRINT *,'If two beam species are given'
                   PRINT *,'they must be "d" and "t" '
!                   print *,'input is ',beam_data%label_injected_ions
                   CALL STOP('Invalid injected species',1)
              ENDIF
          ENDDO
       ELSE !only single species beam
          kbspec =1 
       ENDIF


      !decide on sources based on tangency radii:
      rtangcy(1:beam_data%nbeam) = beam_data%rtcena(1:beam_data%nbeam)

       CALL sort_unique(rtangcy,ks)               !There are ks unique sources
 
      IF(use_P_Nfreya)ks =2 ! value not used in P_Nfreya, only nfreya is limited to 2
      IF(ks .GT. 2)THEN
         PRINT *,'Not set up to handle more than 2 tangency radii'
         PRINT *,'Number of sources specified =',ks
         PRINT *,'with tangency radii :', beam_data%rtcena(:)
         CALL STOP('nfreya_load ',1)
       ENDIF
 



       !check the times we allow only one beam on time and one beam off time
       !for nfreya. OD Fokker Planck time dependentbeam model in Onetwo
       ! is not accounted for here at this time.

       beamon = beam_data%tbona(1)  !if there is only one beam on time
                                    !then the first element of beam_data%tbona
                                    !will have it

       btime = beam_data%tboffa(1) - beam_data%tbona(1)



       !OK, we have ks(<=2) unique sources and kbspec(<=2)
       ! unique beam species, and kbe unique injection
       ! full energies and current fractions
       !and only one beam on time and one beam off time.
       !info is sufficient to construct Freya input:




      IF(beam_data%nbeam_species  .GT. 1)THEN
           fdbeam = fd !fd is # fraction of thermal dueterons
           !in thermal dt mixture ??
         !more than one beam species. Onewo can only handle
         !this case if the beam is 'd' and 't'
         !In this case we will create a " single species"
         !dt beam.
        CALL STOP('dt from nubeam for Freya not implemented',1)
      ELSE !one beam species.
           !get number of distinct sources 
           !(eg perpendicular and/or parallel sources)
         nbeams =0
         DO j=1,beam_data%nbeam

         IF(beam_data%abeams(1) .LT. 2.)THEN
            nameb='h'
            fdbeam        = 0.150e-3    ! isotopic content of d in h  
         ELSE IF(beam_data%abeams(1) .LT. 3.)THEN
            nameb='d'
            fdbeam = 1.0
         ELSE IF(beam_data%abeams(1) .LT. 4.)THEN
            nameb='t'
            fdbeam =0.0
         ELSE 
            PRINT *,'beam species not recognized'
            PRINT *,'A,Z =',beam_data%abeams(1),beam_data%xzbeams(1)
            CALL STOP('unknown beam species',1)
         ENDIF
         ENDDO
      ENDIF


    
      GO TO 10001


      !OBSOLETE CODING:
      nfbeam =1
      nubeam_group(:) = -1 
      nubeam_group(1) = 1 !nubeam_group(i) = j indicates that nubeam beam i is
                          !associated with freya beam j
      first =.TRUE.
      set_diiid = .FALSE.
      match = 0
      CALL set_freya_beam(nfbeam,1,first,set_diiid,match,sfrac2) ! set first freya beam to first nubeam beam
      first =.FALSE.
      DO j=2,beam_data%nbeam      !   loop over no. of beams given in nubeam namelist

         DO k = 1, j-1            !   check previous assigned nfreya  beams for match
            CALL check_beam_attrib(k,j,match)
 
            IF(match .NE. 0 )THEN !nubeam beam j matches previously 
                                  !registerd  nubeam beam k exactly
                                  !or only differs in specification of
                                  !tangency radius.
                                  !nubeam beam k is associated with 
                                  !freya beam group(k)
               CALL set_freya_beam(nfbeam,nubeam_group(k),first,set_diiid, &
                                   match,sfrac2)
            ELSE   ! nubeam beam j does not match nubeam beam k 
               nfbeam = nfbeam+1      !requires new freya beam
               IF(nfbeam .LE. kb)THEN
                 nubeam_group(j) = nfbeam
                 CALL set_freya_beam(nfbeam,j,first,set_diiid,match,sfrac2)
               ELSE
                 PRINT *,'not enough beams available for FREYA'
                 PRINT *,'Manual editing of beams for FREYA required'
                 CALL STOP('nfbeam too large',1)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!---------------------new expanded eams-------------------------------------
10001 CONTINUE
      NFBEAM =0
      first = .TRUE. ; match  = 0 ;       set_diiid = .FALSE.
      DO j=1,beam_data%nbeam
               nfbeam = j       ! use 1:1 mapping requires new freya beam
               IF(nfbeam .LE. kb)THEN
                 nubeam_group(j) = nfbeam
                 CALL set_freya_beam(nfbeam,j,first,set_diiid,match,sfrac2)
               ELSE
                 PRINT *,'not enough beams available for FREYA'
                 PRINT *,'Manual editing of beams for FREYA required'
                 CALL STOP('nfbeam too large',1)
               ENDIF
               first = .FALSE.
      ENDDO
!--------------------------------------------------------------------------
      nbeams = nfbeam        !nbeams is final # of Freya beams

!change sfrac to actual fractions and change bptor to total power:
      
       DO j=1,nbeams
          bptor(j) = sfrac1(j) + sfrac2(j)
          sfrac1(j) = sfrac1(j)/bptor(j)
       ENDDO

       PRINT *,'ebkev =', ebkev
       PRINT *,'bptor =', bptor 
       PRINT *,'nbeams =',nbeams
       print *,'nubeam_group(1:beam_data%nbeam) =',nubeam_group(1:beam_data%nbeam)

       DO j = 1,nbeams
         PRINT *,'fbcur =', fbcur(1:3,j)
         PRINT *,'sfrac1 =',sfrac1(j)
       ENDDO




   !use nfreya etal with  data from ufile:
   IF(ASSOCIATED(beam_data%beam_inject))THEN
      ufx = UBOUND(beam_data%beam_inject,1)
      ufy = UBOUND(beam_data%beam_inject,2) 
      lfx = LBOUND(beam_data%beam_inject,1)
      lfy = LBOUND(beam_data%beam_inject,2)
      PRINT *,'ufx,ufy,lfx,lfy =',ufx,ufy,lfx,lfy
   ENDIF
   IF(ASSOCIATED(beam_data%beam_times))THEN
      uxx = UBOUND(beam_data%beam_times,1)
      lxx = LBOUND(beam_data%beam_times,1)
      PRINT *,'uxx,lxx  =',uxx,lxx
   ENDIF
   IF(ASSOCIATED(beam_data%beam_chan))THEN
      uyx = UBOUND(beam_data%beam_chan,1)
      lyx = LBOUND(beam_data%beam_chan,1)
      PRINT *,'uyx,lyx =',uyx,lyx
   ENDIF




      !deallocate the arrays, they are not used by  nfreya and friends
      IF(ASSOCIATED(beam_data%beam_chan)) DEALLOCATE(beam_data%beam_chan)
      IF(ASSOCIATED(beam_data%beam_inject))DEALLOCATE(beam_data%beam_inject)
      IF(ASSOCIATED(beam_data%beam_times))DEALLOCATE(beam_data%beam_times)





       
      RETURN

 END SUBROUTINE nfreya_load






   SUBROUTINE nblist_read(lun,nbname,ierr)
!  -------------------------------------------------------------------
!  read beam geometry data from beam_data_namelist
!  the data consist of up to 32 sources (each source is a beamline
!  in nubeam) The ordering is not unique apparently. Instead each
!  input "block"  completely specfies the beamline. beamlines that are
!  off are not output in the namelist file and the on  beams
!  are numbered contiguously from 1,2,3..  
!   we can identify
!  left sources by (tangency radii) rtcena(j) = 114.6cm and 
!  right sources by rtcena(j) = 76.2
!  There appears to be no information that allows us to distinguish
!  between beamlines however
!  NBLIST.FOR !   Note: 8 sources in order( HSJ what order?):
!  NBLIST.FOR !   30LEFT, 30RIGHT, 150LEFT, 150RIGHT,
!  NBLIST.FOR !   210LEFT, 210RIGHT, 330LEFT, 330RIGHT

!  ----------------------------------------------------HSJ-11/18/03----
   USE beam_structure, ONLY : neutral_beam
   USE nbnamelist
   USE transp,ONLY : nptcls,nptclf,ndep0,nubeam_nclass,nubeam_dt,   &
                     nubeam0_dt,d_fast_ion,ndifbep
   IMPLICIT NONE
   TYPE(neutral_beam),INTENT(out)  :: nbname
   INTEGER , INTENT (in ) :: lun
   INTEGER , INTENT (out ) :: ierr
   INTEGER ll,j,k,kk
   INTEGER nimp,inta,intmax
 
     ierr = 0
     nbeam = 0
     nimp = 0
     intmax = 8  


     READ(lun,nml = nbdrive_naml,END =10,ERR  =10)

!    READ(lun,nml= nbdrive_naml)
!    set those quantities read from the namelist that are used
!    by Onetwo and are stored in the nbname data structure:
     nbname%abeama(:)  = abeama(1:nbeam)
     nbname%xzbeama(:) = xzbeama(1:nbeam)
     nbname%bmwidra(:) = bmwidra(1:nbeam)
     nbname%bmwidza(:) = bmwidza(1:nbeam) 
     nbname%rtcena(:)  = rtcena(1:nbeam)
     nbname%xlbtna(:)  = xlbtna(1:nbeam)
     nbname%xybsca(:)  = xybsca(1:nbeam)
     nbname%divza(:)   = divza(1:nbeam)
     nbname%divra(:)   = divra(1:nbeam) 
     nbname%foclza(:)  = foclza(1:nbeam)
     nbname%foclra(:)  = foclra(1:nbeam)
     nbname%rapedga(:) = rapedga(1:nbeam)
     nbname%xzpedga(:) = xzpedga(1:nbeam)
     nbname%xlbapa2(:) = xlbapa2(1:nbeam)
     nbname%xlbapa(:)  = xlbapa(1:nbeam)
     nbname%xybapa(:)  = xybapa(1:nbeam)
     nbname%rapedg2(:) = rapedg2(1:nbeam)
     nbname%xzpedg2(:) = xzpedg2(1:nbeam)
     nbname%xbzeta(:)  = xbzeta(1:nbeam) 
     nbname%pinja(:)   = pinja(1:nbeam)
     nbname%einja(:)   = einja(1:nbeam)
     nbname%ffulla(:)  = ffulla(1:nbeam)
     nbname%fhalfa(:)  = fhalfa(1:nbeam) 
     nbname%ntrace(:)  = ntrace(1:nbeam)
     nbname%nbshapa(:) = nbshapa(1:nbeam)
     nbname%nbapsha(:) = nbapsha(1:nbeam)
     nbname%nlco(:)    = nlco(1:nbeam)
     nbname%nbeam      = nbeam
     nbname%nbbcal     = nbbcal
     nbname%nseed      = nseed 
     nubeam_nclass     = nclass
     nbname%nzone_nb   = nzones
     nbname%tbona(:)   = tbona(1:nbeam)
     nbname%tboffa(:)  = tboffa(1:nbeam)
     nubeam0_dt        = nubeam_dt
     d_fast_ion%adiff_a = adiff_a
     d_fast_ion%adiff_0 = adiff_0
     d_fast_ion%adiff_xpin = adiff_xpin
     d_fast_ion%adiff_xpout = adiff_xpout
     d_fast_ion%nkdifb      = nkdifb
!    d_fast_ion%ndifbe was read by subroutine nblist_array_size
     IF(d_fast_ion%ndifbe .gt. 0)then
         IF(ndifbep .lt. d_fast_ion%ndifbe)THEN
            PRINT *,'ERROR, ndifbep = ',ndifbep
            PRINT *,'But requested array is size',d_fast_ion%ndifbe
            CALL STOP(' nblist_read ',1)
         ENDIF
         d_fast_ion%fdifbe(1:ndifbe) = fdifbe(1:ndifbe)
         d_fast_ion%edifbe(1:ndifbe) = edifbe(1:ndifbe)
     ENDIF
     

!     nptcls  these are also read from the namelist 
!             but they are stored in tranp.mod  rather than in nbname
!     nptclf
!     ndep0

!     finally there are some quantities that are accepted by the
!     namelist read for purposes of compatibility
!     only. These quantities are not sued in Onetwo at all.
!     Examples are frac_ions,frac)imp, etc.


     IF(nbeam .GT. nbeamx)THEN
        CALL STOP('nbeam >   nbeamx,tbona,tboffa too small ',1)
     ENDIF

   RETURN


10   ierr =1

   RETURN

   END





   SUBROUTINE nblist_read_parse(lun,nbname,ierr)
!  -------------------------------------------------------------------
!  read beam geometry data from beam_data_namelist
!  the data consist of up to 8 sources (each source is a beamline
!  in nubeam) The ordering is not unique apparently. Instead each
!  input "block"  completely specfies the beamline. beamlines that are
!  off are not output in the namelist file and the ona beams
!  are numbered contiguously from 1,2,3..  
!   we can identify
!  left sources by (tangency radii) rtcena(j) = 114.6cm and 
!  right sources by rtcena(j) = 76.2
!  There appears to be no information that allows us to distinguish
!  between beamlines however
!  NBLIST.FOR !   Note: 8 sources in order( HSJ what order?):
!  NBLIST.FOR !   30LEFT, 30RIGHT, 150LEFT, 150RIGHT,
!  NBLIST.FOR !   210LEFT, 210RIGHT, 330LEFT, 330RIGHT

!  ----------------------------------------------------HSJ-11/18/03----
   USE beam_structure, ONLY : neutral_beam
   USE string_util,ONLY : strip_path,to_upper_case1,to_lower_case1,       &
                          set_r8_1darray_value,set_i4_1darray_value,      &
                          set_L_1darray_value,set_r4_value,               &
                          set_integer_value, set_r4_1darray_value,        &
                          set_r8_value
   USE transp,ONLY : nptcls,nptclf,ndep0,nubeam_nclass,nubeam0_dt,        &
                     d_fast_ion,ndifbep,fdifbe,edifbe
   IMPLICIT NONE
   TYPE(neutral_beam),INTENT(out)  :: nbname
   INTEGER , INTENT (in ) :: lun
   INTEGER , INTENT (out ) :: ierr
   INTEGER ll,j,k,kk,ndifbe
   INTEGER nbeam,nimp,inta,intmax
   CHARACTER(len=256) :: line,line_prev
   LOGICAL no_cont,last_line

   

   !the file actually does not define a formal namelist so we have to 
   !punt to fetch out the data:
   nbeam = 0   ; ndifbe =0
   nimp = 0
   intmax = 8  
   ierr = 0
   line_prev(1:1) = '!'
   last_line = .FALSE.

   DO WHILE(1 .gt.  0)
     no_cont = .TRUE.
     line = line_prev

     IF( .NOT. last_line)  THEN                            
        READ(lun,FMT='(A)',ERR = 15,END=10)line_prev
        go to 11
 10     last_line = .TRUE.
        !check if line just read is a continuation of previous line:
 11     line_prev  = ADJUSTL(line_prev)
        !line_prev = to_upper_case(line_prev)
        CALL to_upper_case1(line_prev)
        ll = LEN_TRIM(line_prev) 
        DO j=1,ll
           IF(line_prev(j:j) == ' ')CYCLE
           IF(line_prev(j:j) == '!')EXIT
           IF(line_prev(j:j) == ';')EXIT
           !first character on line must be '=','+','-' or digit
           !otherwise this is not a continuation line
           IF(line_prev(j:j) == '=' .OR.            &
              line_prev(j:j) == '+' .OR.            &
              line_prev(j:j) == '-' .OR.            &
              line_prev(j:j) == '0' .OR.            &
              line_prev(j:j) == '1' .OR.            &
              line_prev(j:j) == '2' .OR.            &
              line_prev(j:j) == '3' .OR.            &
              line_prev(j:j) == '4' .OR.            &
              line_prev(j:j) == '5' .OR.            &
              line_prev(j:j) == '6' .OR.            &
              line_prev(j:j) == '7' .OR.            &
              line_prev(j:j) == '8' .OR.            &
              line_prev(j:j) == '9')THEN
                  no_cont = .FALSE.
                  PRINT *,'lie_prev(j:j) =',line_prev(j:j)
                  PRINT *,'NOT SET UP TO READ CONTINUATION LINES'
                  PRINT *,'Line encountered is:'
                  PRINT*,line_prev(1:ll)
                  CALL STOP('nblist_read ',1)
            ELSE
               EXIT
            ENDIF
        ENDDO
     ENDIF


     line  = ADJUSTL(line)
     !line = to_upper_case(line)
     CALL to_upper_case1(line)
     ll = LEN_TRIM(line)                    
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
!     IF(k == 0)k=ll                                 !string does not have '!'
     IF(k ==0 ) THEN
         k =ll
     ELSE
         k = k-1
     ENDIF


     kk= INDEX(line(1:k),'ABEAMA')
     IF(kk .NE. 0)                                            &        
        CALL set_r8_1darray_value(k,kk,line,nbname%abeama)

     kk= INDEX(line(1:k),'XZBEAMA')
     IF(kk .NE. 0)                                             &
        CALL set_r8_1darray_value(k,kk,line,nbname%xzbeama)

!     kk= INDEX(line(1:k),'APLASM')
!     IF(kk .NE. 0)                                            &   
!        CALL set_r8_1darray_value(k,kk,line,nbname%aplasm)

!     kk= INDEX(line(1:k),'BACKZ')
!     IF(kk .NE. 0)                                            &   
!        CALL set_r8_1darray_value(k,kk,line,nbname%backz)


     kk= INDEX(line(1:k),'BMWIDRA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%bmwidra)



     kk= INDEX(line(1:k),'BMWIDZA')
     IF(kk .NE. 0)                                            &
          CALL set_r8_1darray_value(k,kk,line,nbname%bmwidza)



     kk= INDEX(line(1:k),'RTCENA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%rtcena)


     kk= INDEX(line(1:k),'XLBTNA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xlbtna)


     kk= INDEX(line(1:k),'XYBSCA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xybsca)


     kk= INDEX(line(1:k),'DIVZA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%divza)

     kk= INDEX(line(1:k),'DIVRA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%divra)


     kk= INDEX(line(1:k),'FOCLZA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%foclza)


     kk= INDEX(line(1:k),'FOCLRA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%foclra)



     kk= INDEX(line(1:k),'RAPEDGA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%rapedga)

     kk= INDEX(line(1:k),'XZPEDGA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xzpedga)


     kk= INDEX(line(1:k),'XLBAPA2')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xlbapa2)



     kk= INDEX(line(1:k),'XLBAPA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xlbapa)



     kk= INDEX(line(1:k),'XYBAPA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xybapa)



     kk= INDEX(line(1:k),'RAPEDG2')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%rapedg2)



     kk= INDEX(line(1:k),'XZPEDG2')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xzpedg2)


!     kk= INDEX(line(1:k),'XZIMPX')
!     IF(kk .NE. 0)                                            &
!         CALL set_r8_1darray_value(k,kk,line,nbname%xzimpx)



!     kk= INDEX(line(1:k),'AIMPX')
!     IF(kk .NE. 0)                                            &
!         CALL set_r8_1darray_value(k,kk,line,nbname%aimpx)


     kk= INDEX(line(1:k),'XBZETA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xbzeta)

     kk= INDEX(line(1:k),'PINJA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%pinja)


     kk= INDEX(line(1:k),'EINJA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%einja)


     kk= INDEX(line(1:k),'FFULLA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%ffulla)

     kk= INDEX(line(1:k),'FHALFA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%fhalfa)


     kk= INDEX(line(1:k),'NTRACE')
     IF(kk .NE. 0)                                            &
         CALL set_i4_1darray_value(k,kk,line,nbname%ntrace)



     kk= INDEX(line(1:k),'NPTCLS')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nptcls)



     kk= INDEX(line(1:k),'NPTCLF')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nptclf)


     kk= INDEX(line(1:k),'NDEP0')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,ndep0)


     kk= INDEX(line(1:k),'NBSHAPA')
     IF(kk .NE. 0)                                            &
         CALL set_i4_1darray_value(k,kk,line,nbname%nbshapa)

    kk= INDEX(line(1:k),'NBAPSHA')
     IF(kk .NE. 0)                                            &
         CALL set_i4_1darray_value(k,kk,line,nbname%nbapsha)


    kk= INDEX(line(1:k),'NLCO')
     IF(kk .NE. 0)                                            &
         CALL set_L_1darray_value(k,kk,line,nbname%nlco)

    
    kk= INDEX(line(1:k),'NBEAM')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nbname%nbeam)

!    kk= INDEX(line(1:k),'NGMAX')
!     IF(kk .NE. 0)                                            &
!         CALL set_integer_value(k,kk,line,nbname%ngmax)

!    kk= INDEX(line(1:k),'NRHIX')
!     IF(kk .NE. 0)                                            &
!         CALL set_integer_value(k,kk,line,nbname%nrhix)


    kk= INDEX(line(1:k),'NBBCAL')
     IF(kk .NE. 0)                                            & 
         CALL set_integer_value(k,kk,line,nbname%nbbcal)
 
    kk= INDEX(line(1:k),'NSEED')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nbname%nseed)
 
    kk= INDEX(line(1:k),'NCLASS')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nubeam_nclass)


 
    kk= INDEX(line(1:k),'NZONES')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nbname%nzone_nb)

    kk= INDEX(line(1:k),'TBONA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%tbona)
 


    kk= INDEX(line(1:k),'TBOFFA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%tboffa)

 

    kk= INDEX(line(1:k),'NUBEAM_DT')
     IF(kk .NE. 0)                                            &
         CALL set_r8_value(k,kk,line,nubeam0_dt)



    kk= INDEX(line(1:k),'ADIFF_0')
     IF(kk .NE. 0)                                            &
         CALL set_r8_value(k,kk,line,d_fast_ion%adiff_0)



    kk= INDEX(line(1:k),'ADIFF_A')
     IF(kk .NE. 0)                                            &
         CALL set_r8_value(k,kk,line,d_fast_ion%adiff_a)




    kk= INDEX(line(1:k),'ADIFF_XPIN')
     IF(kk .NE. 0)                                            &
         CALL set_r8_value(k,kk,line,d_fast_ion%adiff_xpin)



    kk= INDEX(line(1:k),'ADIFF_XPOUT')
     IF(kk .NE. 0)                                            &
         CALL set_r8_value(k,kk,line,d_fast_ion%adiff_xpout)



    kk= INDEX(line(1:k),'NDIFBE')
     IF(kk .NE. 0) then                                      
         CALL set_integer_value(k,kk,line,ndifbe)
         IF( ndifbe .gt. ndifbep )then                        
              PRINT *,'ERROR, ndifbe must be .le. ',ndifbep
              PRINT *,'value of ndifbe specified =',ndifbe
              CALL STOP('nblist_read_parse',1)
         ENDIF
     ENDIF

    kk= INDEX(line(1:k),'NKDIFB')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,d_fast_ion%nkdifb)

     kk= INDEX(line(1:k),'FDIFBE')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,fdifbe)

     kk= INDEX(line(1:k),'EDIFBE')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,edifbe)

    IF(last_line)EXIT
   END DO
   
   IF(d_fast_ion%ndifbe .gt. 0)then
      d_fast_ion%fdifbe(1:d_fast_ion%ndifbe) = fdifbe(1:d_fast_ion%ndifbe)
      d_fast_ion%edifbe(1:d_fast_ion%ndifbe) = edifbe(1:d_fast_ion%ndifbe)
   ENDIF

   RETURN

15 ierr = 1
   RETURN

   END




       SUBROUTINE nblist_array_size(lun,nfast,nbeam,ndifbe,ierr)
!  -------------------------------------------------------------------
!  get values of ngmax,nrhix,nbeam
!  nbeam = # beams
!  nfast = # fast ions ??
!  ----------------------------------------------------HSJ-11/18/03----
   USE string_util,ONLY : to_upper_case1,set_integer_value
   IMPLICIT NONE
   INTEGER,INTENT(OUT) :: nfast,nbeam,ierr,ndifbe

   INTEGER , INTENT (in ) :: lun
   INTEGER ll,j,k,kk
   INTEGER n,inta
   CHARACTER(len=256) :: line

   nbeam = 0 ;  ierr = 0
   DO WHILE(1 .gt. 0)
     READ(lun,FMT='(A)',ERR = 15,END=10)line
     line  = ADJUSTL(line)
     !line = to_upper_case(line)
     CALL to_upper_case1(line)
     ll = LEN_TRIM(line)                    
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk= INDEX(line(1:k),'NBEAM')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nbeam)



     kk= INDEX(line(1:k),'NDIFBE')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,ndifbe)


!     kk= INDEX(line(1:k),'NGMAX')
!     IF(kk .NE. 0)                                            &
!         CALL set_integer_value(k,kk,line,ngmax)

!     kk= INDEX(line(1:k),'NRHIX')
!     IF(kk .NE. 0)                                            &
!         CALL set_integer_value(k,kk,line,nrhix)
   END DO
 
10 CONTINUE
   !PRINT *,'nbeam,nrhix,ngmax =',nbeam,nrhix,ngmax


   nfast = 7 ! number of injected beam species plus number
             ! of fusion product species

   RETURN

15 ierr = 1
   RETURN

   END




   SUBROUTINE set_beam_defaults(nbname,ngmax,nrhix)
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
   USE param,ONLY:            kb
   USE beam_structure, ONLY : neutral_beam
   USE nub,ONLY:              bptor,ebkev,fbcur,mf
   USE string_util,ONLY:      strip_path  
   USE numbrs,ONLY :          nprim
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: ngmax,nrhix
   TYPE(neutral_beam),INTENT(out)  :: nbname

! -----------------------------------------------------------------------
!      DIII_D defaults  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -----------------------------------------------------------------------

   nbname%ngmax                =  nprim
   nbname%nzone_nb             =  mf
   nbname%nbbcal               =  1000000
   nbname%nseed                =  1234567
   nbname%pinja(:)             =  0.0
   nbname%einja(:)             =  0.0       !beam voltages ev
   nbname%ffulla(:)            =  0.75
   nbname%fhalfa(:)            =  0.15
   nbname%rapedg2(:)=0.0
   nbname%bmwidza(:) =24.0              !Ion source half height in cm
   nbname%xybapa(:) =0.0
   nbname%rapedg2(:)=0.0
   nbname%xzpedg2(:) =0.0     !aperture half-height (ignored if circular), cm

   nbname%nbshapa(:) = 1      !default rectangular ion source gcrid shape
   nbname%nbapsh2(:) =1       !shape of optional 2nd aperture (1:rectangle, 2:circle)

   nbname%rtcena(:) = 114.6  !cm   default tangency radius  vcalues
   nbname%nlco(:) = .TRUE.   !default all co injection
   nbname%ntrace(:) = 0      !no trace elements in beams

   nbname%xlbtna(:) = 802.8
   nbname%xybsca(:) =0.0  
   nbname%bmwidra(:) =6.0             !Ion source half width cm
                                      !value of zero means: no 2nd aperture
   nbname%xlbapa(:) = 186.1           !distance, source grid to apertucre, cm
   nbname%xlbapa2(:)= 0.0     !distance, source grid to aperture, cm
                              !value of zero means: no 2nd aperture
   nbname%divra(:) =0.8730000E-02    !Horizontal divergence in radians
   nbname%xzpedga(:) = 24.0
   nbname%divza(:) =0.2270000E-01    ! Vertical divergence in radians 

   nbname%rapedga(:) = 8.85 ! aperture half-width (or radius if circulccar), cm

   nbname%foclza(:)=1000.000
   nbname%foclra(:) =0.1000000E+34
   nbname%nbapsha(:) = 1
   nbname%abeama(:) = 2.   !D
   nbname%xzbeama(:) = 1.  !for all beams
   

   nbname%tbona(:)  = -1.e30
   nbname%tboffa(:) = -2.* ABS(nbname%tbona(:))  ! tboffa < tbona
   nbname%nznbma  = 50
              ! number of pitch angle zones (evenly spaced in vpll/v) in
              ! lab frame fast ion distribution functions.
              ! note: pitch angle grid the same for all species
   nbname%nznbme = 100
              ! maximum number of energy zones in fast ion distribution
              ! functions; actual number and zone spacing can vary from
              ! fast specie to fast specie.
              ! definition of spatial zones:  see xplasma MC grid.

   nbname%ebdmax= 0.0     !maximum energy in distribution function 
                          !for beam ion species in ev , set in 
                          !nbdrive_nubeam_checknaml


   RETURN

   END    SUBROUTINE set_beam_defaults




   SUBROUTINE set_freya_beam(nfb,nnub,first,set_diiid,match,sfrac2)
! --------------------------------------------------------------------
! set freya beam nfb to attributes of nubeam beam nnub
! if match > 0 then append nnub attributes to nfb.
! if match = 0 then create new nfb entry
! if match = -1 then only differenc is in tangency radius of source
! account for this by adjusting sfrac and appending to nfb
! ---------------------------------------------------------------HSJ-03/24/04
 
   USE param,ONLY : kb
   USE transp ,ONLY : beam_data
   USE nub
   USE nub3, ONLY:  izstrp,iexcit
   USE geom,ONLY : machinei
   USE tordlrot,ONLY : angrcple,nbeamtcx
   IMPLICIT NONE
   REAL*8,INTENT(OUT):: sfrac2(:)
   INTEGER,INTENT(in) :: nfb,nnub,match
   INTEGER i
   LOGICAL set_diiid,first
   

   IF(nnub .LE. 0)THEN
      PRINT *,'ERROR, in sub set_freya_beam,nnub =',nnub
      CALL STOP('nnub problem',1)
   ENDIF
 
   IF(first)THEN
       beamon = beam_data%tbona(nnub)
       nsourc= 2
!      hdepsmth =45.
       if(npart == 0)npart = 50000
       npskip=5 
       iexcit = 5
       izstrp(1)  = 1 ; izstrp(2) =1
       relnub = 0.01
       angrcple=1.00
       nbeamtcx=1
       mf = beam_data%nzone_nb
      iterate_beam=.TRUE.      !dont do thi if beam is initially off
   ENDIF
   
   IF(match == 0)THEN  !new entry
      bptor(nfb)     =  beam_data%pinja(nnub)
      ebkev(nfb)     =  beam_data%einja(nnub)/1000.
      fbcur(1,nfb)   =  beam_data%ffulla(nnub)
      fbcur(2,nfb)   =  beam_data%fhalfa(nnub)
      fbcur(3,nfb)   = 1. -fbcur(1,nfb) -fbcur(2,nfb)
      fbcur(3,nfb)   = MAX(0.0D0,fbcur(3,nfb)) !jmp.ibm
      IF(machinei  .EQ. 'diii-d')THEN

        
         IF(ABS(beam_data%rtcena(nnub)-114.5 ) .LE. 5.)THEN
            sfrac2(nfb)    =  beam_data%pinja(nnub) 
         ELSE IF(ABS(beam_data%rtcena(nnub)-74.9 ) .LE. 5.)THEN
            sfrac1(nfb)    =  beam_data%pinja(nnub) 
         ELSE
            PRINT*, 'error in beam tangency radii for DIII-D'
            CALL STOP('beam spec problem',1)
         ENDIF
      ELSE  ! non DIII-D machines
         sfrac1(1) = 1.0
      ENDIF
   ELSE !match =-1,or 1
      IF(ABS(beam_data%rtcena(nnub)-114.5 ) .LE. 5.)THEN
        sfrac2(nfb)    = sfrac2(nfb) +  beam_data%pinja(nnub) 
      ELSE IF(ABS(beam_data%rtcena(nnub)-74.9 ) .LE. 5.)THEN
            sfrac1(nfb)    = sfrac1(nfb) +  beam_data%pinja(nnub) 
      ELSE
            PRINT*, 'error in beam tangency radii for DIII-D'
            CALL STOP('beam spec problem',2)
      ENDIF

    ENDIF

    !ALL OF THE FOLLOWING GEOMETRY ATTRIBUTES HAVE TO BE SET FROM
    !BEAM_DATA%*** QUANTITES. SKIP FOR NOW
    IF(set_diiid)THEN
!     DIII-D beam input.  DEFAULT IS LONG-PULSE-SOURCE SPECIFICATIONS.         
!                                                                              
      naptr = 4                                                                
!                                                                              
      DO i=1,kb                                                                
        anglev(i)    =   0.0                                                   
        angleh(i)    =  19.5                                                   
        bvofset(i)   =   0.0                                                   
        bhofset(i)   =  42.074                                                 
        bleni(i)     = 556.808                                                 
        bcur(i)      = 110.0                                                   
        bptor(i)     =  0.0e6                                                  
        nbshape(i)   = 'rect-lps'                                              
        bheigh(i)    =  48.0                                                   
        bwidth(i)    =  12.0                                                   
        bhdiv(i)     =   0.50       !degrees                                   
        bvdiv(i)     =   1.3        !degrees                                   
        fbcur(1,i)   =   0.7                                                   
        fbcur(2,i)   =   0.2                                                   
        fbcur(3,i)   =   0.1                                                   
        bhfoc(i)     =   1.0d100                                               
        bvfoc(i)     =   1.0d3                                                 
        ebkev(i)     =  75.0                                                   
        sfrac1(i)    =   0.5                                                   
        nashape(1,i) = 's-rect'                                                
        nashape(2,i) = 's-rect'                                                
        nashape(3,i) = 'b-d3d'                                                  
        nashape(4,i) = 'b-circ'                                                
        aheigh(1,i)  =  47.8                                                   
        awidth(1,i)  =  13.8                                                   
        alen(1,i)    = 186.1                                                   
        aheigh(2,i)  =  48.0                                                   
        awidth(2,i)  =  17.7                                                   
        alen(2,i)    = 346.0                                                   
        alen(3,i)    = 449.0                                                   
        awidth(4,i)  =  50.9                                                   
        alen(4,i)    = 500.0                                                   
        blenp(i)     = 539.0                                                   
        rpivot(i)    = 286.6                                                   
        zpivot(i)    =   0.0                                                   
      END DO  
    ENDIF   !set_diiid
!                                                                              


    RETURN
   END    SUBROUTINE set_freya_beam




   SUBROUTINE set_nubeam_misc_defaults
! -----------------------------------------------------------------------------
   USE transp


   IMPLICIT NONE



        nubeam_dt  = .025      ! time interval  between calls to nubeam 
                               ! to do beam deposition and slowing down 
                               !calcualtions

        dt_nubeam_mult = 1.0   ! fraction of dt to use for  nubeam
                               !initialization.  
                               !the time step taken on initialization of
                               !nubeam is dt_nubeam_mult*dt . 
                               !If this time step is too small nubeam will 
                               !fail (the fast ions dont reach the plasma
                               !from the source so no fast ions are 
                               !available). If this happens increase the 
                               !value of dt_numbeam_mult from 0.01 
                               !default) to some larger value 
                               !(but .le.  1.) . The reason for haveing
                               !this multiplier is that nubeam must 
                               !advance the time   whenever it is called.
                               !But on the initialization call Onetwo assumes
                               ! that time is not advanced. We assume that
                               ! dt_nubeam_mult*dt is sufficiently small
                               ! so that this time interval can be ignored.
                               ! dt_nubeam_mult should be left at 1.0 unless
                               ! you are debugging the code
 

          goocon=10.0d0 
                 ! "goose" control -- number of orbit bounces to calculate
                 ! per min([slowing-down-time],[cx-loss-time],...)
                 ! (dimensionless)

         nlfbmflr=.TRUE.
                    ! .TRUE. to accumulate the fast ion distribution 
                    ! (FBM) at the particle position (out on the FLR gyro).
                    ! if .FALSE. accumulate it at the GC.  .FALSE. may
                    ! be better for STs or other weak field large FLR
                    ! devices, because FBM does not retain gyrophase
                    ! information. 


       orbit_losses = .FALSE. 
                   !    same as nlbout in nubeam package,
                   !used by nubeam (nlbout) true if fast ion orbit losses
                   !re to be tallied in file.
                   !ilename will be RUNID(1:LRUNID)//'NB.OUT'
                   !where LRUNID is the nonblank length of RUNID.




         nlbout   =  orbit_losses 
         IF(nlbout)       &
            !get unit no for file which will hold info on orbit losses
            !file open to be done in nubeam package:
         CALL getioun(lunnbx,lunnbx) 
         CALL getioun(lunres,lunres)     !io unit required in nubeam package




        gflr_min=0.0   
                    ! for nlbgflr=.TRUE., if the guiding-center estimate
                    ! of larmor radius, in cm, is less than this, revert
                    ! to the standard guiding-center based treatment (to
                    ! save computational effort).



         nchdvp =3 !orbit integrator selector, see nbspec.dat




! Block: sys
   nubeam_nclass = 0 ! (called nclass in nubeam)
               ! set =1 to activate beam-by-beam full:half:third
               ! energy species resolved profiles of n, Pbe, Pbi,
               ! <v.B>, <(v/vcrit)v.B> (plasma frame) for NCLASS
               ! fast ion interface -- see nclass output block.


! Block: powers    ! beam powers, voltages, and energy fractions
                 ! energy fractions:  for H-isotope beams only
                 !   ffulla(j) = fraction of beam current at full energy
                 !   fhalfa(j) = fraction at half energy
                 !   1-ffulla(j)-fhalfa(j) = fraction at 1/3 energy
                 ! (for non-H beams we assume ffulla(j)=1, fhalfa(j)=0)
                 !
                 ! if f1=ffulla(j), f2=fhalfa(j), f3=1-f1-f2, and
                 ! J = the beam current for beam (j), then, the
                 ! relation  
                 !   einja(j)*J*(f1 + f2/2 +f3/3) = pinja(j)
                 ! is satisfied.
                 !
                 ! generally the powers, voltages and (perhaps) the
                 ! energy fractions should be updated every timestep.

!          RAM  einja(mib;nbeam)=0.0  ! beam voltages (eV)
!          RAM  pinja(mib;nbeam)=0.0  ! beam powers (watts)
!          RAM  ffulla(mib;nbeam)=0.98 ! full energy current fraction
                 ! (a realistic value should be provided)
!          RAM  fhalfa(mib;nbeam)=0.01 ! half energy current fraction
                 ! (a realistic value should be provided)


!  Block: fusion  ! fusion products
         nlfhe3=.FALSE.
               ! .TRUE. for He3 fusion product modeling
         nlfhe4=.FALSE.
               ! .TRUE. for He4 (alpha) fusion product modeling
         nlfst=.FALSE.
               ! .TRUE. for DD triton fusion product modeling
         nlfsp=.FALSE.
               ! .TRUE. for DD proton fusion product modeling
               ! (no recommended-- proton velocities exceed limiting
               ! assumptions in collision operator)
         nlusf3=.FALSE.
               ! .TRUE. to set He3 fusion product source magnitude from
               ! data instead of from past timestep fusion rates
        nlusfa=.FALSE.
               ! .TRUE. to set alpha fusion product source magnitude from
               ! data instead of from past timestep fusion rates
        nlusft=.FALSE.
               ! .TRUE. to set triton fusion product source magnitude from
               ! data instead of from past timestep fusion rates

        nlusfp=.FALSE.
               ! .TRUE. to set proton fusion product source magnitude from
               ! data instead of from past timestep fusion rates

        nlbfpp=.FALSE.
               ! set .TRUE. if beam slowing down is to be modeled using
               ! a fokker planck code, in which case the Monte Carlo code
               ! will be used to calculate deposition only
        plfhe3=100.0  ! (watts) threshhold power above which full Monte Carlo
               ! statistics (NPTCLF particles) are assigned to model the
               ! slowing down of the He3 fusion product ion population
        plfhe4=10000.0  ! (watts) threshhold power above which full Monte Carlo
               ! statistics (NPTCLF particles) are assigned to model the
               ! slowing down of the He4 (alpha) fusion product ion population
        xdatsf3=0.0   ! source magnitude (/sec) used to renormalize He3 fusion
               ! product source profile from previous timestep
               ! only applies if nlusf3.and.nlfhe3 both .TRUE.
        xdatsfa=0.0   ! source magnitude (n/sec) used to renormalize alpha fusion
               ! product source profile from previous timestep
               ! only applies if nlusfa.and.nlfhe4 both .TRUE.

        xdatsft=0.0   ! source magnitude (n/sec) used to renormalize triton fusion
               ! product source profile from previous timestep
               ! only applies if nlusft.and.nlfst both .TRUE.

        xdatsfp=0.0   ! source magnitude (n/sec) used to renormalize proton fusion
               ! product source profile from previous timestep
               ! only applies if nlusfp.and.nlfsp both .TRUE.
        plfst=100.0  ! (watts) threshhold power above which full Monte Carlo
               ! statistics (NPTCLF particles) are assigned to model the
               ! slowing down of the triton fusion product ion population

        plfsp=100.0  ! (watts) threshhold power above which full Monte Carlo
               ! statistics (NPTCLF particles) are assigned to model the
               ! slowing down of the proton fusion product ion population

       nlfatom=.TRUE.
               ! .FALSE. to suppress atomic physics calculations on 
               ! slowing down fusion product ions (recommendation:
               ! leave this switch alone).
       nlbfpp=.FALSE.
               ! set .TRUE. if beam slowing down is to be modeled using
               ! a fokker planck code, in which case the Monte Carlo code
               ! will be used to calculate deposition only

      nptcls=1000 !number of Monte Carlo particles per (normal) beam ion specie
                  ! affects statistics of slowing down calculation
      nptclf=1000 ! number of Monte Carlo particles per fusion product specie
                 ! affects statistics of slowing down calculation
      cxsplt=2.0  ! charge exchange model: splitting factor, tracks per event.

      dxbsmoo=0.05  ! smoothing parameter: half-width of triangular hat
                 ! function convolution smoothing, in x=sqrt(phi/philim)
                 ! applied to all output profiles
                 ! NOTE: if dxbsmoo.gt.0 a spline method will be applied
                 !   to integrated profiles; if dxbsmoo.le.0 a piecewise
                 !   linear method will be applied.  Thus:
                 !      dxbsmoo .le. 0 ==> PBI is a step function
                 !      dxbsmoo .gt. 0 ==> PBI is piecewise quadratic, C1.
                 !   which follows because the integrated profile is 
                 !   differenced to get the local value.

       ndep0=500   ! minimum number of deposition tracks per beam ion specie
                 ! affects statistics of deposition calculations

!       nbbcal=1000000  ! number of interactions for [fast ion] - [fast ion] 
                 ! fusion rate Monte Carlo integrals-- affects statistics
                 ! of this component of fusion product and neutron source
                 ! calculation, set in set_beam_defaults

       wghta=20.0d0  ! weight profile adjustment.  Larger values of WGHTA result
                 ! in more particles with less weight per particle in the core
                 ! and fewer particles with more weight per particle at the
                 ! edge.  
         !!! Larger WGHTA -> better statistics at core, worse at edge
         !!! Smaller WGHTA -> better statistics at edge, worse at core

       cxpcon=20.0d0 ! adjustment on frequency of calls to charge exchange
                 ! operator during slowing down.  Higher number --> more calls
                 ! default value generally OK.
       fppcon=8.0d0  ! adjustment on frequency of calls to collision operator
                 ! during slowing down.  Higher number --> more calls
                 ! default value generally OK.
       dtn_orbit=0.0   ! ***orbit timestep adjustment parameter***
               ! if left to zero
               ! an automatic timestep will be selected and reported-- but
               ! the automatic selection algorithm is not 100% reliable.
               ! CAUTION: if
               ! too large a timestep is used, the code will run faster but
               ! the orbit integrator will do a poor job of conserving 
               ! momentum and energy; particles could spiral in or out
               ! of the system due to this numerical error, changing the
               ! results of the physical simulation.  If too small a timestep
               ! is used, the code will run very slowly.
               ! (seconds)
       xdepmod=1.0       ! anomolous opacity adjustment
               ! this applies a multiplicative factor to all stopping
               ! cross sections.
       nsigexc=0 ! set =1 toturn on approximate excited states correction model
               ! Janev/Boley/Post article, Nucl. Fusion Vol.29, No. 12 (1989).

       nlminsv=.TRUE.
               ! .TRUE. to loop over all impurities when calculating
               ! impurity stopping cross section; .FALSE. to use an
               ! "average Z" single impurity simplifying approximation.

       nlbbcx=.TRUE.
               ! .FALSE. to disable beam-beam charge-exchange and
               ! beam-beam impact-ionization contribution to deposition.
               ! This concerns beam neutrals being stopped by collisions
               ! with slowing down beam ions.
       nlebei=.TRUE.
               ! .TRUE. to include beam energy in approximate
               ! expression for <sigma*v> for electron impact ionization

       dn0out=5.0d11
               ! neutral gas density (/cm3) to assume beyond plasma boundary
               ! (drives "EXTERNAL charge exchange loss")
       nmsigx =2  ! control option: select cross section set for impurity 
               ! stopping of fast neutrals (deposition and cx recapture)
               ! 1: Olson (very ancient), H beams only 
               !    OLSON ET AL., P.R.L. 41,NUM 3,P. 163.
               ! 2: Standard set
               !    Currently:  
               !    R.A. Phaneuf, R.K. Janev, M.S. Pindzola, 
               !    Atomic Data for Fusion, Vol. 5,
               !    Collisions of Carbon and Oxygen Ions with 
               !    Electrons, H, H2, and He, 
               !    ORNL-6090/V5 (5th volume of ORNL-6086) (Feb. 1987).
       xcfanbi=1.0  ! anomolous multiplier: electron drag and energy
                  ! diffusion on beam ions
       xdfanbi=1.0  ! anomolous multiplier: ion drag and energy
                  ! diffusion on beam ions
       xefanbi=1.0  ! anomolous multiplier: pitch angle scattering on
                  ! beam ions

       xcfafus=1.0  ! anomolous multiplier: electron drag and energy
                  ! diffusion on fusion product ions
       xdfafus=1.0  ! anomolous multiplier: ion drag and energy
                  ! diffusion on fusion product ions
       xefafus=1.0  ! anomolous multiplier: pitch angle scattering on
                  ! fusion product ions
       nlbcde=.TRUE.  ! .FALSE. to turn off energy diffusion
       nlbcoh=.TRUE.  ! .FALSE. to suppress acceleration of beam ions by
                    ! toroidal electric field
       nlbcpa=.TRUE.  ! .FALSE. to turn off pitch angle scattering

       nlorbo=.FALSE. ! .TRUE. for "orbit ONLY" mode-- suppress all collisions
                    ! used only in specialized debugging situations.

! Block: flr     ! finite larmor radius corrections
       nlbflr=.TRUE.  ! .TRUE. for standard FLR treatment, based on a mod(B)
                    ! at guiding center and notion of a FLR circle about
                    ! the guiding center, upon which a random gyro-phase
                    ! angle is taken to define the current particle position.
       nlbgflr=.FALSE.  ! .TRUE. for enhanced FLR treatment
                    ! this deforms the FLR trajectory by looking at the
                    ! variation of mod(B) in the guiding center vicinity,
                    ! and is appropriate when the larmor radius is >~ L(B)
                    ! e.g. in low field systems such as STs.  **but** this
                    ! option is extremely expensive computationally!
       gflr_op  =   &
             (/0.,1.e-4,1.e-5,0.,3.,10.,1.,0.2,0.05,1.,2.,0.,0.,0.,0./)
                    ! nlbgflr=.TRUE. expert controls (do not touch)

! Block: adif    ! anomolous diffusion
               ! if this is used, an anomolous diffusivity profile must be
               ! provided via xplasma (see section on input data profiles).
               ! if no profile is provided, the anomolous diffusivity will
               ! be zero and this feature will be inactive.

       nkdifb=3  ! =3: apply to all fast ions;
               ! =1: apply to beam ions only;
               ! =2: apply to fusion product ions only.
! Block: ripple  ! TF field ripple effects.  The available models are crude.
               ! A TFTR project to build a better model ran out of time and
               ! money.

       nrip=0    ! =0: no ripple loss;  =1: simple "loss time" model
               ! =2: adjustable model based on Goldston-White-Boozer and
               !     ripple field strengh. CAUTION-- calibration against
               !     a more exact model is needed.  PPPL expert help should
               !     be sought before using this model

       taurip=0.0  ! must be set if nrip=1: loss time (sec) for trapped ions
       ncoils=0    ! number of TF coils.  Must be set if nrip=2.
       asrd=0.0    ! nrip=2 adjustment factor for fusion product ions
       bsrd=0.0    ! nrip=2 adjustment factor for beam ions.

! Block: outcon  ! code output control options
       nsdbgb=1  ! index to fast ion specie selected for debug output
                ! (artifact-- not currently used).

! Block: fishbone  ! fishbone model (caution-- artifact, may not work)
       nlfbon=.FALSE.  ! .TRUE. to turn on ancient PBX fishbone loss model
               ! "fishbones" were an MHD phenomenon seen on PPPL's PBX
               ! tokamak in the 1980s, which appeared to cause fast ion loss

       fbemin=0.0      ! minimum energy (eV) of affected ions
       fbemax=20.0d4   ! maximum energy (eV) of affected ions
       fvpvmn=-0.1d0   ! minimum vpll/v of affected ions
       fvpvmx= 0.1d0   ! maximum vpll/v of affected ions
       fbltim= 0.05d-3 ! average loss time (seconds) of affected ions
       fshper= 5.0d-3  ! fishbone period (cycle time, seconds)
       fshwid= 1.0d-3  ! fishbone duration (seconds)
       tfshon= -1.0e34     ! fishbone onset time (seconds) default always on
       tfshof= +1.0e34     ! fishbone stop time (seconds) default always on

! Block: box     ! beam-in-box neutral density calculation controls
       nxbox=1   ! number of bins in horizontal direction (x) in front of beam

       nybox=1   ! number of bins in vertical direction (y) in front of beam

       nlbox=100 ! number of bins in the direction (l) of beam injection

       xboxhw=5.0  ! (cm) half-width of entire cartesian box in x direction
                 ! individual bins' x width = 2*xboxhw/nxbox

       yboxhw=5.0  ! (cm) half-height of entire cartesian box in y direction
                 ! individual bins' y width = 2*yboxhw/nybox

       xlbox1=100.0  ! (cm) distance from beam source at which box starts
       xlbox2=500.0  ! (cm) distance from beam source at which box ends
               ! individual bins' l length = (xlbox2-xlbox1)/nlbox

       ndepbox=500
               ! Monte Carlo statistics control:  minimum number of neutral
               ! tracks to use per calculation.

!    Block: minority  ! RF minority speciess
                 ! in the course of a run, RF species can be added
                 ! but once added can never be taken away

      nmini=0     ! number of RF minority species
      
!    RAX  xzmini(mmini;nmini)=0.0  ! Z of RF minority specie (atomic no.)
!    RAX  amini(mmini;nmini)=0.0   ! A of RF minority specie (AMU)

      nthrf=1
                 ! = 1 if minority ions are currently thermal

!  Block: misc    ! miscellaneous

      plhgt=0.0 ! (updown symmetric geometries only) height of plasma 
               ! midplane above/below vacuum vessel midplane
               ! neutral beam component (e.g. ion source, aperture) 
               ! elevations are given with respect to the vacuum vessel
               ! midplane.  (In updown asymmetric cases the vertical
               ! offset is presumed to be given in the MHD equilibrium).
               !   **time dependent quantity**
       nlcprb=.TRUE.
               ! set .FALSE. to suppress compression operator and other
               ! effects of time varying MHD equilibrium on beam ions
               ! (recommendation: leave default, .TRUE., alone).

       edbfac=10.0
               ! debug tripwire -- if an ion reaches edbfax x its initial
                ! deposition energy, kill the run.
       nmcurb=1  ! set to 0 to suppress fast ion driven current.
               ! nmcurb=1 or 2 gives an aspect-ratio approximate 
               ! neoclassical shielding calculation for the current;
               ! nmcurb=3 gives a simple "Spitzer" shielding calculation.
               ! If in the future a fully general geometry shielding
               ! calculation is available this will likely be offered as
               ! nmcurb=4 (but not yet).

       nlfdep=.FALSE.
               ! set .TRUE. to compute midplane index oriented deposition
               ! distribution function in the lab frame using the "fbm"
               ! grids for pitch and energy, with particles mapped to
               ! their outer midplane crossing (a few particles which
               ! never cross the midplane at R > Raxis will be missed!).
               ! setting this switch causes a specialized, rarely used
               ! output data structure to be created
       nlsym2b =.FALSE. ! reset by isym below     
               ! set to true  for updown symmetery
       nzones  = 2* inzri !max # zones for fast ion calcs
                     !must be a multiple of inzri
  
!  Block: saw     ! sawtooth model controls
       ! note: on every fast ion step, sawtooth mixing is computed, but
       ! not applied.  If in fact a sawtooth occurs, a call to "sawnbi"
       ! applies the sawtooth to the Monte Carlo population.
       !
       ! The "sawnbi" call cannot change the fast ion density and energy
       ! density profiles in the driver code; instead, these profiles are
       ! computed as a part of each regular fast ion timestep.  It is up
       ! to the driver code to save these profiles 

      nlsawb=.TRUE.  ! if .TRUE. and there is a sawtooth, beam ions are mixed.
      nlsawf=.TRUE.  ! if .TRUE. and there is a sawtooth, fusion product ions
               ! are mixed.
!   Block ?? 
!     rhmin(mj, miz,mibs)
        ! DENSITY OF ICRF MINORITY ION SPECIES (N/CM**3)
!     den0mn(mig, mj, mis)
        ! NEUTRAL DENSITY OF HYDROGENIC SPECIES DUE TO SOURCE
!   en0mn(mig, mj, mis)
        ! NEUTRAL temperature OF HYDROGENIC SPECIES DUE TO SOURCE
!   omg0mn(mig, mj, mis)
           ! NEUTRAL rotation velocity OF HYDROGENIC SPECIES DUE TO SOURCE
        nlntmj= .FALSE.   !use time smoothed majority temperature ?? (
           ! (not defined in transp documentation )
!       tmjsm(mj)
           ! time smoothed majority temperature.


!       rhi(mj, miz)
            ! TOTAL IMPURITY DENSITY AT BDY/ZONE J
!       xzimpj(mj)
        ! Average impurity atomic number by zones
!       aimpj(mj)
        ! Average impurity atomic weight by zones
        !onetwo note: aimpj is r* on nj grid, converted
        !to nubeam requirements in set_nbprofs

!      nlfbmfpp(mibs)
!      nsbeams     not defined in transp

!      mibs = 7:15
              ! max total number of fast ion species including beam ions,
              ! fusion product ions, and RF minority ions
!      xzbeams(mibs)
!      abeams(mibs)
!      nlfprod(mibs)
!      rhbs(mj, miz, mibs)
        ! BEAM DENSITY SPECIES IBS - CURRENT TI STEP
!      omegag(mj)
        ! ROTATION ANGULAR VELOCITY AT TGEO (RADIANS/SEC)
        !Onetwo note: omegag(xp_nj) is r*8, converted to
        !required nubeam form in sub set_nbprofs
!      phiprg(mj)
        ! radial electrostatic potential ???? (incorrectly
        ! defined in nbspec.dat)
        !Onetwo note: phiprg(xp_nj) is r*8, converted to
        !required nubeam form in sub set_nbprofs
!      curt(mj)
        ! AREA INTEGRAL TOROIDAL CURRENT, AMPS, TO BDY J
        !Onetwo note: curt(xp_nj) is r*8, converted to
        !required nubeam form in sub set_nbprofs
!      vpoh(mj)
        ! LOOP VOLTAGE as used for POH calculation, zone J
        !Onetwo note: vpoh(xp_nj) is r*8, converted to
        !required nubeam form in sub set_nbprofs

! ng    !number of gases

!      rhob(mj, mig, miz)
        ! BACKGROUND PLASMA DENSITY AT BDY/ZONE J
        !Onetwo note: rhob  is directly replace with en,see
        !sub set_nbprofs
!      rhix(mj, miz, mimpt)
        ! MULTIPLE IMPURITY DENSITY AT BDY/ZONE J USED IN CODE


   RETURN
   END    SUBROUTINE set_nubeam_misc_defaults






SUBROUTINE nbdrive_nubeam_checknaml(ierr)

  USE transp
  USE io,ONLY : lun_msgs => ncrt
  IMPLICIT NONE
  REAL*8, PARAMETER :: CZERO = 0.0d0
  INTEGER, INTENT(out) :: ierr

  !  perform additional namelist checks.  Call a NUBEAM routine to
  !  pre-fetch the fast ion distribution function grids.

  !  The following checks are performed:
  !
  !  a) no duplicates in thermal species lists
  !  b) for each beam injected specie a matching thermal specie must
  !     be present
  !  c) for each fusion product specie a matching thermal specie is
  !     not required, but a warning will be issued.

  !
  !  call ancillary NUBEAM routines to check ion species lists; return 
  !  indices to specific isotopes, total no. of fast ion species, error 
  !  flags if duplicate species are given, etc.

  !  the NUBEAM calls in this routine do not affect the NUBEAM module's
  !  internal state-- no NUBEAM f90 modules or COMMONs are used-- 
  !  and they can be called at any time.  They provide information
  !  which is available both to NUBEAM and to NUBEAM driver codes, 
  !  and which depends only on the routines' input arguments.

  !-------------------------------
  INTEGER iwarn, ierr_p, i, j, iA, iZ, iA2, iZ2, imatch, ilin,ngmax,nsfast
  REAL*8 zelin,zeinjmax
  !-------------------------------




  ierr = 0
  ngmax = beam_data%ngmax

  !-------------------------------------------------------------
  !  check thermal ion species list...
  !  input: species list; output: indices, error/warning flags
  !--------------------------------------------------------------
  DO i=1,ngmax-1
     iA =  beam_data%aplasm(i) + 0.1d0
     iZ =  beam_data%backz(i) + 0.1d0
     DO j=i+1,ngmax
        iA2 = beam_data%aplasm(j) + 0.1d0
        iZ2 = beam_data%backz(j) + 0.1d0
        IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) ierr=ierr+1
     ENDDO
  ENDDO

  IF(ierr.NE.0) THEN
     WRITE(lun_msgs,*) ' ?duplicate species in thermal species list:'
     WRITE(lun_msgs,*) '  beam_data%aplasm(1:ngmax) = ',beam_data%aplasm(1:ngmax)
     WRITE(lun_msgs,*) '   beam_data%backZ(1:ngmax) = ', beam_data%backZ(1:ngmax)
     CALL STOP( 'nbdrive_nubeam_checknaml',1)
  ENDIF

  !  do same for impurities

  DO i=1,beam_data%nrhix-1
     iA = beam_data%aimpx(i)
     iZ = beam_data%xzimpx(i)
     DO j=i+1,beam_data%nrhix
        iA2 = beam_data%aimpx(j)
        iZ2 = beam_data%xzimpx(j)
        IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) ierr=ierr+1
     ENDDO
  ENDDO

  IF(ierr.NE.0) THEN
     WRITE(lun_msgs,*) ' ?duplicate species in impurity species list:'
     WRITE(lun_msgs,*) '  Aimpx(1:beam_data%nrhix) = ', beam_data%Aimpx(1:beam_data%nrhix)
     WRITE(lun_msgs,*) '  xZimpx(1:beam_data%nrhix) = ',beam_data%xZimpx(1:beam_data%nrhix)
     CALL STOP( 'nbdrive_nubeam_checknaml',2)
  ENDIF

  !------------------------------------------------------------------
  !  build list of fast species from beam species...
  !  NOTE beam_data%abeams and beam_data%xzbeams are not passed to nubeam
  !  Instead beam_data%abeama and beam_data%xzbeama are passed
  !  This is due to the fact that nubeam constructs abeams and xzbeams
  !  Internally. We keep the local results in beam_data so that we may compare 
  !  with the output  of nubeam in sub load_12.
  !-------------------------------------------------------------------
 
  nsfast=0
  DO i=1,beam_data%nbeam              !each beam may be a different species
     iA=beam_data%abeama(i)           !but beams are also distinguished
                                      !
     iZ=beam_data%xzbeama(i)
     imatch=0
     DO j=1,nsfast
        iA2=beam_data%abeams(j)
        iZ2=beam_data%xzbeams(j)
        IF(iZ2.GT.2) THEN
           WRITE(lun_msgs,*) ' ?nbdrive: beams must be isotopes of H or He'
           ierr = ierr+1
           imatch = 99
           EXIT
        ELSE IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) THEN
           imatch=j
           EXIT
        ENDIF
     ENDDO

     IF(imatch.EQ.0) THEN

        nsfast = nsfast + 1
        beam_data%abeams(nsfast) = beam_data%abeama(i)
        beam_data%xzbeams(nsfast) = beam_data%xzbeama(i)

        !  require that corresponding thermal specie exist

        imatch = 0
        DO j=1,ngmax
           iA2 = beam_data%aplasm(j) + 0.1d0
           iZ2 = beam_data%backz(j) + 0.1d0
           IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) THEN
              imatch=j
              EXIT
           ENDIF
        ENDDO

        IF(imatch.EQ.0) THEN
           WRITE(lun_msgs,*) ' ?nbdrive_nubeam_checknaml: beam no. ',i
           WRITE(lun_msgs,*) '  injects specie A=',beam_data%abeama(i),' Z=',beam_data%xzbeama(i)
           WRITE(lun_msgs,*) '  but there is no such thermal specie:'
           WRITE(lun_msgs,*) '  Aplasm(1:ngmax) = ',beam_data%aplasm(1:ngmax)
           WRITE(lun_msgs,*) '   backZ(1:ngmax) = ', beam_data%backZ(1:ngmax)
           ierr = ierr+1
           EXIT
        ENDIF
     ENDIF
  ENDDO


  IF(ierr.GT.0) THEN
     WRITE(lun_msgs,*) &
          ' NUBEAM nbi_datchkb detected problem with beam species.'
     CALL STOP('nbdrive_nubeam_checknaml',1)
  ENDIF

   beam_data%nbeam_species = nsfast     ! number of unique  injected beam species
                                        ! regardless of beam number

  !-------------------------------
  !  check fusion product species
  !  also compute array of maximum energies for distribution functions.
  zeinjmax = MAXVAL(beam_data%einja(1:MAX(1,beam_data%nbeam)))
  beam_data%ebdmaxa = CZERO
  IF(beam_data%ebdmax .LE. 0.0)beam_data%ebdmax = 1.1*zeinjmax
  beam_data%ebdmaxa(1:nsfast) = beam_data%ebdmax   ! max energy for BEAM species (eV) from namelist
   beam_data%nznbmea = 0
   beam_data%nznbmea(1:nsfast) = beam_data%nznbme   ! no. of energy zones for BEAM species dist fcns

    beam_data%nlfprod(:) = 0            ! fusion product species flags

  !  input: flags for specific fusion product species
  !         information arrays for beam species
  ! output: information arrays extended to cover fusion product species
  !         indices for specific fusion product isotopes

  IF(nlfsp) THEN
     WRITE(lun_msgs,*) ' ?nbdrive_nubeam_checknaml: nlfsp = .TRUE.'
     WRITE(lun_msgs,*) '  NUBEAM needs upgrade to be able to follow'
     WRITE(lun_msgs,*) '  fusion product protons.'
     ierr = ierr + 1
  ENDIF

  !  note setting of maximum energy for fusion product distribution
  !  functions (eV) using NUBEAM utility call "r8_nbfusn_emax(iZ,iA,emax)"

  IF(nlfst) THEN
     nsfast = nsfast + 1
     beam_data%xzbeams(nsfast)=1
     beam_data%abeams(nsfast)=3
     CALL r8_nbfusn_emax(1,3,beam_data%ebdmaxa(nsfast))
                        ! safely above the 1.01 MeV fusion T birth energy
      beam_data%nznbmea(nsfast)=beam_data%nznbme
      beam_data%nlfprod(nsfast)= 1
  ENDIF

  IF(nlfhe3) THEN
     nsfast = nsfast + 1
     beam_data%xzbeams(nsfast)=2
     beam_data%abeams(nsfast)=3
     CALL r8_nbfusn_emax(2,3,beam_data%ebdmaxa(nsfast))
                        ! safely above the 0.82 MeV fusion He3 birth energy
      beam_data%nznbmea(nsfast)=beam_data%nznbme
      beam_data%nlfprod(nsfast)= 1
  ENDIF

  IF(nlfhe4) THEN
     nsfast = nsfast + 1
     beam_data%xzbeams(nsfast)=2
     beam_data%abeams(nsfast)=4
     CALL r8_nbfusn_emax(2,4,beam_data%ebdmaxa(nsfast))
                        ! safely above the alpha (He4) birth energy
                        ! for any fusion reaction
     beam_data%nznbmea(nsfast)=beam_data%nznbme
     beam_data%nlfprod(nsfast)= 1
  ENDIF

   beam_data%nfusion_species = nsfast - beam_data%nbeam_species





  !-------------------------------
  !  there is also an r8_nbi_datckrf call for RF minorities -- not
  !    currently used in nbdrive.
  !-------------------------------
  !  allocate energy grids for the fast ion distributions; call a NUBEAM
  !  routine to fill these in.  Although not necessary to make NUBEAM work,
  !  this allows the NUBEAM caller to see the distribution function grids
  !  ahead of time.
  !    note that the above calls to "r8_nbfusn_emax(iZ,iA,emax)" to set
  !    the max energy.

  !  in nbdrive, all species distribution functions have the same number
  !  of energy grid points, but different maximum energies.  In general,
  !  it is allowed for the number of energy grid points to vary depending
  !  on fast specie, however.

  ALLOCATE(beam_data%efbm(beam_data%nznbme,nsfast),beam_data%efbmb(beam_data%nznbme+1,nsfast))
  beam_data%efbm=CZERO; beam_data%efbmb=CZERO
  DO i=1,nsfast
     CALL r8_nbi_efbm(beam_data%nlfprod(i), beam_data%nznbmea(i),beam_data%ebdmaxa(i), &   ! ** NUBEAM CALL **
          ilin,zelin,beam_data%efbmb(1: beam_data%nznbmea(i)+1,i), &
          beam_data%efbm(1: beam_data%nznbmea(i),i))
     ! (r8_)nbi_efbm outputs:
     ! efbm -- grid zone centers
     ! efbmb -- grid boundaries
     ! ilen,zelin -- index and energy to which grid is linearly (rather
     ! than logarithmically) spaced.
  ENDDO




  !-------------------------------
  !  the following is for the benefit of nbdrive output routines.
  !  it is placed here because the fast ion species list has now
  !  been defined.
  !-------------------------------
  !  now gen some labels for species
  
  IF(ierr.EQ.0) THEN
     ALLOCATE(beam_data%label_thions(beam_data%ngmax))
     DO i=1,beam_data%ngmax
        CALL zlabel( beam_data%backz(i), beam_data%aplasm(i),'thermal plasma ion',beam_data%label_thions(i))
     ENDDO

     ALLOCATE(beam_data%label_impions(beam_data%nrhix))
     DO i=1,beam_data%nrhix
        CALL zlabel( beam_data%xzimpx(i), beam_data%aimpx(i),'impurity plasma ion',beam_data%label_impions(i))
     ENDDO

     ALLOCATE(beam_data%label_fastions(nsfast))
     DO i=1,nsfast
        IF(beam_data%nlfprod(i) .gt. 0 ) THEN
           CALL zlabel( beam_data%xzbeams(i), beam_data%abeams(i),'fusion product ion', &
                beam_data%label_fastions(i))
        ELSE
           CALL zlabel( beam_data%xzbeams(i), beam_data%abeams(i),'beam ion', &
                beam_data%label_fastions(i))
        ENDIF
     ENDDO

     ALLOCATE(beam_data%label_beamions(beam_data%nbeam_species))
     i=0
     DO j=1,nsfast
         IF(beam_data%nlfprod(j) .le. 0 )THEN
           i=i+1
           if(i .gt. beam_data%nbeam_species)       &
                CALL STOP("checknaml species count",1)

           IF(beam_data%abeams(j) == 1.0 )THEN
              beam_data%label_beamions(i) = 'h' !note lower case for onetwo
           ELSE IF(beam_data%abeams(j) == 2.0)THEN
              beam_data%label_beamions(i) = 'd'
           ELSE
              CALL STOP("checknaml species type",1)
           ENDIF
        ENDIF
     ENDDO
  ENDIF

  IF(ierr.EQ.0) THEN
     WRITE(lun_msgs,*) ' ---> nbdrive_nubeam_checknaml: success.'
  ENDIF

  CONTAINS
    SUBROUTINE zlabel(zz,aa,suffix,label)

      USE periodic_table_mod

      !  generate labels for ion species

      REAL*8, INTENT(in) :: zz,aa   ! Z & A of species
      CHARACTER*(*), INTENT(in) :: suffix   ! label suffix
      CHARACTER*(*), INTENT(out) :: label   ! generated label.

      REAL*4 aaa
      INTEGER iz,ia

      CHARACTER*12 prefix

      !----------------

      aaa = SNGL(aa)
      iz = zz + 0.1
      ia = aa + 0.1

      prefix='?'
      IF(iz.EQ.1) THEN
         IF(ia.EQ.1) THEN
            prefix='H'
         ELSE IF(ia.EQ.2) THEN
            prefix='D'
         ELSE IF(ia.EQ.3) THEN
            prefix='T'
         ENDIF
      ELSE IF(iz.EQ.2) THEN
         IF(ia.EQ.3) THEN
            prefix='He3'
         ELSE IF(ia.EQ.4) THEN
            prefix='He4'
         ENDIF
      ELSE

         prefix = to_periodic_table(iz,aaa,-1,0)

      ENDIF

      label = TRIM(prefix)//' '//TRIM(suffix)

    END SUBROUTINE zlabel
END SUBROUTINE nbdrive_nubeam_checknaml








    SUBROUTINE ufile_bypass(read_ufile)
! ------------------------------------------------------------
!   subroutine assumes beam_data_namelist was read and beam 
!   data sructure was loaded. Here we check the quantities 
!   nbeam, tbona(1:nbeam),tboffa(1:nbeam),
!   pinja(1:nbeam),einja(1:nbeam), ffulla(1:nbeam) and
!   fhalfa(1:nbeam)
!   for valid input. If this input is valid for all beams then set flag 
!   read_ufile = .false., indicating that the ufile should not be processed.
!   Otherwise read_ufile = .true. is returned and the ufile will have to
!   be read to get the inputs.
!
!   beam_data%switch_on_time, 
!
! -------------------------------------------------------------HSJ 12/03/03
    USE transp,ONLY : beam_data,use_nubeam,nubeam_restart,nubeam_dtmin
    USE solcon,ONLY : time0,timmax
    USE io,    ONLY : nout,ncrt,nitre
    IMPLICIT NONE
    INTEGER j,k,l,ll,nt,nchan,nx,ierr
    REAL*8 tb,tmin
    REAL *8,DIMENSION(:),ALLOCATABLE :: dummy
    LOGICAL,INTENT(OUT) :: read_ufile
    INTEGER kk,km,i,kmm


   INTERFACE local
     SUBROUTINE sort_unique(a,km)
       REAL *8, INTENT(INOUT),DIMENSION(:) :: a
       INTEGER, INTENT(OUT) :: km
     END SUBROUTINE sort_unique
   END INTERFACE !jmp.ibm

   INTERFACE local2 !jmp.ibm
     SUBROUTINE REALLOCATE1DD(A,KM)  
       IMPLICIT NONE
       REAL *8, DIMENSION(:),POINTER :: a
       INTEGER, INTENT(IN) :: km
      END SUBROUTINE REALLOCATE1DD
   END INTERFACE

    read_ufile = .FALSE.

    !allow  0.0 as valid pinj
    !allow fhalfa = 0.0 for neg ion injection
    DO j =1,beam_data%nbeam
       IF (      beam_data%pinja(j) .LT. 0.0           &
            .OR. beam_data%einja(j) .LE. 0.0           &
            .OR. beam_data%ffulla(j) .LE. 0.0          )THEN
           read_ufile = .TRUE.
           EXIT
       ENDIF
       IF( beam_data%tbona(j) .GE.  beam_data%tboffa(j))THEN
          read_ufile = .TRUE.
          write(ncrt,110)j,beam_data%tbona(j),beam_data%tboffa(j)
          write(nout,110)j,beam_data%tbona(j),beam_data%tboffa(j)
          write(nitre,110)j,beam_data%tbona(j),beam_data%tboffa(j)
110       format(2x,'Error for beam number ',i5,/,2x,'tbona > tboffa',/, &
                   ' required, rbona ,tboffa =',2(2x,1pe12.6))
          EXIT
       ENDIF

    ENDDO

    IF(read_ufile)RETURN              !data not valid


    !WE HAVE valid data, load beam_data structure:
    nchan = beam_data%nbeam * 4               !4 data elements per beam


    IF(ASSOCIATED( beam_data%beam_chan))THEN
                DEALLOCATE(beam_data%beam_chan)
                ALLOCATE(beam_data%beam_chan(nchan))
                DO j=1,nchan
                   beam_data%beam_chan(j) = j
                ENDDO
    ENDIF



    !check each tbona,tboffa 
    !the values in tbona,tboffa may not be unique
    !we turn them into unique arrays here by loading
    !the unique values into switch_on_time and switch_off_time

    k = SIZE(beam_data%tbona)
    IF(SIZE(beam_data%switch_on_time).LT. k )                    &
       CALL REALLOCATE1DD(beam_data%switch_on_time,k)
    beam_data%switch_on_time(:) = beam_data%tbona(:)
    !sort switch on times,eliminate duplicate times:
    CALL sort_unique(beam_data%switch_on_time,km)
    CALL REALLOCATE1DD(beam_data%switch_on_time,KM)
    !PRINT *,'switch_on =',beam_data%switch_on_time


    k = SIZE(beam_data%tboffa)
    IF(SIZE(beam_data%switch_off_time).LT. K)                     &
       CALL REALLOCATE1DD(beam_data%switch_off_time,k)
    beam_data%switch_off_time(:) = beam_data%tboffa(:)
    !sort switch on times,eliminate duplicate times:
    CALL sort_unique(beam_data%switch_off_time,km)
    k=SIZE(beam_data%switch_off_time)
    IF( km .lt. k)                                                  &
           CALL REALLOCATE1DD(beam_data%switch_off_time,KM)
    !PRINT *,'switch_off =',beam_data%switch_off_time


    IF(ASSOCIATED( beam_data%beam_times)) DEALLOCATE(beam_data%beam_times)
    IF(ASSOCIATED(beam_data%beam_inject)) DEALLOCATE(beam_data%beam_inject)
    tmin = MINVAL(beam_data%switch_on_time)
!    if(nubeam_restart .gt. 0)then
!       tmin = time0
!       DO j=1,SIZE(beam_data%switch_on_time)
!          IF(beam_data%switch_on_time(j) .lt. tmin)THEN
!             beam_data%switch_on_time(j) = tmin
!          ENDIF
!       ENDDO
!    endif
    !if beam comes on right at time0 reset beam on time to 
    !allow for small beam risetime (THIS is not for physics
    !reasons. It is to incorporate linear in time interpolation
    !of beam power. note that beam_power_rise_time =0.0 is valid only
    !if tmin = time0).
    IF(tmin .LT. time0 .and. nubeam_restart == 0)THEN                                !beamon before time0
       PRINT *,'beam is not allowed to be on before time0'
       PRINT *,'time0, beam on time =',time0,tmin
       CALL STOP('ufile_bypass, beamon < time0',1)
    ELSE                                                    !beamon after time0
       IF(tmin-time0 .LT. beam_data%beam_power_rise_time)THEN
          beam_data%beam_power_rise_time = 0.5 * (tmin-time0)
          beam_data%beam_power_rise_time =                        &
              MAX(beam_data%beam_power_rise_time,nubeam_dtmin)
       ENDIF
       k = SIZE(beam_data%switch_off_time)+SIZE(beam_data%switch_on_time)
       ALLOCATE(dummy(k))
       DUMMY(1:SIZE(beam_data%switch_on_time))=  beam_data%switch_on_time(:)
       DUMMY(SIZE(beam_data%switch_on_time)+1:k) = beam_data%switch_off_time(:)

       CALL sort_unique(dummy,km)

       nx = 3*km + 2

       ALLOCATE(beam_data%beam_times(nx))
       ALLOCATE(beam_data%beam_switch_times(SIZE(dummy)))
       l=1
       ll =0
       beam_data%beam_times(l) = time0
       DO j=1,SIZE(dummy)
          tb = dummy(j)
          IF(tb .GE. timmax)EXIT
          DO k=1,SIZE(beam_data%switch_on_time)
             IF(tb .EQ. beam_data%switch_on_time(k))THEN
                l=l+1
                beam_data%beam_times(l) = tb - beam_data%beam_power_rise_time
                l=l+1
                beam_data%beam_times(l) = tb

                l=l+1                             !3/13
                beam_data%beam_times(l) = tb +  beam_data%beam_power_rise_time

                ll =ll+1
                beam_data%beam_switch_times(ll) = tb
             ENDIF
          ENDDO
          beam_data%beam_power_decay_time = beam_data%beam_power_rise_time
          DO k=1,SIZE(beam_data%switch_off_time)
             IF(tb .EQ. beam_data%switch_off_time(k))THEN
                l=l+1
                beam_data%beam_times(l) = tb
                ll =ll+1
                beam_data%beam_switch_times(ll) = tb
                l=l+1
                beam_data%beam_times(l) = tb + beam_data%beam_power_decay_time
                l=l+1                              !3/13
                beam_data%beam_times(l) = tb + 2.*beam_data%beam_power_decay_time

             ENDIF
          ENDDO
       ENDDO 
       l=l+1
       beam_data%beam_times(l) = timmax
       if(l .lt. SIZE(beam_data%beam_times))          &
                     CALL REALLOCATE1DD(beam_data%beam_times,l)
      !duplicate in beam_data%beam_times  may still exist due to duplicates in 
      !beam_data%switch_off_time and beam_data%switch_on_time. Eliminate them
      !now:
    CALL sort_unique(beam_data%beam_times,km)
    iF(km .lt. SIZE(beam_data%beam_times))          &
                     CALL REALLOCATE1DD(beam_data%beam_times,km)


       ALLOCATE(beam_data%beam_inject(km,4*beam_data%nbeam))
       DO j=1,SIZE(beam_data%beam_times)
          DO k=1,beam_data%nbeam
             beam_data%beam_inject(j,beam_data%nbeam+k)   = beam_data%einja(k)
             beam_data%beam_inject(j,2*beam_data%nbeam+k) = beam_data%ffulla(k)
             beam_data%beam_inject(j,3*beam_data%nbeam+k) = beam_data%fhalfa(k)
             IF(beam_data%beam_times(j) .LT. beam_data%tbona(k)      &
                     -0.5 * beam_data%beam_power_rise_time           &
                       .OR. beam_data%beam_times(j) .GT.             &
                      0.5*beam_data%beam_power_decay_time +          &
                                                beam_data%tboffa(k))THEN
                beam_data%beam_inject(j,k) = 0.0
             ELSE
                beam_data%beam_inject(j,k) =  beam_data%pinja(k)
             ENDIF
          ENDDO
       ENDDO
    ENDIF


    DEALLOCATE(dummy)
    DEALLOCATE(beam_data%switch_off_time)
    DEALLOCATE(beam_data%switch_on_time)

    CALL  nbdrive_nubeam_checknaml(ierr)

    IF( .NOT. use_nubeam) CALL nfreya_load  
 
    RETURN
    END  SUBROUTINE ufile_bypass        





   
    SUBROUTINE  REALLOCATE1DD(A,KM)  
! REALLOCATE 1DIM DOUBLE ARRAY A
! CHANGE SIZE OF ARRAY A TO KM
      IMPLICIT NONE
      REAL *8, DIMENSION(:),POINTER :: a
      REAL *8, DIMENSION(:),ALLOCATABLE :: dummy
      INTEGER, INTENT(IN) ::  km
      INTEGER sizea
      sizea  = SIZE(a)
      ALLOCATE(dummy(sizea))
      dummy(:) = a(:)
!      DEALLOCATE(a)   !deallocastion of pointer not allowed in lf95
                       !according to f90 nadbook because a is a pointer array
                       !we do not have to deallocate it.
                       !(pg 181,Adams ,et al "Fortran 95 ahndbook")
      ALLOCATE(a(km))
      IF(KM .GT. sizea)then   !increase size of a
         a(1:sizea) = dummy(1:sizea)
         a(sizea+1:km) =0.0d0
      ELSE                     !decrease size of a
         a(1:km) = dummy(1:km)
      ENDIF
      DEALLOCATE(dummy)
      RETURN
      END





   
    SUBROUTINE  REALLOCATE1DD0(A,KM)  
! REALLOCATE 1DIM DOUBLE ARRAY A, starting at 0
! CHANGE SIZE OF ARRAY A TO KM + 1
      IMPLICIT NONE
      REAL *8,  DIMENSION(:),POINTER :: a
      REAL *8,DIMENSION(:),ALLOCATABLE :: dummy
      INTEGER km,sizea
      sizea  = SIZE(a)
      ALLOCATE(dummy(0:sizea-1))
      dummy(:) = a(:)
      DEALLOCATE(a)
      ALLOCATE(a(0:km))
      IF(KM .GT. sizea)then
         a(0:sizea-1) = dummy(0:sizea-1)
         a(sizea:km) =0.0d0
      ELSE
         a(0:km-1) = dummy(0:km-1)
      ENDIF
      DEALLOCATE(dummy)
      RETURN
      END



   
    SUBROUTINE  REALLOCATE1DS(A,KM) 
!-------------- I see little value in function Overloading  --- HSJ 
!
! REALLOCATE 1DIM SINGLE  ARRAY A
! CHANGE SIZE OF ARRAY A TO KM
      IMPLICIT NONE
      REAL *4,  DIMENSION(:),POINTER :: a
      REAL *4,DIMENSION(:),ALLOCATABLE :: dummy
      INTEGER km,sizea
      sizea  = SIZE(a)
      ALLOCATE(dummy(sizea))
      dummy(:) = a(:)
      DEALLOCATE(a)
      ALLOCATE(a(km))
      IF(KM .GT. sizea)then
         a(1:sizea) = dummy(1:sizea)
         a(sizea+1:km) =0.0
      ELSE
         a(1:km) = dummy(1:km)
      ENDIF
      DEALLOCATE(dummy)
      RETURN
      END




    SUBROUTINE check_beam_attrib(k,j,match)
! --------------------------------------------------------
! check attributes of beam j against those of beam k
! set match =1 IF perfect match
! set match =-1   IF ONLY difference is in
! tangency radius of source for the two beams 
! otherwise RETURN match = 0
! ------------------------------------------------------------
    USE transp,     ONLY : beam_data
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k,j
    REAL *8 ::  dtmatch = 0.001
    INTEGER match

    match = 1                                         !assume perfect match

    IF(ABS(beam_data%rtcena(j) -  beam_data%rtcena(k)) .GT. 0.01) match = -1


    IF(ABS(beam_data%tbona(j) - beam_data%tbona(k)) .GT. dtmatch) match =0
    IF(ABS(beam_data%tboffa(j) - beam_data%tboffa(k)) .GT. dtmatch) match =0


    IF(ABS(beam_data%abeama(j) -  beam_data%abeama(k)) .GT. 0.01)match =0
    IF(ABS(beam_data%xzbeama(j) -  beam_data%xzbeama(k)) .GT. 0.01)match =0

    IF(ABS(beam_data%aplasm(j) -  beam_data%aplasm(k)) .GT. 0.01)match =0
    IF(ABS(beam_data%backz(j) -  beam_data%backz(k)) .GT. 0.01)match =0

    IF(ABS(beam_data%bmwidra(j) -  beam_data%bmwidra(k)) .GT. 0.01)match =0
    IF(ABS(beam_data%bmwidza(j) -  beam_data%bmwidza(k)) .GT. 0.01)match =0

    IF(ABS(beam_data%xlbtna(j) -  beam_data%xlbtna(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%xybsca(j) -  beam_data%xybsca(k)) .GT. 0.01)match = 0

    IF(ABS(beam_data%divza(j) -  beam_data%divza(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%divra(j) -  beam_data%divra(k)) .GT. 0.01)match = 0


    IF(ABS(beam_data%foclza(j) -  beam_data%foclza(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%foclra(j) -  beam_data%foclra(k)) .GT. 0.01)match = 0


    IF(ABS(beam_data%rapedga(j) -  beam_data%rapedga(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%xzpedga(j) -  beam_data%xzpedga(k)) .GT. 0.01)match = 0


    IF(ABS(beam_data%xlbapa2(j) -  beam_data%xlbapa2(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%xlbapa(j) -  beam_data%xlbapa(k)) .GT. 0.01)match = 0


    IF(ABS(beam_data%xybapa(j) -  beam_data%xybapa(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%rapedg2(j) -  beam_data%rapedg2(k)) .GT. 0.01)match = 0

    IF(ABS(beam_data%xzpedg2(j) -  beam_data%xzpedg2(k)) .GT. 0.01)match = 0

    IF(ABS(beam_data%xzimpx(j) -  beam_data%xzimpx(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%aimpx(j) -  beam_data%aimpx(k)) .GT. 0.01)match = 0

    IF(ABS(beam_data%xbzeta(j) -  beam_data%xbzeta(k)) .GT. 0.01)match = 0

    IF(ABS(beam_data%pinja(j) -  beam_data%pinja(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%einja(j) -  beam_data%einja(k)) .GT. 0.01)match = 0

    IF(ABS(beam_data%ffulla(j) -  beam_data%ffulla(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%fhalfa(j) -  beam_data%fhalfa(k)) .GT. 0.01)match = 0


    IF(ABS(beam_data%nbshapa(j) -  beam_data%nbshapa(k)) .GT. 0.01)match = 0
    IF(ABS(beam_data%nbapsha(j) -  beam_data%nbapsha(k)) .GT. 0.01)match = 0


!    IF(beam_data%nlco(j) .NE.  beam_data%nlco(k))match = 0
    IF(beam_data%nlco(j))THEN
       IF(.NOT. beam_data%nlco(k))match =0
    ELSEIF( .NOT. beam_data%nlco(j))THEN
       IF(beam_data%nlco(k))match =0
    ENDIF
    IF(ABS(beam_data%nbapsha(j) -  beam_data%nbapsha(k)) .GT. 0.01)match = 0

    END SUBROUTINE check_beam_attrib





    SUBROUTINE    set_nubeam_impure
! -----------------------------------------------------------------
! define nrhix, the number of impurity species
! and associated charge,mass no
! xzimpx,aimpx
! ------------------------------------------------------------------
      USE numbrs, ONLY:  nimp
      USE transp, ONLY : beam_data
      USE ions ,ONLY : namei
      INTEGER nrhix,j
      

       IF(beam_data%nrhix .GE. nimp)THEN
             beam_data%nrhix = nimp
       ELSE
          PRINT *,'beam_data%nrhix =',beam_data%nrhix
          PRINT *,'nrhix =',beam_data%nrhix
          CALL STOP('nrhix problem',1)
       ENDIF

!impurity ions:


        DO j=1,beam_data%nrhix 
           SELECT CASE(TRIM(namei(j)))   !Name   of ith impurity ion species  
              CASE( 'he')                ! helium     
                   beam_data%xzimpx(j) = 2.
                   beam_data%aimpx(j)  = 4.002
              CASE(  'c')                !carbon                         
                   beam_data%xzimpx(j) = 6.
                   beam_data%aimpx(j)  = 12.011
                   beam_data%aimpx(j)  = 12.000
              CASE('o' )                 ! oxygen 
                   beam_data%xzimpx(j) = 8.
                   beam_data%aimpx(j)  = 15.9994
              CASE( 'si')                !silicon 
                   beam_data%xzimpx(j) = 14.            
                   beam_data%aimpx(j)  = 28.0855       
              CASE('ar')                 ! argon 
                   beam_data%xzimpx(j) = 18.
                   beam_data%aimpx(j)  = 39.948
              CASE( 'ti')                ! titanium                      
                   beam_data%xzimpx(j) = 22.            
                   beam_data%aimpx(j)  = 47.88
              CASE('cr')                 !  chromium                    
                   beam_data%xzimpx(j) = 24.
                   beam_data%aimpx(j)  = 51.9961
              CASE('fe')                 ! iron 
                   beam_data%xzimpx(j) = 26.    
                   beam_data%aimpx(j)  = 55.847
              CASE('ni')                 ! nickel 
                   beam_data%xzimpx(j) = 28.              
                   beam_data%aimpx(j)  = 58.6934          
              CASE('kr')                 ! krypton 
                   beam_data%xzimpx(j) = 36.       
                   beam_data%aimpx(j)  = 83.80
              CASE('mo')                 !molybdenum  
                   beam_data%xzimpx(j) = 42.    
                   beam_data%aimpx(j)  = 95.94
              CASE('w')                  ! tungsten   
                   beam_data%xzimpx(j) = 74.
                   beam_data%aimpx(j)  = 183.84
              CASE DEFAULT
                   CALL STOP('sub init,atom .no not known',0)
           END SELECT
        ENDDO

    END  SUBROUTINE set_nubeam_impure



    SUBROUTINE set_nubeam_thermal
! -----------------------------------------------------------------
! define ngmax, the number of thermal primary ion species and
! associated charge,mass no
! backz and aplasm
! ------------------------------------------------------------------
      USE numbrs, ONLY:  nprim
      USE transp, ONLY : beam_data
      USE ions ,ONLY : namep,atw,atomno
      IMPLICIT NONE
      INTEGER ngmax,j

      ngmax =0
      DO j=1,nprim
         IF (namep(j) .NE. 'dt')  THEN
            ngmax = ngmax+1
            beam_data%backz(ngmax)  = atomno(j)
            beam_data%aplasm(ngmax) = atw(j)
 
         ELSE
            ngmax = ngmax+1                      !deuterium
            beam_data%backz(ngmax) =1.
            beam_data%aplasm(ngmax) = 2.0
            ngmax = ngmax+1                      !tritium
            beam_data%backz(ngmax) =1.
            beam_data%aplasm(ngmax) = 3.0       
         ENDIF
       ENDDO
       IF(beam_data%ngmax .GE. ngmax)THEN
             beam_data%ngmax = ngmax
       ELSE
          PRINT *,'beam_data%ngmax =',beam_data%ngmax
          PRINT *,'ngmax =',ngmax
          CALL STOP('ngmax problem',1)
       ENDIF


      RETURN
    END  SUBROUTINE set_nubeam_thermal
