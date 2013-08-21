!
  MODULE beam_structure
     USE nrtype,                          ONLY : DP,SP,I4B

     IMPLICIT NONE 


     TYPE,PUBLIC ::  neutral_beam
        INTEGER                            :: nbeam,ngmax,nzone_nb,nfast
        INTEGER                            :: nrhix,nznbma,nznbme,nbbcal
        INTEGER                            :: nseed,nbeam_species
        INTEGER                            :: nfusion_species
        LOGICAL ::            data_allocated
        CHARACTER*40, DIMENSION(:),POINTER :: label_thions   ! thermal plasma
        CHARACTER*40, DIMENSION(:),POINTER :: label_fastions
        CHARACTER*40, DIMENSION(:),POINTER :: label_impions
        CHARACTER*8,  DIMENSION(:),POINTER :: label_beamions
        CHARACTER*32, DIMENSION(:),POINTER :: beam_id
        INTEGER,  DIMENSION(:),POINTER  :: ntrace
        INTEGER,  DIMENSION(:),POINTER  :: nbapsh2
        INTEGER,  DIMENSION(:),POINTER  :: nbshapa
        INTEGER,  DIMENSION(:),POINTER  :: nbapsha
        INTEGER,  DIMENSION(:),POINTER  :: nznbmea
        INTEGER,  DIMENSION(:),POINTER  :: nlfprod
        INTEGER,  DIMENSION(:),POINTER :: istrtbm
        INTEGER,  DIMENSION(:),POINTER :: iendbm
        INTEGER,  DIMENSION(:),POINTER :: beamlet_active
        INTEGER,  DIMENSION(:),POINTER :: beamlet_cp
        LOGICAL,  DIMENSION(:),POINTER  :: nlco
        LOGICAL                         :: nlbdat
        LOGICAL                         :: nlbdat_set
        REAL(DP),  DIMENSION(:),POINTER :: backz
        REAL(DP),  DIMENSION(:),POINTER :: aplasm
        REAL(DP),  DIMENSION(:),POINTER :: xzbeama
        REAL(DP),  DIMENSION(:),POINTER :: abeama
        REAL(DP),  DIMENSION(:),POINTER :: rtcena
        REAL(DP),  DIMENSION(:),POINTER :: xlbtna
        REAL(DP),  DIMENSION(:),POINTER :: xybsca
        REAL(DP),  DIMENSION(:),POINTER :: bmwidra
        REAL(DP),  DIMENSION(:),POINTER :: bmwidza
        REAL(DP),  DIMENSION(:),POINTER :: divza
        REAL(DP),  DIMENSION(:),POINTER :: divra
        REAL(DP),  DIMENSION(:),POINTER :: foclza
        REAL(DP),  DIMENSION(:),POINTER :: foclra
        REAL(DP),  DIMENSION(:),POINTER :: rapedga
        REAL(DP),  DIMENSION(:),POINTER :: xzpedga
        REAL(DP),  DIMENSION(:),POINTER :: xlbapa2
        REAL(DP),  DIMENSION(:),POINTER :: xlbapa
        REAL(DP),  DIMENSION(:),POINTER :: xybapa
        REAL(DP),  DIMENSION(:),POINTER :: rapedg2
        REAL(DP),  DIMENSION(:),POINTER :: xzpedg2
        REAL(DP),  DIMENSION(:),POINTER :: xzimpx
        REAL(DP),  DIMENSION(:),POINTER :: aimpx
        REAL(DP),  DIMENSION(:),POINTER :: xbzeta
        REAL(DP),  DIMENSION(:),POINTER :: pinja      
        REAL(DP),  DIMENSION(:),POINTER :: einja  
        REAL(DP),  DIMENSION(:),POINTER :: ffulla
        REAL(DP),  DIMENSION(:),POINTER :: fhalfa
        REAL(DP),  DIMENSION(:),POINTER :: abeams
        REAL(DP),  DIMENSION(:),POINTER :: xzbeams
        REAL(DP),  DIMENSION(:),POINTER :: ebdmaxa
        REAL(DP),  DIMENSION(:),POINTER :: xrapoffa
        REAL(DP),  DIMENSION(:),POINTER :: xzapoffa
        REAL(DP),  DIMENSION(:),POINTER :: xrapoff2
        REAL(DP),  DIMENSION(:),POINTER :: xzapoff2
!
        REAL(DP),  DIMENSION(:,:),POINTER :: efbm,efbmb
        REAL(DP),  DIMENSION(:),POINTER   :: beam_times  ! ufile
        REAL(DP),  DIMENSION(:),POINTER   :: beam_times_prev

#ifdef NFREYA    ! for P_Nfreya use DP
        REAL(DP),  DIMENSION(:,:),POINTER :: beam_inject ! ufile
        REAL(DP),  DIMENSION(:),POINTER   :: beam_chan   ! ufile
#endif

#ifdef ONETWO    ! for Onetwo use SP
        REAL(SP),  DIMENSION(:,:),POINTER :: beam_inject ! ufile
        REAL(SP),  DIMENSION(:),POINTER   :: beam_chan   ! ufile
#endif
        REAL(DP),  DIMENSION(:),POINTER   :: tbona,tboffa,tbonac,tboffac
        REAL(DP),  DIMENSION(:),POINTER   :: switch_on_time,switch_off_time
        REAL(DP),  DIMENSION(:),POINTER   :: switch_on_off_delay
        REAL(DP),  DIMENSION(:),POINTER   :: beam_switch_times
        REAL(DP),  DIMENSION(:),POINTER   :: beam_switch_times_prev
        REAL(DP),  DIMENSION(:),POINTER   :: beam_sw
        REAL(DP)   pwf_tot_intg 
        REAL(SP)   ebdmax  ! both Onetwo and P_Nreya
        REAL(DP)   beam_power_rise_time 
        REAL(DP)   beam_power_decay_time 
    END TYPE neutral_beam
!
!
    CONTAINS
!

    SUBROUTINE beam_structure_allocate(beam_struct,nbeam_tr,nrhix,ngmax,nfast)
! ------------------------------------------------------------------------
!  allocate space for the data structure name_beam
!  INPUT 
!     nbeam_tr # beams
!     nrhix    # impurity species
!     ngmax    # total species
!     nfast    # fast species
! --------------------------------------------------------HSJ mod 3/23/11-
    USE allocate_err,                         ONLY : allocate_error

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbeam_tr,nrhix,ngmax,nfast
    INTEGER istat
    TYPE(neutral_beam),INTENT(inout):: beam_struct
!
!
         beam_struct%data_allocated = .TRUE.
         beam_struct%nbeam = nbeam_tr
         beam_struct%nfast = nfast
         beam_struct%nrhix = nrhix
         beam_struct%ngmax = ngmax
         beam_struct%beam_power_rise_time = 5.e-4
         beam_struct%beam_power_decay_time = beam_struct%beam_power_rise_time
         ALLOCATE (beam_struct%ntrace(1:nbeam_tr),STAT = istat)
           IF(istat .NE. 0)                                              &
             CALL allocate_error("ntrace,sub beam_structure_allocate",0,istat)
!
         ALLOCATE (beam_struct%nbapsh2(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("nbapsh2,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%nbshapa(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
            CALL allocate_error("nbshapa,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%nbapsha(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("nbapsha,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%nlco(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("nlco,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%backz(1:ngmax),STAT = istat)
         IF(istat .NE. 0)                                                &
             CALL allocate_error("backz,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%aplasm(1:ngmax),STAT = istat)
         IF(istat .NE. 0)                                                &
              CALL allocate_error("aplasm,sub beam_structure_allocate",0,istat)

!
         ALLOCATE (beam_struct%xzbeama(1:nbeam_tr),STAT = istat)
           IF(istat .NE. 0)                                              &
             CALL allocate_error("xzbeama,sub beam_structure_allocate",0,istat)
!
         ALLOCATE (beam_struct%abeama(1:nbeam_tr),STAT = istat)
           IF(istat .NE. 0)                                              &
             CALL allocate_error("abeama,sub beam_structure_allocate",0,istat)
!
         ALLOCATE (beam_struct%rtcena(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("rtcena,sub beam_structure_allocate",0,istat)

         ALLOCATE (beam_struct%xrapoffa(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xrapoffa,sub beam_structure_allocate",0,istat)

         ALLOCATE (beam_struct%xzapoffa(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xzapoffa,sub beam_structure_allocate",0,istat)



         ALLOCATE (beam_struct%xrapoff2(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xrapoff2,sub beam_structure_allocate",0,istat)

         ALLOCATE (beam_struct%xzapoff2(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xzapoff2,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%xlbtna(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xlbtna,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%xybsca(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xybsca,sub beam_structure_allocate",0,istat)
!
 
             ! elevation of center of beam source grid above vacuum vessel
             ! midplane "Z=0". (negative if located below midplane). cm.
             ! default value, 0.0, means source lies on the midplane.
             ! in an axisymmetric system, rtcena, xlbtna, xybsca suffice 
             ! to specify the source grid position.
         ALLOCATE (beam_struct%bmwidra(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("bmwidra,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%bmwidza(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("bmwidza,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%divra(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("divra,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%divza(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("divza,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%foclra(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("foclra,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%foclza(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("foclza,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%rapedga(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("rapedga,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%xzpedga(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xzpedga,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%xlbapa2(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xlbapa2,sub beam_structure_allocate",0,istat)
!
         ALLOCATE (beam_struct%xlbapa(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xlbapa,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%xzpedg2(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xzpedg2,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%rapedg2(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("rapedg2,sub beam_structure_allocate",0,istat)
!
         ALLOCATE (beam_struct%xybapa(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xybapa,sub beam_structure_allocate",0,istat)
!
             ! aperture half-width (or radius if circular), cm
             ! value of zero means: no 2nd aperture
!
         ALLOCATE (beam_struct%xzpedg2(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xzpedg2,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%xzimpx(1:nrhix),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xzimpx,sub beam_structure_allocate",0,istat)
!
         ALLOCATE (beam_struct%aimpx(1:nrhix),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("aimpx,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%xbzeta(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("aimpx,sub beam_structure_allocate",0,istat)
!
         ALLOCATE (beam_struct%pinja(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("pinja,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%ffulla(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("ffulla,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%fhalfa(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("fhalfa,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%einja(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("einja,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%tbona(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("tbona,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%tboffa(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("tboffa,sub beam_structure_allocate",0,istat)
!
!
!

        ALLOCATE (beam_struct%istrtbm(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("istrtbm,sub init",0,istat)
             beam_struct%istrtbm(:) = 0

         ALLOCATE (beam_struct%iendbm(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("iendbm,sub init",0,istat)
             beam_struct%iendbm(:) = 0

         ALLOCATE (beam_struct%beamlet_active(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("beamlet_active,sub init",0,istat)
             beam_struct%beamlet_active(:) = - 100
        ALLOCATE (beam_struct%beamlet_cp(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("beamlet_cp,sub init",0,istat)
             beam_struct%beamlet_cp(:) = beam_struct%beamlet_active(:)

!
!
!
         ALLOCATE (beam_struct%tbonac(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("tbonac,sub init",0,istat)
!
!
         ALLOCATE (beam_struct%tboffac(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("tboffa,sub init",0,istat)




         ALLOCATE (beam_struct%switch_on_time(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("switch_on_time,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%switch_off_time(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("switch_off_time,sub beam_structure_allocate",0,istat)

         ALLOCATE (beam_struct%beam_id(1:nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("beam_id ,sub beam_structure_allocate",0,istat)
!
!
!
!
! fast ion quantities (not necessarily related directly to beam)
!
!
!
         ALLOCATE (beam_struct%abeams(1:nfast),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("abeams,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%xzbeams(1:nfast),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("xzbeams,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%ebdmaxa(1:nfast),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("ebdmaxa,sub beam_structure_allocate",0,istat)
!
!
!
         ALLOCATE (beam_struct%nznbmea(1:nfast),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("nznbmea,sub beam_structure_allocate",0,istat)
!
!
         ALLOCATE (beam_struct%nlfprod(1:nfast),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("nlfprod,sub beam_structure_allocate",0,istat)
!
!
!
!  NOTE: beam_struct%beam_times and beam_struct%beam_inject are not
!        allocated here because we dont know the size at this point.
         ALLOCATE (beam_struct%beam_chan(4*nbeam_tr),STAT = istat)
        IF(istat .NE. 0)                                                 &
             CALL allocate_error("beam_chan,sub beam_structure_allocate",0,istat)
    RETURN
    END SUBROUTINE beam_structure_allocate
!
!
  END MODULE beam_structure
