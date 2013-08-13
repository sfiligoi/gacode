      MODULE ions_gcnmp   
      USE nrtype ,   ONLY : DP,I4B
      IMPLICIT NONE
!    name_size = 8 changed to 10 5/2/11 to accomadte onetwo P_nfreya namelist writes
! change it back to 8(caused statefile problems)
!    name_size = 8 must match size in ions.f90  !!!
      INTEGER,PARAMETER,PUBLIC :: name_size = 8,nimp_ml= 11,nprim_ml = 6
!      INTEGER,PARAMETER,PUBLIC :: name_size = 10,nimp_ml= 11,nprim_ml = 6
      INTEGER(I4B),PUBLIC :: nion,nprim,nprimp1,nneu,nimp
      INTEGER(I4B),PUBLIC :: ntrp,ntri  ! number of primary,impurity ions run in simulation
      CHARACTER(len = name_size) ,PUBLIC, DIMENSION(:), POINTER :: namep
      CHARACTER(len = name_size) ,PUBLIC, DIMENSION(:), POINTER :: namei
      CHARACTER(len = name_size) ,PUBLIC, DIMENSION(nimp_ml)    :: namei_ml ! set in set_ion_prop
      CHARACTER(len = name_size) ,PUBLIC, DIMENSION(nprim_ml)   :: namep_ml ! set in set_ion_prop
!     namei_ml is master list of impurities. If this list is changed then corresponding
!     items in get_charge_state must bemodified as well
      CHARACTER(len = name_size) ,PUBLIC, DIMENSION(:), POINTER :: namen
!      CHARACTER(len = name_size) ,PUBLIC, DIMENSION(:), POINTER :: nameb
      CHARACTER(len = name_size) ,PUBLIC, DIMENSION(:), POINTER :: namefus ! currentlly not set
!

      REAL(DP),PUBLIC, DIMENSION(nimp_ml)        :: atwi_ml     
      REAL(DP),PUBLIC, DIMENSION(nprim_ml)       :: atwp_ml    
      REAL(DP),PUBLIC, DIMENSION(nimp_ml)        :: atomnoi_ml
      REAL(DP),PUBLIC, DIMENSION(nprim_ml)       :: atomnop_ml

      REAL(DP),PUBLIC, DIMENSION(:),   POINTER   :: atw     !(nion)
      REAL(DP),PUBLIC, DIMENSION(:),   POINTER   :: atomno  !(nion)
      REAL(DP),PUBLIC, DIMENSION(:),   POINTER   :: dmassden 
      REAL(DP),PUBLIC, DIMENSION(:,:), POINTER   :: dzdtim  !(nj, kimp)
      REAL(DP),PUBLIC, DIMENSION(:,:), POINTER   :: z,z_back!(nj,kion)
      REAL(DP),PUBLIC, DIMENSION(:,:), POINTER   :: zsq     !(nj,kion)
      REAL(DP),PUBLIC, DIMENSION(:,:), POINTER   :: dzdte   !(nj, kion)
      REAL(DP),PUBLIC, DIMENSION(:),   POINTER   :: zeff    !(nj)
      REAL(DP),PUBLIC, DIMENSION(:),   POINTER   :: tot_primary_ion_den    !(nj)
      REAL(DP),PUBLIC, DIMENSION(:),   POINTER   :: tot_thermal_ion_den    !(nj)
      REAL(DP),PUBLIC, DIMENSION(:),   POINTER   :: ni_sc
      REAL(DP),PUBLIC                            :: rfatmwt,fd_thermal
      INTEGER(I4B),PUBLIC, DIMENSION(:),POINTER  :: fi_index ! index of fast ion species in namep

      INTEGER(I4B),PUBLIC, DIMENSION(:),POINTER  :: beam_index   ! index of nameb(i) in nameb_ml
      INTEGER(I4B),PUBLIC, DIMENSION(:),POINTER  :: nprim_index  ! index of namep(i) in namep_ml
      INTEGER(I4B),PUBLIC, DIMENSION(:),POINTER  :: nimp_index   ! index of namei(i) in namei_ml

      INTEGER(I4B),PUBLIC                        :: d_index  ! index of deuterium in namep
      INTEGER(I4B),PUBLIC                        :: t_index  !          tritium
      INTEGER(I4B),PUBLIC                        :: he_index !          he 
      INTEGER(I4B),PUBLIC                        :: dt_index !          fd_thermal dt mixture
      INTEGER(I4B),PUBLIC                        :: h_index 

      DATA namei_ml /'ar',  'c ',  'cr',  'fe',  'he',  'kr',  'mo',  'ni',   &
                     'o ',  'si',  'w ' /
      DATA atwi_ml / 40._DP,12._DP,52._DP,56._Dp,4._DP,84._DP,96._DP,59._Dp, &
                     16._DP, 28._DP,184._DP /
      DATA atomnoi_ml/18._DP,6._DP,24._Dp,26._DP,2._DP,36._DP,42._DP,28._DP,  &
                      8._DP,14._DP,74._DP /


      DATA namep_ml / 'h','d','t','dt','he3','he4' / ! this list is currently only used in nfreya
                                                     ! the order of items in the lists 
                                                     ! must be maintained
      DATA atwp_ml /   1._DP,2._DP,3._DP,2.5_DP,3._Dp,4._DP/ ! note fudge for dt assumes 50/50
      DATA atomnop_ml/ 1._DP,1._DP,1._DP,1._DP,2._Dp,2._DP/



      END MODULE ions_gcnmp
