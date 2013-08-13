  MODULE gcnmp_version
    USE nrtype,                    ONLY : DP, I4B
    IMPLICIT NONE
    CHARACTER(LEN= 24) :: gcnmp_ver = 'GCNMP version 0.99'
    REAL(DP)  :: gcnmp_ver_no = 0.99
    !DATA gcnmp_ver /'GCNMP version 0.98' /
    ! v98 includes tglf from git hub (tglf v1.93?)
    ! v99 includes tglf from git hub (taken on Nov 2012) (tglf v1.94?) 
    ! thisv98,99 version have parallel  grid/wave vector implementation
! NOTE: python plot programs  check for 3 fields in gcnmp_ver
 
  END MODULE gcnmp_version
