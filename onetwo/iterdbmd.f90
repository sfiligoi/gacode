
   MODULE iterdbmd
   USE nrtype, ONLY :   I4B ,DP
   IMPLICIT NONE
!
! --- INCLUDE file for ITER database-related quantities
!
!     iterdb : inone (namelist 3) input param
!     niterdb: Fortran I/O unit number for iterdbfilename
!     iterdsc: inone (namelist 3) input parm
!     irwflag: inone (namelist 3) input parm
!     iterdbfilename: Fortran variable stores name of ITER database file
!     (actual name is assigned at start of cray101.f)
! --- other quantities related to ITER database are in file psig.i
!
      INTEGER(I4B) iterdb, niterdb, iterdsc, irwflag ,statefile_type
      character*128       iterdbfilename,statefile_name,iterdb_outpt
      LOGICAL  create_GCNMP_input,initialize_from_statefile,create_XPTOR_input
      REAL(DP) smth_mhd_parms !jmp.ibm
!
   END MODULE iterdbmd
