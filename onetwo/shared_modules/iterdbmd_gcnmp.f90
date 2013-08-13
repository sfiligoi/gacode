
   MODULE iterdbmd_gcnmp
   USE nrtype, ONLY :   I4B ,DP
   USE solcon_gcnmp, ONLY :   std_file_size
   IMPLICIT NONE
!
! --- INCLUDE file for ITER database-related quantities
!
!     iterdb : inone (namelist 3) input param
!     niterdb: Fortran I/O unit number for iterdb_file_name
!     iterdsc: inone (namelist 3) input parm
!     irwflag: inone (namelist 3) input parm
!     iterdb_file_name: Fortran variable stores name of ITER database file
!     (actual name is assigned at start of cray101.f)
! --- other quantities related to ITER database are in file psig.i
!
      INTEGER(I4B) iterdsc, irwflag,ncid
      CHARACTER(LEN = std_file_size)       iterdb_file_name,iterdb_outpt
      CHARACTER(len =3) idb_sufix
      
!
   END MODULE iterdbmd_gcnmp
