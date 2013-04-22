!-----------------------------------------------------------------------
! file hdf_psapi
! hdf_api module
!  Very generic module meant for writing HDF5 files with particular
!  attributes.  In case we want to convert over to C, it should make
!  this easier.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! code organization for hdf_psapi.
!-----------------------------------------------------------------------
! 4.  vshdf5_initprovenance(vsprov)
! 1 . write_provenance(gInId,grName,prin,h5err)
!-----------------------------------------------------------------------
! module hdf_api
!-----------------------------------------------------------------------
  module hdf5_psapi
  use hdf5_api
!-----------------------------------------------------------------------
! Provenance data for the codes.  See:
!  https://ice.txcorp.com/trac/vizschema/wiki/ProvenanceMetaData
! Some of these are hard to do in fortran
!-----------------------------------------------------------------------
  type provenancedata
   !att vstype = "runinfo"
   character(len=30) :: software        ! code name
   character(len=30) :: swVersion       ! software version
   character(len=30) :: swRevision      ! revision number from subversion
   character(len=30) :: vsVersion       ! vs compliance version
   character(len=30) :: fcCompiler      ! fortran compiler
   character(len=30) :: fcCompilerVersion      ! fc version
   character(len=30) :: fcCompilerFlags        ! fc flags
   character(len=30) :: buildHost              ! build host
   character(len=30) :: buildHostType          ! host type
   character(len=30) :: buildDate              ! Build date
   character(len=30) :: runHost                ! run host
   character(len=30) :: runHostType            ! run host type
   character(len=30) :: user                   ! user
   character(len=30) :: runDate                ! run date
   character(len=30) :: commandLine            ! commandline
  end type
  contains
!-----------------------------------------------------------------------
! subprogram 0.1 vshdf5_initprovenance
! Initialize these variables to default values.
!-----------------------------------------------------------------------
  subroutine vshdf5_initprovenance(vsprov)
  type(provenancedata), intent(inout) :: vsprov
  
  vsprov%software=" "    ! code name
  vsprov%swVersion=" "   ! software version
  vsprov%swRevision=" "  ! revision number from subversion
  vsprov%vsVersion="3.0.0"   ! VS compliance version
  vsprov%fcCompiler=" "  ! fortran compiler
  vsprov%fcCompilerVersion=" "  ! FC version
  vsprov%fcCompilerFlags=" "    ! FC flags
  vsprov%buildHost=" "          ! Build host
  vsprov%buildHostType=" "      ! Host type
  vsprov%buildDate=" "          ! Build host
  vsprov%runHost=" "            ! Date of build
  vsprov%runHostType=" "        ! Run Host type
  vsprov%runDate=" "            ! run date
  vsprov%user=" "               ! user
  vsprov%commandLine=" "        ! commandline
  return
  end subroutine vshdf5_initprovenance
!-----------------------------------------------------------------------
! subprogram 20. write_provenance
! Make a group that contains the provenance data.  See:
!  https://ice.txcorp.com/trac/vizschema/wiki/ProvenanceMetaData
!-----------------------------------------------------------------------
  subroutine write_provenance(gInId,grName,prin,h5in,h5err)
  integer(HID_T), intent(in) :: gInId
  character*(*), intent(in) :: grName
  type(provenanceData), intent(in) :: prin
  type(hdf5InOpts), intent(in) :: h5in
  type(hdf5ErrorType), intent(inoUT) :: h5err
  return
  end subroutine write_provenance

  end module hdf5_psapi

