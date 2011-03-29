!-----------------------------------------------------------------------
!     file hdf_api
!     hdf_api module
!      Very generic module meant for writing HDF5 files with particular
!      attributes.  In case we want to convert over to C, it should make
!      this easier.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization for hdf_api.
!-----------------------------------------------------------------------
!     0.  check_dims
!     1.  vshdf5_fcinit
!     2.  open_oldh5file
!     3.  open_newh5file
!     4.  close_h5file
!     5.  open_group
!     6.  make_group
!     7.  close_group
!     8.  test_group
!     9.  get_nmembers
!     10. make_reference
!     11. write_attribute_ch_sc
!     12. write_attribute_ch_vec
!     13. write_attribute_int_sc
!     14. write_attribute_int_vec
!     15. write_attribute_rl_sc
!     16. write_attribute_rl_vec
!     17. dump_h5in_attributes
!     18. dump_h5_1d
!     19. dump_h5_2d
!     20. dump_h5_3d
!     20. dump_h5_4d
!     21. read_dims
!     22. read_h5_1d
!     23. read_h5_2d
!     24. read_h5_3d
!-----------------------------------------------------------------------
!     module hdf_api
!-----------------------------------------------------------------------
      MODULE hdf5_api
      USE HDF5
      IMPLICIT NONE
      CHARACTER(5), PARAMETER, PRIVATE :: h5fortranApiVersion="1.0"
      INTEGER, PARAMETER, PRIVATE :: i4=SELECTED_INT_KIND(9)
      INTEGER, PARAMETER, PRIVATE :: i8=SELECTED_INT_KIND(18)
      INTEGER, PARAMETER, PRIVATE :: r4=SELECTED_REAL_KIND(6,37)
      INTEGER, PARAMETER, PRIVATE :: r8=SELECTED_REAL_KIND(13,307)
!-----------------------------------------------------------------------
!     Input parameters to control the attributes and how written out
!     The method used to determine whether to write the character
!     variables are written out is something like the following:
!        IF(TRIM(h5in%vsMD)/="") THEN
!           CALL write_attribute(dset_id,"vsMD",h5in%vsMD,errval)
!        ENDIF
!     This needs to be paid attention to if writing out several
!     variables in a row.
!
!     MESH DISCUSSION:
!     The key to making the vsSchema work is to associate fields with
!     their meshes.  The way to point the mesh it use the mesh member
!     the derived type below.  It is also used to define the mesh by
!     prepending the mesh type from visSchemawith "mesh-"; e.g.,
!      h5in%mesh="mesh-structured" defines a mesh and
!      h5in%mesh="/coreMesh"       points to where the mesh is defined.
!     For valid types of meshes, details on how the multi-domain
!     specification works, and centering issues, see the visSchema wiki:
!      https://ice.txcorp.com/trac/vizschema/wiki/
!
!     WRD_TYPE:  Valid wrd_types are:
!      H5T_NATIVE_REAL
!      H5T_NATIVE_DOUBLE
!      H5T_NATIVE_INTEGER
!      H5T_NATIVE_INT16   !?
!-----------------------------------------------------------------------
!     IMPORTANT:::
!     It is very important to have codes use this API use the initvars
!     When adding variables to this derived type, make sure they are 
!      initializing correctly.
!-----------------------------------------------------------------------
      TYPE hdf5InOpts
         INTEGER(HID_T) :: wrd_type
         LOGICAL :: doTranspose        
         LOGICAL :: verbose               ! Whether to write verbose output
         LOGICAL :: debug                 ! Write even more verbose output for debugging
         LOGICAL :: pIO                   ! Whether parallel I/O is used
         LOGICAL :: wrVsTime              ! Whether to write time to attribute
         LOGICAL :: typeConvert           ! Whether to demote the type
         LOGICAL :: unitConvert           ! Whether to write vsUnitCnv
         INTEGER(SELECTED_INT_KIND(9)) :: comm     ! Communicator associated w/ open file
         INTEGER(SELECTED_INT_KIND(9)) :: info     ! mpi_info: set to MPI_INFO_NULL if unsure
         CHARACTER(LEN=30) :: mesh             ! See above
         CHARACTER(LEN=30) :: units            ! Units
         CHARACTER(LEN=30) :: vsCentering      ! How to center variables on mesh
         CHARACTER(LEN=30) :: vsMD             ! Multidomain variable
         REAL(r8) :: vsTime                    ! Bassi requires r8
         INTEGER(HID_T) :: vsStep=0            ! Step # associated with time
         REAL(r8) :: vsUnitCnv=1.              ! conversion factor
         ! DOUBLE PRECISION :: vsTime          ! Different b/c time dependent
      END TYPE
!-----------------------------------------------------------------------
!     Example of how to use h5err type after an fcapi call.
!     IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
!-----------------------------------------------------------------------
      TYPE hdf5ErrorType
         LOGICAL :: errBool
         CHARACTER(64) :: errorMsg
      END TYPE
!-----------------------------------------------------------------------
!     Provenance data for the codes.  See:
!      https://ice.txcorp.com/trac/vizschema/wiki/ProvenanceMetaData
!     Some of these are hard to do in fortran
!-----------------------------------------------------------------------
      TYPE provenanceData
       !Att vsType = "runInfo"
       CHARACTER(LEN=30) :: software        ! code name
       CHARACTER(LEN=30) :: swVersion       ! software version
       CHARACTER(LEN=30) :: swRevision      ! revision number from subversion
       CHARACTER(LEN=30) :: vsVersion       ! VS compliance version
       CHARACTER(LEN=30) :: fcCompiler      ! fortran compiler
       CHARACTER(LEN=30) :: fcCompilerVersion      ! FC version
       CHARACTER(LEN=30) :: fcCompilerFlags        ! FC flags
       CHARACTER(LEN=30) :: buildHost              ! Build host
       CHARACTER(LEN=30) :: buildHostType          ! Host type
       CHARACTER(LEN=30) :: runHost                ! Run host
       CHARACTER(LEN=30) :: runHostType            ! Run Host type
       CHARACTER(LEN=30) :: user                   ! user
       CHARACTER(LEN=30) :: runDate                ! run date
       CHARACTER(LEN=30) :: commandLine            ! commandline
      END TYPE
!-----------------------------------------------------------------------
!     subprogram name interfaces
!-----------------------------------------------------------------------
      INTERFACE write_attribute
        MODULE PROCEDURE write_attribute_ch_sc,write_attribute_ch_vec &
                       ,write_attribute_rl_sc,write_attribute_rl_vec &
                       ,write_attribute_int_sc,write_attribute_int_vec
      END INTERFACE
      INTERFACE dump_h5
        MODULE PROCEDURE  dump_h5_int,dump_h5_rldbl, &
                         dump_h5_1d,dump_h5_2d,dump_h5_3d,dump_h5_4d, &
                         dump_rl_1d,dump_rl_2d,dump_rl_3d,dump_rl_4d
      END INTERFACE
      ! Dump components into separate groups
      INTERFACE read_h5
        MODULE PROCEDURE read_h5_1d,read_h5_2d,read_h5_3d
      END INTERFACE
      INTERFACE read_attribute
        MODULE PROCEDURE read_attribute_rl_sc
      END INTERFACE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram -1. write mismatched dimension errors
!-----------------------------------------------------------------------
      SUBROUTINE check_dims(dims, fdims, errval)
      INTEGER(HSIZE_T), DIMENSION(:), INTENT(IN) :: dims, fdims
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
         INTEGER i
         DO i = 1, SIZE(dims)
           IF (dims(i) /= fdims(i)) THEN
             write(*, *) "ERROR: dims (", dims, ")"
             write(*, *) "  /=  fdims (", fdims, ")"
             write(*, *) "fdims = dims in the file (use h5ls)"
             write(*, *) "dims  = dims allocated for array to read"
             errval%errorMsg = 'ERROR: dims /= fdims'
             errval%errBool = .true.
             RETURN
           ENDIF
         ENDDO
         errval%errBool = .false.
      END SUBROUTINE check_dims

!-----------------------------------------------------------------------
!     function 1. h5accessMethod
!     Return something the correct hdf5 access method given a more
!     memoral names
!-----------------------------------------------------------------------
      FUNCTION h5accessMethod(access_method)
      INTEGER(HID_T) :: h5accessMethod
      CHARACTER(*), INTENT(IN) :: access_method
      SELECT CASE(access_method)
      CASE("overwr")
         h5accessMethod=H5F_ACC_TRUNC_F                   ! Overwrite file
      CASE("rdwr")
         h5accessMethod=H5F_ACC_RDWR_F                    ! Read-write
      CASE("rdonly")
         h5accessMethod=H5F_ACC_RDONLY_F                  ! Read only
      END SELECT
      RETURN
      END FUNCTION h5accessMethod

!-----------------------------------------------------------------------
!     subprogram 0. vshdf5_fcinit
!     Open fortran hdf5 and set open/close parameters.
!-----------------------------------------------------------------------
      SUBROUTINE vshdf5_fcinit()
      ! integer(hid_t) :: err
      integer :: err
      ! write(*, *) "vshdf5_fcinit: entered"
      CALL h5dont_atexit_f(err)
      ! write(*, *) "vshdf5_fcinit: h5dont_atexit_f returned."
      CALL h5open_f(err)
      ! write(*, *) "vshdf5_fcinit: h5open_f returned."
      ! write(*, *) "vshdf5_fcinit: leaving."
      RETURN
      END SUBROUTINE vshdf5_fcinit

!-----------------------------------------------------------------------
!     subprogram 0.1 vshdf5_inith5vars
!     Initialize these variables to default values.
!-----------------------------------------------------------------------
      SUBROUTINE vshdf5_inith5vars(h5in, h5err)
      TYPE(hdf5InOpts), INTENT(INOUT) :: h5in
      TYPE(hdf5ErrorType), INTENT(INOUT) :: h5err
!-----------------------------------------------------------------------
!     Defaults
!-----------------------------------------------------------------------
      ! According to xlf:
      ! (E) Null literal string is not permitted.  A single blank is assumed.
      ! So these should be single blanks to avoid warnings
      h5in%wrd_type=H5T_NATIVE_DOUBLE
      h5in%vsCentering=" "
      h5in%doTranspose=.false.
      h5in%verbose=.false.
      h5in%debug=.false.
      h5in%pIO=.false.
      h5in%wrVsTime=.false.
      h5in%typeConvert=.false.
      h5in%unitConvert=.false.
      h5in%mesh =  " "
      h5in%units =  " "
      h5in%vsCentering =  " "
      h5in%vsMD =  " "
      h5err%errBool = .false.
      h5err%errorMsg =  " "
      RETURN
      END SUBROUTINE vshdf5_inith5vars

!-----------------------------------------------------------------------
!     subprogram 0.1 vshdf5_initprovenance
!     Initialize these variables to default values.
!-----------------------------------------------------------------------
      SUBROUTINE vshdf5_initprovenance(vsprov)
      TYPE(provenanceData), INTENT(INOUT) :: vsprov
      
      vsprov%software=" "    ! code name
      vsprov%swVersion=" "   ! software version
      vsprov%swRevision=" "  ! revision number from subversion
      vsprov%vsVersion="3.0.0"   ! VS compliance version
      vsprov%fcCompiler=" "  ! fortran compiler
      vsprov%fcCompilerVersion=" "  ! FC version
      vsprov%fcCompilerFlags=" "    ! FC flags
      vsprov%buildHost=" "          ! Build host
      vsprov%buildHostType=" "      ! Host type
      vsprov%runHost=" "            ! Run host
      vsprov%runHostType=" "        ! Run Host type
      vsprov%user=" "               ! user
      vsprov%runDate=" "            ! run date
      vsprov%commandLine=" "        ! commandline
      RETURN
      END SUBROUTINE vshdf5_initprovenance
!-----------------------------------------------------------------------
!     subprogram 1. open_oldh5file
!     Open file for writing and write file attributes
!     Create the group for the independent variables at this stage
!-----------------------------------------------------------------------
      SUBROUTINE open_oldh5file(fname,fileId,fdesc,rootGid,h5in,h5err)
      CHARACTER(*), INTENT(IN) :: fname,fdesc
      INTEGER(HID_T), INTENT(OUT) :: fileId,rootGid
      TYPE(hdf5ErrorType), INTENT(INOUT) :: h5err
      TYPE(hdf5InOpts), INTENT(INOUT) :: h5in
      INTEGER(HID_T) :: access_mode
      INTEGER(HID_T) :: plist_id       ! Property list identifier
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: info
      INTEGER(HID_T) :: error

      LOGICAL :: file_exists
      if (h5in%verbose) then
        write(*, *) " open_oldh5file: entered."
      endif
!-----------------------------------------------------------------------
!     Create and open the file
!     The if statements seem to work better to give the expected
!      behavior for read-write for already open file.
!-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(fname),EXIST=file_exists)
      IF (.NOT. file_exists) THEN
         h5err%errorMsg = 'ERROR: file does not exist: '//fname
         if (h5in%verbose) then
           write(*, *) 'ERROR: file does not exist: ', fname
         endif
         h5err%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Setup file access property list with parallel I/O access.
!    jrc 1jul09: Need not have parallel I/O access when only one proc is
!    reading
!-----------------------------------------------------------------------
#ifdef __MPI
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, h5in%comm, h5in%info, error)
! jrc 12jul09: I think these should be set outside.  A parallel
! program may not want to use parallel I/O, e.g., if it is doing I/O
! on only one rank.  Scott K, what do you think?
      h5in%pIO=.true.
#endif

!-----------------------------------------------------------------------
!     Open file -- collectively if needed
!-----------------------------------------------------------------------
      if (h5in%verbose) then
        write(*, *) "open_oldh5file: calling h5fopen_f."
      endif
#ifdef __MPI
      if (h5in%pIO) then
        CALL h5fopen_f(fname,H5F_ACC_RDONLY_F,fileId,error, &
                      access_prp=plist_id)
      else
#else
      CALL h5fopen_f(fname, H5F_ACC_RDONLY_F, fileId, error)
#endif
#ifdef __MPI
      endif
#endif
      if (h5in%verbose) then
        write(*, *) "open_oldh5file: h5fopen_f returned."
      endif

!-----------------------------------------------------------------------
!     Grab the root group id which is created by default
!-----------------------------------------------------------------------
      CALL h5gopen_f(fileId,"/",rootGid,error)
      IF (error==FAIL) THEN
         h5err%errorMsg = 'ERROR: Error grabbing root ID: '//fname
         h5err%errBool = .true.
         RETURN
      ENDIF
      h5err%errBool = .false.
      RETURN
      END SUBROUTINE open_oldh5file

!-----------------------------------------------------------------------
!     subprogram 1. open_newh5file
!     Open file for writing and write file attributes
!     Create the group for the independent variables at this stage
!-----------------------------------------------------------------------
      SUBROUTINE open_newh5file(fname,fileId,fdesc,rootGid,h5in,h5err)
      CHARACTER(*), INTENT(IN) :: fname,fdesc
      INTEGER(HID_T), INTENT(OUT) :: fileId,rootGid
      TYPE(hdf5InOpts), INTENT(INOUT) :: h5in
      TYPE(hdf5ErrorType), INTENT(INOUT) :: h5err
      INTEGER(HID_T) :: access_mode
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER(HID_T) :: plist_id       ! Property list identifier
      INTEGER :: error, info

      LOGICAL :: file_exists
!-----------------------------------------------------------------------
!    Setup file access property list with parallel I/O access.
!-----------------------------------------------------------------------
#ifdef __MPI
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        CALL h5pset_fapl_mpio_f(plist_id, h5in%comm, h5in%info, error)
        h5in%pIO=.true.
#endif
!-----------------------------------------------------------------------
!     Open file -- collectively if needed
!     The if statements seem to work better to give the expected
!      behavior for read-write for already open file.
!-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(fname),EXIST=file_exists)
      IF (.NOT. file_exists) THEN
         ! Always create with over-write to avoid errors
         access_mode=h5accessMethod("overwr")
#ifdef __MPI
         CALL h5fcreate_f(TRIM(fname),access_mode,fileId,error, &
                      access_prp=plist_id)
#else
         CALL h5fcreate_f(TRIM(fname),access_mode,fileId,error)
#endif
      ELSE
         OPEN(UNIT=999,FILE=TRIM(fname),FORM='UNFORMATTED', &
           POSITION='REWIND',STATUS='REPLACE')
         CLOSE(UNIT=999)
         access_mode=h5accessMethod("overwr")
#ifdef __MPI
         CALL h5fcreate_f(TRIM(fname),access_mode,fileId,error, &
                      access_prp=plist_id)
#else
         CALL h5fcreate_f(TRIM(fname),access_mode,fileId,error)
#endif
      ENDIF
      IF (error==FAIL) THEN
         h5err%errorMsg = 'ERROR: Error opening file: '//fname
         h5err%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Grab the root group id which is created by default
!-----------------------------------------------------------------------
      CALL h5gopen_f(fileId,"/",rootGid,error)
      IF (error==FAIL) THEN
         h5err%errorMsg = 'ERROR: Error grabbing root ID: '//fname
         h5err%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Put in things which describe the api
!-----------------------------------------------------------------------
!      CALL write_attribute(rootGid,"vsFcVERSION","1.0",h5err)
!-----------------------------------------------------------------------
!     Put in file description when creating file
!-----------------------------------------------------------------------
      CALL write_attribute(rootGid,"Description",fdesc,h5err)
!-----------------------------------------------------------------------
      h5err%errBool = .false.
      RETURN
      END SUBROUTINE open_newh5file

!-----------------------------------------------------------------------
!     subprogram 2. close_h5file
!     Close the file associated with fileId.
!-----------------------------------------------------------------------
      SUBROUTINE close_h5file(fileId,root_id,h5err)
      INTEGER, INTENT(IN) :: fileId,root_id
      TYPE(hdf5ErrorType), INTENT(INOUT) :: h5err
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     Close the root group and file
!-----------------------------------------------------------------------
      CALL h5gclose_f(root_id, error)
      CALL h5fclose_f(fileId, error)
      IF (error==FAIL) THEN
         h5err%errorMsg = 'ERROR: Error in close_h5file'
         h5err%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      h5err%errBool = .false.
      RETURN
      END SUBROUTINE close_h5file

!-----------------------------------------------------------------------
!     subprogram 3. open_group
!     Open a group in a safe way
!-----------------------------------------------------------------------
      SUBROUTINE open_group(inid,gname,gid,errval)
      CHARACTER(*), INTENT(IN) :: gname
      INTEGER(HID_T), INTENT(IN) :: inid
      INTEGER(HID_T), INTENT(OUT) :: gid
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     Open group
!-----------------------------------------------------------------------
      CALL h5gopen_f(inid,gname,gid,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Error opening group: '//gname
         errval%errBool = .true.
         RETURN
      ENDIF
      errval%errBool = .false.
      RETURN
      END SUBROUTINE open_group

!-----------------------------------------------------------------------
!     subprogram 3. make_group
!     Create a group in a safe way
!-----------------------------------------------------------------------
      SUBROUTINE make_group(inid,gname,gid,meshtitle,errval)
      CHARACTER(*), INTENT(IN) :: gname
      CHARACTER*(*), INTENT(IN) :: meshtitle
      INTEGER(HID_T), INTENT(IN) :: inid
      INTEGER(HID_T), INTENT(OUT) :: gid
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error

      CHARACTER(64) :: msg
!-----------------------------------------------------------------------
!     Create group
!-----------------------------------------------------------------------
      !CALL h5gopen_f(inid,gname,gid,error)
      ! If it failed, most likely it doesn't exist
      CALL h5gcreate_f(inid,gname,gid,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Error opening group: '//gname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
!      CALL write_attribute(gid,'CLASS','GROUP',errval)
!      CALL write_attribute(gid,'VERSION','1.0',errval)
!       CALL write_attribute(gid,'TITLE',TRIM(meshtitle),errval)
      IF(LEN_TRIM(meshtitle)>0) THEN
       CALL write_attribute(gid,'vsType','mesh',errval)
       CALL write_attribute(gid,'vsKind','structured',errval)
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE make_group

!-----------------------------------------------------------------------
!     subprogram 4. close_group
!     Close a group in a safe way
!-----------------------------------------------------------------------
      SUBROUTINE close_group(gname,inid,errval)
      CHARACTER(*), INTENT(IN) :: gname
      INTEGER(HID_T), INTENT(IN) :: inid
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     Close group
!-----------------------------------------------------------------------
      CALL h5gclose_f(inid,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Error closing group: '//TRIM(gname)
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE close_group

!-----------------------------------------------------------------------
!     subprogram 5. test_group
!     See if group exists
!-----------------------------------------------------------------------
      SUBROUTINE test_group(inid,gname,group_exists,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER(*), INTENT(IN) :: gname
      LOGICAL, INTENT(OUT) :: group_exists
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
      INTEGER(HID_T) :: gid
!-----------------------------------------------------------------------
!     Determine whether group exists by trying to opening it and testing
!      error message
!-----------------------------------------------------------------------
      CALL h5gopen_f(inid,gname,gid,error)
      IF (error==FAIL) THEN
          group_exists=.FALSE.
      ELSE
          group_exists=.TRUE.
          CALL h5gclose_f(gid, error)
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE test_group

!-----------------------------------------------------------------------
!     subprogram 6. get_nmember
!     Get the number of members of a group
!-----------------------------------------------------------------------
      SUBROUTINE get_nmembers(inid,gname,nmembers,errval)
      CHARACTER(*), INTENT(IN) :: gname
      INTEGER(HID_T), INTENT(IN) :: inid
      INTEGER(HID_T), INTENT(OUT) :: nmembers
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     Determine whether group exists by trying to opening it and testing
!      error message
!-----------------------------------------------------------------------
      CALL h5gn_members_f(inid,gname,nmembers, error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Error in get_nmembers for'//gname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE get_nmembers

!-----------------------------------------------------------------------
!     subprogram 7. make reference
!     Simplify referencing of one object (source) to another (target)
!-----------------------------------------------------------------------
      SUBROUTINE make_reference(inid,in_gid,sname,tname,errval)
      INTEGER(HID_T), INTENT(IN) :: inid,in_gid
      CHARACTER*(*), INTENT(IN) :: sname,tname
      TYPE(hdf5ErrorType) :: errval
      INTEGER :: error

      INTEGER(HSIZE_T), DIMENSION(1) :: dimsr= (/4/)
      INTEGER :: trank = 1
      INTEGER(HID_T) :: type_id      ! Attribute Dataspace identifier
      INTEGER(HID_T) :: tgt_id,tgt_sid       ! Attribute identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
!-----------------------------------------------------------------------
!     See refobjexample.f90
!-----------------------------------------------------------------------
      !
      ! Create dataspace and dataset to store references to the objects
      !
      CALL h5screate_simple_f(trank, dimsr, tgt_sid, error)
      CALL h5dcreate_f(inid,tname,H5T_STD_REF_OBJ,tgt_sid,tgt_id,error)
      !
      ! Create a datatype and store in the file
      !
      CALL h5tcopy_f(H5T_NATIVE_REAL, type_id, error)
      CALL h5tcommit_f(inid, "MyType", type_id, error)

      errval%errBool = .false.
      RETURN
      END SUBROUTINE make_reference

!-----------------------------------------------------------------------
!     subprogram 10. write_attribute_ch_sc
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      SUBROUTINE write_attribute_ch_sc(inid,aname,attribute,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname,attribute
      TYPE(hdf5ErrorType) :: errval

      INTEGER :: arank=1
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
      IF (LEN_TRIM(attribute)==0) THEN
         errval%errBool = .false.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1
      attrlen=LEN_TRIM(attribute)

      ! Create the data space for the time attribute.
      CALL h5screate_f(H5S_SCALAR_F, aspace_id, error)
      !CALL h5screate_simple_f(arank, adims, aspace_id, error)

      ! Create datatype for the attribute.
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attrlen, error)

      ! Create dataset attribute for the group
      CALL h5acreate_f(inid, aname, atype_id, aspace_id,attr_id, error)

      IF (error==FAIL) THEN
         errval%errorMsg = 'Cannot create attribute '//aname//attribute
         errval%errBool = .true.
         RETURN
      ELSE
        ! Write the attribute data.
        CALL h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
        ! Close the attribute.
        CALL h5aclose_f(attr_id, error)
      ENDIF

      ! Close the dataspace.
      CALL h5sclose_f(aspace_id,error)

      errval%errBool = .false.
      RETURN
      END SUBROUTINE write_attribute_ch_sc

!-----------------------------------------------------------------------
!     subprogram 11. write_attribute_ch_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      SUBROUTINE write_attribute_ch_vec(inid,aname,attribute,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      CHARACTER*(*), DIMENSION(:), INTENT(IN) :: attribute
      TYPE(hdf5ErrorType) :: errval

      INTEGER :: arank=1
      INTEGER(i4) :: i
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims  ! Attribute dimension
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
      INTEGER(SIZE_T) :: ati
!-----------------------------------------------------------------------
!     If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
      IF (LEN_TRIM(attribute(1))==0) THEN
         errval%errBool = .false.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = SIZE(attribute)
      adims(:) = (/data_dims(1)/)
      attrlen=0
      DO i=1,data_dims(1)
         ! attrlen=MAX(attrlen,LEN(attribute(i)))
         ati = LEN(attribute(i))
         attrlen=MAX(attrlen,ati)
      ENDDO

      ! Create the data space for the time attribute.
      CALL h5screate_simple_f(arank, adims, aspace_id, error)

      ! Create datatype for the attribute.
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attrlen, error)

      ! Create dataset attribute for the group
      CALL h5acreate_f(inid, aname, atype_id, aspace_id,attr_id, error)

      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Can not create attribute '//aname
         errval%errBool = .true.
         RETURN
      ELSE
        ! Write the attribute data.
        CALL h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
        ! Close the attribute.
        CALL h5aclose_f(attr_id, error)
      ENDIF

      ! Close the dataspace.
      CALL h5sclose_f(aspace_id,error)

      errval%errBool = .false.
      RETURN
      END SUBROUTINE write_attribute_ch_vec

!-----------------------------------------------------------------------
!     subprogram 12. write_attribute_int_sc
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      SUBROUTINE write_attribute_int_sc(inid,aname,attribute,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      INTEGER(i4), INTENT(IN) :: attribute
      TYPE(hdf5ErrorType) :: errval

      INTEGER :: arank=1
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1

      ! Create the data space for the time attribute.
      CALL h5screate_simple_f(arank, adims, aspace_id, error)

      ! Create dataset attribute for the group
      CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
      CALL h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)

      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Can not create attribute '//aname
         errval%errBool = .true.
         RETURN
      ELSE
        ! Write the attribute data.
        CALL h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
        ! Close the attribute.
        CALL h5aclose_f(attr_id, error)
      ENDIF

      ! Close the dataspace.
      CALL h5sclose_f(aspace_id,error)

      errval%errBool = .false.
      RETURN
      END SUBROUTINE write_attribute_int_sc

!-----------------------------------------------------------------------
!     subprogram 13. write_attribute_int_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      SUBROUTINE write_attribute_int_vec(inid,aname,attribute,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      INTEGER(i4), DIMENSION(:), INTENT(IN) :: attribute
      TYPE(hdf5ErrorType) :: errval

      INTEGER :: arank=1
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims  ! Attribute dimension
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = SIZE(attribute)
      adims(:) = (/data_dims(1)/)

      ! Create the data space for the attribute.
      CALL h5screate_simple_f(arank, adims, aspace_id, error)

      ! Create dataset attribute for the group
      CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
      CALL h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)


      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Can not create attribute '//aname
         errval%errBool = .true.
         RETURN
      ELSE
        ! Write the attribute data.
        CALL h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
        ! Close the attribute.
        CALL h5aclose_f(attr_id, error)
      ENDIF

      ! Close the dataspace.
      CALL h5sclose_f(aspace_id,error)

      errval%errBool = .false.
      RETURN
      END SUBROUTINE write_attribute_int_vec

!-----------------------------------------------------------------------
!     subprogram 13. write_attribute_rl_sc
!-----------------------------------------------------------------------
      SUBROUTINE write_attribute_rl_sc(inid,aname,attribute,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      REAL(r8), INTENT(IN) :: attribute
      TYPE(hdf5ErrorType), INTENT(OUT) :: errval

      INTEGER :: arank=1
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1

      ! Create the data space for the time attribute.
      CALL h5screate_f(H5S_SCALAR_F, aspace_id, error)

      ! Create dataset attribute for the group

      CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
      CALL h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)

      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Can not create attribute '//aname
         errval%errBool = .true.
         RETURN
      ELSE
        ! Write the attribute data.
!-PRE        IF (h5in%typeConvert) THEN
!-PRE          CALL h5awrite_f(attr_id, atype_id, REAL(attribute,r4), data_dims, error)
!-PRE        ELSE
          CALL h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
!-PRE        ENDIF
        ! Close the attribute.
        CALL h5aclose_f(attr_id, error)
      ENDIF

      ! Close the dataspace.
      CALL h5sclose_f(aspace_id,error)

      errval%errBool = .false.
      RETURN
      END SUBROUTINE write_attribute_rl_sc

!-----------------------------------------------------------------------
!     subprogram 14. write_attribute_rl_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      SUBROUTINE write_attribute_rl_vec(inid,aname,attribute,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      REAL(r8), DIMENSION(:), INTENT(IN) :: attribute
      type(hdf5ErrorType) :: errval
      INTEGER :: arank=1
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: adims  ! Attribute dimension
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = SIZE(attribute)
      adims(:) = (/data_dims(1)/)

      ! Create the data space for the attribute.
      CALL h5screate_simple_f(arank, adims, aspace_id, error)

      ! Create dataset attribute for the group
      CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
      CALL h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)

      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Can not create attribute '//aname
         errval%errBool = .true.
         RETURN
      ELSE
        ! Write the attribute data.
!-PRE        IF(h5in%typeConvert) THEN
!-PRE          CALL h5awrite_f(attr_id, atype_id, REAL(attribute,r4), data_dims, error)
!-PRE        ELSE
          CALL h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
!-PRE        ENDIF
        ! Close the attribute.
        CALL h5aclose_f(attr_id, error)
      ENDIF

      ! Close the dataspace.
      CALL h5sclose_f(aspace_id,error)

      errval%errBool = .false.
      RETURN
      END SUBROUTINE write_attribute_rl_vec

!-----------------------------------------------------------------------
!     subprogram 20. make_mesh_group
!     Defines a mesh group that points to other variables that define
!      the actual mesh
!-----------------------------------------------------------------------
      SUBROUTINE make_mesh_group(gInId,gridId,h5in,meshName,&
                  meshKind,axis0,axis1,axis2,transform,trName,errval)
      INTEGER(HID_T), INTENT(IN) :: gInId
      INTEGER(HID_T), INTENT(INOUT) :: gridId
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      CHARACTER*(*), INTENT(IN) :: meshname,axis0,axis1,axis2
      CHARACTER*(*), INTENT(IN) :: meshKind,transform,trName
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     Open the group
!-----------------------------------------------------------------------
      CALL make_group(gInId, meshName, gridId,"",errval)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL write_attribute(gridId,'vsType',"mesh",errval)
      CALL write_attribute(gridId,'vsKind',meshKind,errval)
      CALL write_attribute(gridId,'vsAxis0',axis0,errval)
      CALL write_attribute(gridId,'vsAxis1',axis1,errval)
      CALL write_attribute(gridId,'vsAxis2',axis2,errval)
      IF(LEN_TRIM(transform)>0) THEN
         CALL write_attribute(gridId,'vsTransform',transform,errval)
         CALL write_attribute(gridId,'vsTransformedMesh',trName,errval)
      ENDIF
      IF(LEN_TRIM(h5in%vsCentering)>0) THEN
        CALL write_attribute(gridId,'vsCentering',h5in%vsCentering,&
                                  errval)
      ENDIF
!-----------------------------------------------------------------------
!     vsMD: Multidomain cabilities
!-----------------------------------------------------------------------
      IF(LEN_TRIM(h5in%vsMD)>0) THEN
         CALL write_attribute(gridId,"vsMD",h5in%vsMD,errval)
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE make_mesh_group
!-----------------------------------------------------------------------
!     subprogram 20. make_time_group
!     Make a group that contains the time data.  See:
!        https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
      SUBROUTINE make_time_group(gInId,grName,h5in,h5err)
      INTEGER(HID_T), INTENT(IN) :: gInId
      CHARACTER*(*), INTENT(IN) :: grName
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType), INTENT(INOUT) :: h5err
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER(HID_T) :: timeId
      INTEGER :: error
!-----------------------------------------------------------------------
!     Open the group
!-----------------------------------------------------------------------
      CALL make_group(gInId, grName, timeId,"",h5err)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL write_attribute(timeId,'vsType',"time",h5err)
      CALL write_attribute(timeId,'vsStep',h5in%vsStep,h5err)
      CALL write_attribute(timeId,'vsTime',h5in%vsTime,h5err)
      IF(h5in%unitConvert) THEN
         CALL write_attribute(timeId,'vsUnitConvert',h5in%vsUnitCnv,h5err)
      ENDIF
      CALL close_group(grName,timeId,h5err)
      RETURN
      END SUBROUTINE make_time_group
!-----------------------------------------------------------------------
!     subprogram 20. make_provenance_group
!     Make a group that contains the provenance data.  See:
!      https://ice.txcorp.com/trac/vizschema/wiki/ProvenanceMetaData
!-----------------------------------------------------------------------
      SUBROUTINE make_provenance_group(gInId,grName,prin,h5err)
      INTEGER(HID_T), INTENT(IN) :: gInId
      CHARACTER*(*), INTENT(IN) :: grName
      TYPE(provenanceData), INTENT(IN) :: prin
      TYPE(hdf5ErrorType), INTENT(INOUT) :: h5err
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER(HID_T) :: prId
      INTEGER :: error
!-----------------------------------------------------------------------
!     Open the group
!-----------------------------------------------------------------------
      CALL make_group(gInId, grName, prId,"",h5err)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes for provenance data
!  Att vsType = "runInfo"                             // Required string
!  Att vsSoftware = "VORPAL"                          // Required string
!  Att vsSwVersion = "3.1.0"                          // Required string
!  Att vsVsVersion = "2.1.0"                          // Required string
!  Att vsSwRevision = "r10168"                        // Optional string
!  Att vsCxxCompiler = "g++"                          // Optional string
!  Att vsCxxCompilerVersion = "4.3.0"                 // Optional string
!  Att vsCxxCompilerFlags = "-O2 -g"                  // Optional string
!  Att vsBuildHost = "boron.txcorp.com"               // Optional string
!  Att vsBuildHostType = "x86_64-unknown-linux-gnu"   // Optional string
!  Att vsRunHost = "boron.txcorp.com"                 // Optional string
!  Att vsUser = "cary"                                // Optional string
!  Att vsRunDate = "Mon Nov 10 09:06:14 MST 2008"     // Optional string
!  Att vsCommandLine = "vorpal -i run.in"             // Optional string
!-----------------------------------------------------------------------
      CALL write_attribute(prId,'vsType',"runInfo",h5err)
      IF(LEN_TRIM(prin%software)>0)  &
       CALL write_attribute(prId,'vsSoftware',prin%software,h5err)
      IF(LEN_TRIM(prin%swVersion)>0)  &
        CALL write_attribute(prId,'vsSwVersion',prin%swVersion,h5err)
      IF(LEN_TRIM(prin%swRevision)>0)  &
        CALL write_attribute(prId,'vsSwRevision',prin%swRevision,h5err)
      IF(LEN_TRIM(prin%vsVersion)>0)  &
        CALL write_attribute(prId,'vsVsVersion',prin%vsVersion,h5err)
      IF(LEN_TRIM(prin%fcCompiler)>0)  &
        CALL write_attribute(prId,'vsFcCompiler',prin%fcCompiler,h5err)
      IF(LEN_TRIM(prin%fcCompilerVersion)>0)  &
        CALL write_attribute(prId,'vsFcCompilerVersion',prin%fcCompilerVersion,h5err)
      IF(LEN_TRIM(prin%fcCompilerFlags)>0)  &
        CALL write_attribute(prId,'vsFcCompilerFlags',prin%fcCompilerFlags,h5err)
      IF(LEN_TRIM(prin%buildHost)>0)  &
        CALL write_attribute(prId,'vsBuildHost',prin%buildHost,h5err)
      IF(LEN_TRIM(prin%buildHostType)>0)  &
        CALL write_attribute(prId,'vsBuildHostType',prin%buildHostType,h5err)
      IF(LEN_TRIM(prin%runHost)>0)  &
        CALL write_attribute(prId,'vsRunHost',prin%runHost,h5err)
      IF(LEN_TRIM(prin%runHostType)>0)  &
        CALL write_attribute(prId,'vsRunHostType',prin%runHostType,h5err)
      IF(LEN_TRIM(prin%user)>0)  &
        CALL write_attribute(prId,'vsUser',prin%user,h5err)
      IF(LEN_TRIM(prin%runHost)>0)  &
        CALL write_attribute(prId,'vsRunHost',prin%runHost,h5err)
      IF(LEN_TRIM(prin%runDate)>0)  &
        CALL write_attribute(prId,'vsRunDate',prin%runDate,h5err)
      IF(LEN_TRIM(prin%commandLine)>0)  &
        CALL write_attribute(prId,'vsCommandLine',prin%commandLine,h5err)
      CALL close_group(grName,prId,h5err)
      RETURN
      END SUBROUTINE make_provenance_group
!-----------------------------------------------------------------------
!     subprogram 20. dump_h5in_attributes
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      SUBROUTINE dump_h5in_attributes(dset_id,h5in,h5err)
      INTEGER(HID_T), INTENT(IN) :: dset_id
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType), INTENT(INOUT) :: h5err
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      IF(LEN_TRIM(h5in%mesh)>0) THEN
       IF(h5in%mesh(1:5)=="mesh-") THEN
         IF(h5in%debug) WRITE(*,*) 'Writing vsType attributes'
         CALL write_attribute(dset_id,'vsType',"mesh",h5err)
         IF(h5in%debug) WRITE(*,*) 'Writing vsMesh attributes',h5in%mesh(6:)
         CALL write_attribute(dset_id,'vsKind',h5in%mesh(6:),h5err)
       ELSE
         IF(h5in%debug) WRITE(*,*) 'Writing vsType attributes'
         CALL write_attribute(dset_id,'vsType',"variable",h5err)
         IF(h5in%debug) WRITE(*,*) 'Writing vsMesh attributes',h5in%mesh
         CALL write_attribute(dset_id,'vsMesh',h5in%mesh,h5err)
         IF(LEN_TRIM(h5in%vsCentering)>0) THEN
           IF(h5in%debug) WRITE(*,*) 'Writing vsCentering attributes'
           CALL write_attribute(dset_id,'vsCentering',h5in%vsCentering, &
                                  h5err)
         ENDIF
       ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     vsMD: Multidomain cabilities
!-----------------------------------------------------------------------
      IF(LEN_TRIM(h5in%vsMD)>0) THEN
         IF(h5in%debug) WRITE(*,*) 'Writing vsMD attributes',h5in%vsMD
         CALL write_attribute(dset_id,"vsMD",h5in%vsMD,h5err)
      ENDIF
!-----------------------------------------------------------------------
!     Label the units
!-----------------------------------------------------------------------
      IF(LEN_TRIM(h5in%units)>0) THEN
         IF(h5in%debug) WRITE(*,*) 'Writing units attributes',h5in%units
         CALL write_attribute(dset_id,"units",h5in%units,h5err)
      ENDIF
!-----------------------------------------------------------------------
!     If we have the ability to write the conversion factors than do so
!-----------------------------------------------------------------------
      IF(h5in%unitConvert) THEN
         CALL write_attribute(dset_id,'vsUnitConvert',h5in%vsUnitCnv,h5err)
      ENDIF
!-----------------------------------------------------------------------
      IF(h5in%debug) WRITE(*,*) "Returning from h5in_attributes"
      RETURN
      END SUBROUTINE dump_h5in_attributes
!-----------------------------------------------------------------------
!     subprogram 20. dump_h5_int
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      SUBROUTINE dump_h5_int(inid,aname,value,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      INTEGER(i4), INTENT(IN) :: value
      TYPE(hdf5ErrorType) :: errval
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
      INTEGER(HID_T) :: dspace_id, dset_id
      INTEGER(HID_T) :: plist_id       ! Property list identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: dims=0
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname, &
                      H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,value,dims,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Data set write failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_h5_int
!-----------------------------------------------------------------------
!     subprogram 20. dump_h5
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      SUBROUTINE dump_h5_rldbl(inid,aname,value,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      REAL(r8), INTENT(IN) :: value
      TYPE(hdf5ErrorType) :: errval
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
      INTEGER(HID_T) :: dspace_id, dset_id
      INTEGER(HID_T) :: plist_id       ! Property list identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: dims=0
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      IF(h5in%typeConvert) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(value,r4),dims,error)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,value,dims,error)
      ENDIF
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Data set write failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_h5_rldbl
!-----------------------------------------------------------------------
!     subprogram 20. dump_h5_1d
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      SUBROUTINE dump_h5_1d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      REAL(r8), DIMENSION(:), INTENT(IN) :: array
      TYPE(hdf5ErrorType) :: errval
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
      INTEGER :: rank
      INTEGER(HID_T) :: dspace_id, dset_id
      INTEGER(HID_T) :: plist_id       ! Property list identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: dims
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 1;            dims(:) = (/SIZE(array,1)/)
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
!#ifdef __MPI
!       CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error, &
!                     xfer_prp = plist_id)
!#else
      IF(h5in%typeConvert) THEN
       CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(array,r4),dims,error)
      ELSE
       CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
      ENDIF
!#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Data set write failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_h5_1d
!-----------------------------------------------------------------------
!     subprogram 21. dump_h5_2d
!     Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
      SUBROUTINE dump_h5_2d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER(*), INTENT(IN) :: aname
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error

      INTEGER(HID_T) dspace_id, rank, dset_id, plist_id
      INTEGER(HSIZE_T) :: dims(2)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 2
      IF(h5in%doTranspose) THEN
         dims(2) = SIZE(array,1);  dims(1) = SIZE(array,2)
      ELSE
         dims(1) = SIZE(array,1);  dims(2) = SIZE(array,2)
      ENDIF
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
      IF(h5in%typeConvert) THEN
       IF(h5in%doTranspose) THEN
         CALL h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(REAL(array,r4)),dims, &
                       error,xfer_prp = plist_id)
       ELSE
         CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(array,r4),dims, &
                       error,xfer_prp = plist_id)
       ENDIF
      ELSE
       IF(h5in%doTranspose) THEN
         CALL h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                       error,xfer_prp = plist_id)
       ELSE
         CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                       error,xfer_prp = plist_id)
       ENDIF
      ENDIF
#else
      IF(h5in%typeConvert) THEN
       IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(REAL(array,r4)),dims, &
                       error)
       ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(array,r4),dims,error)
       ENDIF
      ELSE
       IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                       error)
       ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
       ENDIF
      ENDIF
#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Writing data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_h5_2d

!-----------------------------------------------------------------------
!     subprogram 22. dump_h5_3d
!-----------------------------------------------------------------------
      SUBROUTINE dump_h5_3d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER(*), INTENT(IN) :: aname
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error, i,j
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: tmparray

      INTEGER(HID_T) dspace_id, rank, dset_id, plist_id
      INTEGER(HSIZE_T) :: dims(3)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 3
      IF(h5in%doTranspose) THEN
       dims(3)=SIZE(array,1);dims(2)=SIZE(array,2);dims(1)=SIZE(array,3)
       ALLOCATE(tmparray(dims(1),dims(2),dims(3)))
       DO i=1,dims(1); DO j=1,dims(2)
          tmparray(i,j,:)=array(:,j,i)
       ENDDO; ENDDO
      ELSE
       dims(1)=SIZE(array,1);dims(2)=SIZE(array,2);dims(3)=SIZE(array,3)
      ENDIF
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      IF(h5in%debug) WRITE(*,*) 'Calling h5screate_simple_f', dims
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      IF(h5in%debug) WRITE(*,*) 'Calling h5dcreate_simple_f', h5in%wrd_type
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
      IF(h5in%typeConvert) THEN
       IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(tmparray,r4),dims, &
                       error,xfer_prp = plist_id)
        DEALLOCATE(tmparray)
       ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(array,r4),dims, &
                       error,xfer_prp = plist_id)
       ENDIF
      ELSE
       IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                       error,xfer_prp = plist_id)
        DEALLOCATE(tmparray)
       ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                       error,xfer_prp = plist_id)
       ENDIF
      ENDIF
#else
      IF(h5in%typeConvert) THEN
       IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(tmparray,r4),dims,error)
        DEALLOCATE(tmparray)
       ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(array,r4),dims,error)
       ENDIF
      ELSE
       IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
        DEALLOCATE(tmparray)
       ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
       ENDIF
      ENDIF
#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Writing data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_h5_3d
!-----------------------------------------------------------------------
!     subprogram 23. dump_h5_4d
!-----------------------------------------------------------------------
      SUBROUTINE dump_h5_4d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER(*), INTENT(IN) :: aname
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error, i,j,k
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: tmparray

      INTEGER(HID_T) dspace_id, rank, dset_id, plist_id
      INTEGER(HSIZE_T) :: dims(4)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 4
      IF(h5in%doTranspose) THEN
       dims(4)=SIZE(array,1);  dims(2)=SIZE(array,3)
       dims(3)=SIZE(array,2);  dims(1)=SIZE(array,4)
       write(*,*) dims
       ALLOCATE(tmparray(dims(1),dims(2),dims(3),dims(4)))
       DO i=1,dims(1); DO j=1,dims(2); DO k=1,dims(3)
          tmparray(i,j,k,:)=array(:,k,j,i)
      ENDDO; ENDDO; ENDDO
      ELSE
       dims(1)=SIZE(array,1);  dims(2)=SIZE(array,2)
       dims(3)=SIZE(array,3);  dims(4)=SIZE(array,4)
      ENDIF
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
      IF(h5in%typeConvert) THEN
        IF(h5in%doTranspose) THEN
          CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(tmparray,r4),dims, &
                         error,xfer_prp = plist_id)
          DEALLOCATE(tmparray)
        ELSE
          CALL h5dwrite_f(dset_id,h5in%wrd_type,REAL(array,r4),dims, &
                         error,xfer_prp = plist_id)
        ENDIF
      ELSE
        IF(h5in%doTranspose) THEN
          CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                         error,xfer_prp = plist_id)
          DEALLOCATE(tmparray)
        ELSE
          CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                         error,xfer_prp = plist_id)
        ENDIF
      ENDIF
#else
      IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
        DEALLOCATE(tmparray)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
      ENDIF
#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Writing data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_h5_4d
!-----------------------------------------------------------------------
!     subprogram 20. dump_rl_1d
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      SUBROUTINE dump_rl_1d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER*(*), INTENT(IN) :: aname
      REAL(r4), DIMENSION(:), INTENT(IN) :: array
      TYPE(hdf5ErrorType) :: errval
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error
      INTEGER :: rank
      INTEGER(HID_T) :: dspace_id, dset_id
      INTEGER(HID_T) :: plist_id       ! Property list identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: dims
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 1;            dims(:) = (/SIZE(array,1)/)
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
!#ifdef __MPI
!       CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error, &
!                     xfer_prp = plist_id)
!#else
       CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
!#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Data set write failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_rl_1d
!-----------------------------------------------------------------------
!     subprogram 21. dump_rl_2d
!     Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
      SUBROUTINE dump_rl_2d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER(*), INTENT(IN) :: aname
      REAL(r4), DIMENSION(:,:), INTENT(IN) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error

      INTEGER(HID_T) dspace_id, rank, dset_id, plist_id
      INTEGER(HSIZE_T) :: dims(2)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 2
      IF(h5in%doTranspose) THEN
         dims(2) = SIZE(array,1);  dims(1) = SIZE(array,2)
      ELSE
         dims(1) = SIZE(array,1);  dims(2) = SIZE(array,2)
      ENDIF
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
      IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                       error,xfer_prp = plist_id)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                       error,xfer_prp = plist_id)
      ENDIF
#else
      IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                       error)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
      ENDIF
#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Writing data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_rl_2d

!-----------------------------------------------------------------------
!     subprogram 22. dump_rl_3d
!-----------------------------------------------------------------------
      SUBROUTINE dump_rl_3d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER(*), INTENT(IN) :: aname
      REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error, i,j
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: tmparray

      INTEGER(HID_T) dspace_id, rank, dset_id, plist_id
      INTEGER(HSIZE_T) :: dims(3)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 3
      IF(h5in%doTranspose) THEN
       dims(3)=SIZE(array,1);dims(2)=SIZE(array,2);dims(1)=SIZE(array,3)
       ALLOCATE(tmparray(dims(1),dims(2),dims(3)))
       DO i=1,dims(1); DO j=1,dims(2)
          tmparray(i,j,:)=array(:,j,i)
       ENDDO; ENDDO
      ELSE
       dims(1)=SIZE(array,1);dims(2)=SIZE(array,2);dims(3)=SIZE(array,3)
      ENDIF
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
      IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                       error,xfer_prp = plist_id)
        DEALLOCATE(tmparray)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                       error,xfer_prp = plist_id)
      ENDIF
#else
      IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
        DEALLOCATE(tmparray)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
      ENDIF
#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Writing data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_rl_3d
!-----------------------------------------------------------------------
!     subprogram 23. dump_rl_4d
!-----------------------------------------------------------------------
      SUBROUTINE dump_rl_4d(inid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: inid
      CHARACTER(*), INTENT(IN) :: aname
      REAL(r4), DIMENSION(:,:,:,:), INTENT(IN) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error, i,j,k
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: tmparray

      INTEGER(HID_T) dspace_id, rank, dset_id, plist_id
      INTEGER(HSIZE_T) :: dims(4)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 4
      IF(h5in%doTranspose) THEN
       dims(4)=SIZE(array,1);  dims(3)=SIZE(array,2)
       dims(3)=SIZE(array,2);  dims(1)=SIZE(array,4)
       ALLOCATE(tmparray(dims(1),dims(2),dims(3),dims(4)))
       DO i=1,dims(1); DO j=1,dims(2); DO k=1,dims(3)
          tmparray(i,j,k,:)=array(:,k,j,i)
      ENDDO; ENDDO; ENDDO
      ELSE
       dims(1)=SIZE(array,1);  dims(2)=SIZE(array,2)
       dims(3)=SIZE(array,3);  dims(4)=SIZE(array,4)
      ENDIF
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      CALL h5screate_simple_f(rank,dims,dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Create the data set.
!     Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
      CALL h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Create data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       IF (error==FAIL) THEN
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          RETURN
       ENDIF
#endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
      IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                       error,xfer_prp = plist_id)
        DEALLOCATE(tmparray)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                       error,xfer_prp = plist_id)
      ENDIF
#else
      IF(h5in%doTranspose) THEN
        CALL h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
        DEALLOCATE(tmparray)
      ELSE
        CALL h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
      ENDIF
#endif
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Writing data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      CALL dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
      CALL h5sclose_f(dspace_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data space failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE dump_rl_4d




!-----------------------------------------------------------------------
!     subprogram 40. read_dims
!     Read the dimensions of dataset associated with aname: 1d array
!-----------------------------------------------------------------------
      SUBROUTINE read_dims(dset_id,aname,dims,errval)
      INTEGER(HID_T), INTENT(IN) :: dset_id
      CHARACTER*(*), INTENT(IN) :: aname
      INTEGER(HSIZE_T), DIMENSION(:), INTENT(INOUT) :: dims
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error

      INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: maxdims
      INTEGER(HID_T) dspace_id
!-----------------------------------------------------------------------
!     Get dataset's dataspace handle.
!-----------------------------------------------------------------------
      ALLOCATE(maxdims(SIZE(dims)))
      CALL h5dget_space_f(dset_id, dspace_id, error)
!-----------------------------------------------------------------------
!     Get dataspace's dimensinons.
!-----------------------------------------------------------------------
      CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE read_dims

!-----------------------------------------------------------------------
!     subprogram 41. read_h5_1d
!     Read simple data set: 1d array
!-----------------------------------------------------------------------
      SUBROUTINE read_h5_1d(fid,aname,array,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: fid
      CHARACTER*(*), INTENT(IN) :: aname
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: array
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER(HSIZE_T), DIMENSION(1) :: dims, fdims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error

      INTEGER(HID_T) dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) ' Reading 1d array: ', aname
      CALL h5dopen_f(fid, aname, dset_id, error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Find data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      dims(1)=SIZE(array)
      call read_dims(dset_id,aname,fdims,errval)
      call check_dims(dims,fdims, errval)
      if (errval%errBool) RETURN
      CALL h5dread_f(dset_id,h5in%wrd_type,array,dims,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Reading data set failed for '//aname
         errval%errBool = .true.
         CALL h5dclose_f(dset_id,error)
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE read_h5_1d

!-----------------------------------------------------------------------
!     subprogram 42. read_h5_2d
!     Read simple data set: 2d array
!-----------------------------------------------------------------------
      SUBROUTINE read_h5_2d(fid,aname,array,h5in,errval)
      INTEGER, INTENT(IN) :: fid
      CHARACTER*(*), INTENT(IN) :: aname
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER(HSIZE_T), DIMENSION(2) :: dims, fdims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error, i
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: tmparray

      INTEGER(HID_T) dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) ' Reading 2d array: ', aname
      CALL h5dopen_f(fid, aname, dset_id, error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Find data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Check dims
!-----------------------------------------------------------------------
      IF(.NOT. h5in%doTranspose) THEN
         dims(1)=SIZE(array,1); dims(2)=SIZE(array,2)
      ELSE
         dims(1)=SIZE(array,2); dims(2)=SIZE(array,1)
      ENDIF
      call read_dims(dset_id,aname,fdims,errval)
      call check_dims(dims,fdims, errval)
      IF (errval%errBool) THEN
        CALL h5dclose_f(dset_id,error)
        RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      IF(.NOT. h5in%doTranspose) THEN
         CALL h5dread_f(dset_id,h5in%wrd_type,array,dims,error)
      ELSE
         ALLOCATE(tmparray(dims(1),dims(2)))
         CALL h5dread_f(dset_id,h5in%wrd_type,tmparray,dims,error)
         DO i=1,dims(1);      array(:,i)=tmparray(i,:);     ENDDO
         DEALLOCATE(tmparray)
      ENDIF
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Read data set failed for '//aname
         errval%errBool = .true.
         CALL h5dclose_f(dset_id,error)
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE read_h5_2d

!-----------------------------------------------------------------------
!     subprogram 43. read_h5_3d
!     Read simple data set: 3d array
!-----------------------------------------------------------------------
      SUBROUTINE read_h5_3d(fid,aname,array,h5in,errval)
      INTEGER, INTENT(IN) :: fid
      CHARACTER*(*), INTENT(IN) :: aname
      REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: array
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER(HSIZE_T), DIMENSION(3) :: dims, fdims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error, i,j
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: tmparray

      INTEGER dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) ' Reading 3d array: ', aname
      CALL h5dopen_f(fid, TRIM(aname), dset_id, error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Find data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Check dims
!-----------------------------------------------------------------------
      IF(.NOT. h5in%doTranspose) THEN
         dims(1)=SIZE(array,1)
         dims(2)=SIZE(array,2)
         dims(3)=SIZE(array,3)
      ELSE
         dims(1)=SIZE(array,3)
         dims(2)=SIZE(array,2)
         dims(3)=SIZE(array,1)
      ENDIF
      call read_dims(dset_id,aname,fdims,errval)
      !call check_dims(dims,fdims, errval)
      IF (errval%errBool) THEN
        CALL h5dclose_f(dset_id,error)
        RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      IF(.NOT. h5in%doTranspose) THEN
         CALL h5dread_f(dset_id,h5in%wrd_type,array,dims,error)
      ELSE
         ALLOCATE(tmparray(dims(1),dims(2),dims(3)))
         CALL h5dread_f(dset_id,h5in%wrd_type,tmparray,dims,error)
         DO i=1,dims(1); DO j=1,dims(2)
          array(:,j,i)=tmparray(i,j,:)
         ENDDO; ENDDO
         DEALLOCATE(tmparray)
      ENDIF
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Read data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      CALL h5dclose_f(dset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE read_h5_3d
!-----------------------------------------------------------------------
!     subprogram 41. read_attribute_rl_sc
!     Read real scalar attribute
!-----------------------------------------------------------------------
      SUBROUTINE read_attribute_rl_sc(fid,aname,val,h5in,errval)
      INTEGER(HID_T), INTENT(IN) :: fid
      CHARACTER*(*), INTENT(IN) :: aname
      TYPE(hdf5InOpts), INTENT(IN) :: h5in
      double precision, INTENT(INOUT) :: val
      TYPE(hdf5ErrorType), INTENT(INOUT) :: errval
      INTEGER(HSIZE_T), DIMENSION(1) :: dims, fdims
      INTEGER,PARAMETER :: FAIL=-1
      INTEGER :: error

      INTEGER(HID_T) aset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      IF(h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      CALL h5aopen_name_f(fid, aname, aset_id, error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Find data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=1
      CALL h5aread_f(aset_id,h5in%wrd_type,val,dims,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Reading data set failed for '//aname
         errval%errBool = .true.
         CALL h5aclose_f(aset_id,error)
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      CALL h5aclose_f(aset_id,error)
      IF (error==FAIL) THEN
         errval%errorMsg = 'ERROR: Close data set failed for '//aname
         errval%errBool = .true.
         RETURN
      ENDIF
!-----------------------------------------------------------------------
      errval%errBool = .false.
      RETURN
      END SUBROUTINE read_attribute_rl_sc
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE hdf5_api
