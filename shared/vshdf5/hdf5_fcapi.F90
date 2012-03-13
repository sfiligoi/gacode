!-----------------------------------------------------------------------
! file hdf_api
! hdf_api module
!  Very generic module meant for writing HDF5 files with particular
!  attributes.  In case we want to convert over to C, it should make
!  this easier.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! code organization for hdf_api.
!-----------------------------------------------------------------------
! 1.  check_dims(dims, fdims, errval)
! 2.  vshdf5_fcinit()
! 3.  vshdf5_inith5vars(h5in, h5err)
! 4.  open_h5file(openmethod,fname,fileid,fdesc,rootgid,h5in,h5err)
! 5.  open_oldh5file(fname,fileid,fdesc,rootgid,h5in,h5err)
! 6.  open_newh5file(fname,fileId,fdesc,rootGid,h5in,h5err)
! 7.  close_h5file(fileId,root_id,h5err)
! 8.  open_group(inid,gname,gid,errval)
! 9.  make_group(inid,gname,gid,meshtitle,errval)
! 10. close_group(gname,inid,errval)
! 11. test_group(inid,gname,group_exists,errval)
! 12. get_nmembers(inid,gname,nmembers,errval)
! 13. make_reference(inid,in_gid,sname,tname,errval)
! 14. make_mesh_group(gInId,gridId,h5in,meshName,&
! 15. make_time_group(gInId,grName,h5in,h5err)
!
! 20. write_attribute_ch_sc(inid,aname,attribute,errval)
! 21. write_attribute_ch_vec(inid,aname,attribute,errval)
! 22. write_attribute_int_sc(inid,aname,attribute,errval)
! 23. write_attribute_int_vec(inid,aname,attribute,errval)
! 24. write_attribute_rl_sc(inid,aname,attribute,errval)
! 25. write_attribute_rl_vec(inid,aname,attribute,errval)
!
! 30. dump_h5in_attributes(dset_id,h5in,h5err)
! 31. dump_h5_int(inid,aname,value,h5in,errval)
! 32. dump_h5_rldbl(inid,aname,value,h5in,errval)
! 33. dump_h5_1d(inid,aname,array,h5in,errval)
! 34. dump_h5_2d(inid,aname,array,h5in,errval)
! 35. dump_h5_3d(inid,aname,array,h5in,errval)
! 36. dump_h5_4d(inid,aname,array,h5in,errval)
! 37. dump_rl_1d(inid,aname,array,h5in,errval)
! 38. dump_rl_2d(inid,aname,array,h5in,errval)
! 39. dump_rl_3d(inid,aname,array,h5in,errval)
! 39. dump_rl_4d(inid,aname,array,h5in,errval)
!
! 40. add_h5_int(inid,aname,value,h5in,errval)
! 42. add_h5_dbl(inid,aname,array,h5in,errval)
! 42. add_h5_int_1d(inid,aname,array,h5in,errval)
! 42. add_h5_1d(inid,aname,array,h5in,errval)
! 43. add_h5_2d(inid,aname,array,h5in,errval)
!
! 50. read_dims(dset_id,aname,dims,errval)
! 51. read_h5_1d(fid,aname,array,h5in,errval)
! 52. read_h5_2d(fid,aname,array,h5in,errval)
! 53. read_h5_3d(fid,aname,array,h5in,errval)
! 54. read_attribute_rl_sc(fid,aname,val,h5in,errval)
!-----------------------------------------------------------------------
! module hdf_api
!-----------------------------------------------------------------------
  module hdf5_api
  use hdf5
  implicit none
  character(5), parameter, private :: h5fortranapiversion="1.0"
  integer, parameter, private :: i4=selected_int_kind(9)
  integer, parameter, private :: i8=selected_int_kind(18)
  integer, parameter, private :: r4=selected_real_kind(6,37)
  integer, parameter, private :: r8=selected_real_kind(13,307)
!-----------------------------------------------------------------------
! Input parameters to control the attributes and how written out
! The method used to determine whether to write the character
! variables are written out is something like the following:
!    if(TRIM(h5in%vsMD)/="") then
!       call write_attribute(dset_id,"vsMD",h5in%vsMD,errval)
!    endif
! This needs to be paid attention to if writing out several
! variables in a row.
!
! MESH DISCUSSION:
! The key to making the vsSchema work is to associate fields with
! their meshes.  The way to point the mesh it use the mesh member
! the derived type below.  It is also used to define the mesh by
! prepending the mesh type from visSchemawith "mesh-"; e.g.,
!  h5in%mesh="mesh-structured" defines a mesh and
!  h5in%mesh="/coreMesh"       points to where the mesh is defined.
! For valid types of meshes, details on how the multi-domain
! specification works, and centering issues, see the visSchema wiki:
!  https://ice.txcorp.com/trac/vizschema/wiki/
!
! WRD_TYPE:  Valid wrd_types are:
!  H5T_NATIVE_REAL
!  H5T_NATIVE_DOUBLE
!  H5T_NATIVE_integer
!  H5T_NATIVE_INT16   !?
!-----------------------------------------------------------------------
! IMPORTANT:::
! It is very important to have codes use this API use the initvars
! When adding variables to this derived type, make sure they are 
!  initializing correctly.
!-----------------------------------------------------------------------
  type hdf5inopts
     integer(hid_t) :: wrd_type
     logical :: dotranspose        
     logical :: verbose               ! whether to write verbose output
     logical :: debug                 ! write even more verbose output for debugging
     logical :: pio                   ! whether parallel i/o is used
     logical :: wrvstime              ! whether to write time to attribute
     logical :: typeconvert           ! whether to demote the type
     logical :: unitconvert           ! whether to write vsunitcnv
     integer(selected_int_kind(9)) :: comm     ! communicator associated w/ open file
     integer(selected_int_kind(9)) :: info     ! mpi_info: set to mpi_info_null if unsure
     character(len=30) :: mesh             ! see above
     character(len=30) :: units            ! units
     character(len=30) :: vscentering      ! how to center variables on mesh
     character(len=30) :: vsmd             ! multidomain variable
     real(r8) :: vstime                    ! bassi requires r8
     integer(hid_t) :: vsstep=0            ! step # associated with time
     real(r8) :: vsunitcnv=1.              ! conversion factor
     ! double precision :: vstime          ! different b/c time dependent
  end type
!-----------------------------------------------------------------------
!     Example of how to use h5err type after an fcapi call.
!     if(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
!-----------------------------------------------------------------------
  type hdf5errortype
     logical :: errbool
     character(64) :: errormsg
  end type
!-----------------------------------------------------------------------
! subprogram name interfaces
!-----------------------------------------------------------------------
  interface write_attribute
    module procedure write_attribute_ch_sc,write_attribute_ch_vec &
                   ,write_attribute_rl_sc,write_attribute_rl_vec &
                   ,write_attribute_int_sc,write_attribute_int_vec
  end interface
  interface dump_h5
    module procedure  dump_h5_int,dump_h5_rldbl, &
                     dump_h5_1d,dump_h5_2d,dump_h5_3d,dump_h5_4d, &
                     dump_rl_1d,dump_rl_2d,dump_rl_3d,dump_rl_4d
  end interface
  ! This is like dump but does an append
  interface add_h5
    module procedure  add_h5_int, add_h5_dbl, &
       add_h5_int_1d, add_h5_1d, add_h5_2d, &
       add_h5_3d, add_h5_4d
  end interface
  ! dump components into separate groups
  interface read_h5
    module procedure read_h5_1d,read_h5_2d,read_h5_3d
  end interface
  interface read_attribute
    module procedure read_attribute_rl_sc
  end interface

  contains
!-----------------------------------------------------------------------
! subprogram -1. write mismatched dimension errors
!-----------------------------------------------------------------------
  subroutine check_dims(dims, fdims, errval)
  integer(hsize_t), dimension(:), intent(in) :: dims, fdims
  type(hdf5errortype), intent(inout) :: errval
     integer i
     do i = 1, size(dims)
       if (dims(i) /= fdims(i)) then
         write(*, *) "error: dims (", dims, ")"
         write(*, *) "  /=  fdims (", fdims, ")"
         write(*, *) "fdims = dims in the file (use h5ls)"
         write(*, *) "dims  = dims allocated for array to read"
         errval%errormsg = 'error: dims /= fdims'
         errval%errbool = .true.
         return
       endif
     enddo
     errval%errbool = .false.
  end subroutine check_dims

!-----------------------------------------------------------------------
! function 1. h5accessMethod
! Return something the correct hdf5 access method given a more
! memoral names
!-----------------------------------------------------------------------
  function h5accessmethod(access_method)
  integer(hid_t) :: h5accessmethod
  character(*), intent(in) :: access_method
  select case(access_method)
  case("overwr")
     h5accessmethod=h5f_acc_trunc_f                   ! overwrite file
  case("rdwr")
     h5accessmethod=h5f_acc_rdwr_f                    ! read-write
  case("rdonly")
     h5accessmethod=h5f_acc_rdonly_f                  ! read only
  end select
  return
  end function h5accessmethod

!-----------------------------------------------------------------------
! subprogram 0. vshdf5_fcinit
! Open fortran hdf5 and set open/close parameters.
!-----------------------------------------------------------------------
  subroutine vshdf5_fcinit()
  ! integer(hid_t) :: err
  integer :: err
  ! write(*, *) "vshdf5_fcinit: entered"
  call h5dont_atexit_f(err)
  ! write(*, *) "vshdf5_fcinit: h5dont_atexit_f returned."
  call h5open_f(err)
  ! write(*, *) "vshdf5_fcinit: h5open_f returned."
  ! write(*, *) "vshdf5_fcinit: leaving."
  return
  end subroutine vshdf5_fcinit

!-----------------------------------------------------------------------
! subprogram 0.1 vshdf5_inith5vars
! Initialize these variables to default values.
!-----------------------------------------------------------------------
  subroutine vshdf5_inith5vars(h5in, h5err)
  type(hdf5inopts), intent(inout) :: h5in
  type(hdf5errortype), intent(inout) :: h5err
!-----------------------------------------------------------------------
! Defaults
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
  return
  end subroutine vshdf5_inith5vars
!-----------------------------------------------------------------------
! subprogram 0. open_h5file
! Open file for writing and write file attributes
! This is just a nice wrapper for open_newh5file and open_oldh5file where
! we just specify the openmethod.  
! Separate subroutines kept because for simplicity in debugging
!-----------------------------------------------------------------------
  subroutine open_h5file(openmethod,fname,fileid,fdesc,rootgid,h5in,h5err)
  character(*), intent(in) :: fname,fdesc,openmethod
  integer(hid_t), intent(out) :: fileid,rootgid
  type(hdf5errortype), intent(inout) :: h5err
  type(hdf5inopts), intent(inout) :: h5in
  integer(hid_t) :: access_mode
  integer(hid_t) :: plist_id       ! property list identifier
  integer,parameter :: fail=-1
  integer :: info
  integer(hid_t) :: error

  logical :: file_exists
  if (h5in%verbose) then
    write(*, *) " open_h5file: entered."
  endif
!-----------------------------------------------------------------------
! Create and open the file
! The if statements seem to work better to give the expected
!  behavior for read-write for already open file.
!-----------------------------------------------------------------------
  select case(openmethod)
  case('overwr')
    call open_newh5file(fname,fileid,fdesc,rootgid,h5in,h5err)
  case('append')
    call open_oldh5file(fname,fileid,fdesc,rootgid,h5in,h5err)
  case default
    h5err%errorMsg = "open_h5file called with incorrect openmethod"
    h5err%errBool = .true.
    return
  end select
  return
  end subroutine open_h5file
!-----------------------------------------------------------------------
! subprogram 1. open_oldh5file
! Open file for writing and write file attributes
! Create the group for the independent variables at this stage
!-----------------------------------------------------------------------
  subroutine open_oldh5file(fname,fileid,fdesc,rootgid,h5in,h5err)
  character(*), intent(in) :: fname,fdesc
  integer(hid_t), intent(out) :: fileid,rootgid
  type(hdf5errortype), intent(inout) :: h5err
  type(hdf5inopts), intent(inout) :: h5in
  integer(hid_t) :: access_mode
  integer(hid_t) :: plist_id       ! property list identifier
  integer,parameter :: fail=-1
  integer :: info
  integer(hid_t) :: error

  logical :: file_exists
  if (h5in%verbose) then
    write(*, *) " open_oldh5file: entered."
  endif
!-----------------------------------------------------------------------
! Create and open the file
! The if statements seem to work better to give the expected
!  behavior for read-write for already open file.
!-----------------------------------------------------------------------
  inquire(file=trim(fname),exist=file_exists)
  if (.not. file_exists) then
     h5err%errormsg = 'error: file does not exist: '//fname
     if (h5in%verbose) then
       write(*, *) 'error: file does not exist: ', fname
     endif
     h5err%errbool = .true.
     return
  endif
!-------------------------------------------------------------------
! Setup file access property list with parallel I/O access.
! jrc 1jul09: Need not have parallel I/O access when only one proc is
! reading
!-----------------------------------------------------------------------
#ifdef __MPI
   call h5pcreate_f(h5p_file_access_f, plist_id, error)
   call h5pset_fapl_mpio_f(plist_id, h5in%comm, h5in%info, error)
! jrc 12jul09: i think these should be set outside.  a parallel
! program may not want to use parallel i/o, e.g., if it is doing i/o
! on only one rank.  scott k, what do you think?
   h5in%pio=.true.
#endif

!-----------------------------------------------------------------------
! Open file -- collectively if needed
!-----------------------------------------------------------------------
  if (h5in%verbose) then
    write(*, *) "open_oldh5file: calling h5fopen_f."
  endif
#ifdef __MPI
  if (h5in%pIO) then
    call h5fopen_f(fname,H5F_ACC_RDWR_F,fileId,error, &
                  access_prp=plist_id)
  else
#else
  call h5fopen_f(fname, H5F_ACC_RDWR_F, fileId, error)
#endif
#ifdef __MPI
  endif
#endif
  if (h5in%verbose) then
    write(*, *) "open_oldh5file: h5fopen_f returned."
  endif

!-----------------------------------------------------------------------
! Grab the root group id which is created by default
!-----------------------------------------------------------------------
  call h5gopen_f(fileId,"/",rootGid,error)
  if (error==FAIL) then
     h5err%errorMsg = 'ERROR: Error grabbing root ID: '//fname
     h5err%errBool = .true.
     return
  endif
  h5err%errBool = .false.
  return
  end subroutine open_oldh5file

!-----------------------------------------------------------------------
! subprogram 1. open_newh5file
! Open file for writing and write file attributes
! Create the group for the independent variables at this stage
!-----------------------------------------------------------------------
  subroutine open_newh5file(fname,fileId,fdesc,rootGid,h5in,h5err)
  character(*), intent(in) :: fname,fdesc
  integer(HID_T), intent(out) :: fileId,rootGid
  TYPE(hdf5InOpts), intent(inoUT) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: h5err
  integer(HID_T) :: access_mode
  integer,parameter :: FAIL=-1
  integer(HID_T) :: plist_id       ! Property list identifier
  integer :: error, info

  LOGICAL :: file_exists
!-----------------------------------------------------------------------
!    Setup file access property list with parallel I/O access.
!-----------------------------------------------------------------------
#ifdef __MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, h5in%comm, h5in%info, error)
    h5in%pIO=.true.
#endif
!-----------------------------------------------------------------------
! Open file -- collectively if needed
! The if statements seem to work better to give the expected
!  behavior for read-write for already open file.
!-----------------------------------------------------------------------
  INQUIRE(FILE=TRIM(fname),EXIST=file_exists)
  if (.NOT. file_exists) then
     ! Always create with over-write to avoid errors
     access_mode=h5accessMethod("overwr")
#ifdef __MPI
     call h5fcreate_f(TRIM(fname),access_mode,fileId,error, &
                  access_prp=plist_id)
#else
     call h5fcreate_f(TRIM(fname),access_mode,fileId,error)
#endif
  else
     OPEN(UNIT=999,FILE=TRIM(fname),FORM='UNFORMATTED', &
       POSITION='REWIND',STATUS='REPLACE')
     CLOSE(UNIT=999)
     access_mode=h5accessMethod("overwr")
#ifdef __MPI
     call h5fcreate_f(TRIM(fname),access_mode,fileId,error, &
                  access_prp=plist_id)
#else
     call h5fcreate_f(TRIM(fname),access_mode,fileId,error)
#endif
  endif
  if (error==FAIL) then
     h5err%errorMsg = 'ERROR: Error opening file: '//fname
     h5err%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Grab the root group id which is created by default
!-----------------------------------------------------------------------
  call h5gopen_f(fileId,"/",rootGid,error)
  if (error==FAIL) then
     h5err%errorMsg = 'ERROR: Error grabbing root ID: '//fname
     h5err%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Put in things which describe the api
!-----------------------------------------------------------------------
!  call write_attribute(rootGid,"vsFcVERSION","1.0",h5err)
!-----------------------------------------------------------------------
! Put in file description when creating file
!-----------------------------------------------------------------------
  call write_attribute(rootGid,"Description",fdesc,h5err)
!-----------------------------------------------------------------------
  h5err%errBool = .false.
  return
  end subroutine open_newh5file

!-----------------------------------------------------------------------
! subprogram 2. close_h5file
! Close the file associated with fileId.
!-----------------------------------------------------------------------
  subroutine close_h5file(fileId,root_id,h5err)
  integer, intent(in) :: fileId,root_id
  TYPE(hdf5ErrorType), intent(inoUT) :: h5err
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! Close the root group and file
!-----------------------------------------------------------------------
  call h5gclose_f(root_id, error)
  call h5fclose_f(fileId, error)
  if (error==FAIL) then
     h5err%errorMsg = 'ERROR: Error in close_h5file'
     h5err%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  h5err%errBool = .false.
  return
  end subroutine close_h5file

!-----------------------------------------------------------------------
! subprogram 3. open_group
! Open a group in a safe way
!-----------------------------------------------------------------------
  subroutine open_group(inid,gname,gid,errval)
  character(*), intent(in) :: gname
  integer(HID_T), intent(in) :: inid
  integer(HID_T), intent(out) :: gid
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer,parameter :: FAIL=-1
  integer :: error
!-------------------------------------------------------------------
! Open group
!-----------------------------------------------------------------------
  call h5gopen_f(inid,gname,gid,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Error opening group: '//gname
     errval%errBool = .true.
     return
  endif
  errval%errBool = .false.
  return
  end subroutine open_group

!-----------------------------------------------------------------------
! subprogram 3. make_group
! Create a group in a safe way
!-----------------------------------------------------------------------
  subroutine make_group(inid,gname,gid,meshtitle,errval)
  character(*), intent(in) :: gname
  character*(*), intent(in) :: meshtitle
  integer(HID_T), intent(in) :: inid
  integer(HID_T), intent(out) :: gid
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error

  character(64) :: msg
!-----------------------------------------------------------------------
! Create group
!-----------------------------------------------------------------------
  !call h5gopen_f(inid,gname,gid,error)
  ! If it failed, most likely it doesn't exist
  call h5gcreate_f(inid,gname,gid,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Error opening group: '//gname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
!  call write_attribute(gid,'CLASS','GROUP',errval)
!  call write_attribute(gid,'VERSION','1.0',errval)
!   call write_attribute(gid,'TITLE',TRIM(meshtitle),errval)
  if(len_trim(meshtitle)>0) then
   call write_attribute(gid,'vsType','mesh',errval)
   call write_attribute(gid,'vsKind','structured',errval)
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine make_group

!-----------------------------------------------------------------------
! subprogram 4. close_group
! Close a group in a safe way
!-----------------------------------------------------------------------
  subroutine close_group(gname,inid,errval)
  character(*), intent(in) :: gname
  integer(HID_T), intent(in) :: inid
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! Close group
!-----------------------------------------------------------------------
  call h5gclose_f(inid,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Error closing group: '//TRIM(gname)
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine close_group

!-----------------------------------------------------------------------
! subprogram 5. test_group
! See if group exists
!-----------------------------------------------------------------------
  subroutine test_group(inid,gname,group_exists,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: gname
  LOGICAL, intent(out) :: group_exists
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer,parameter :: FAIL=-1
  integer :: error
  integer(HID_T) :: gid
!-----------------------------------------------------------------------
! Determine whether group exists by trying to opening it and testing
!  error message
!-----------------------------------------------------------------------
  call h5gopen_f(inid,gname,gid,error)
  if (error==FAIL) then
      group_exists=.FALSE.
  else
      group_exists=.TRUE.
      call h5gclose_f(gid, error)
  endif
!-------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine test_group

!-----------------------------------------------------------------------
! subprogram 6. get_nmember
! Get the number of members of a group
!-----------------------------------------------------------------------
  subroutine get_nmembers(inid,gname,nmembers,errval)
  character(*), intent(in) :: gname
  integer(HID_T), intent(in) :: inid
  integer(HID_T), intent(out) :: nmembers
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! Determine whether group exists by trying to opening it and testing
!  error message
!-----------------------------------------------------------------------
  call h5gn_members_f(inid,gname,nmembers, error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Error in get_nmembers for'//gname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine get_nmembers

!-----------------------------------------------------------------------
! subprogram 7. make reference
! Simplify referencing of one object (source) to another (target)
!-----------------------------------------------------------------------
  subroutine make_reference(inid,in_gid,sname,tname,errval)
  integer(HID_T), intent(in) :: inid,in_gid
  character*(*), intent(in) :: sname,tname
  TYPE(hdf5ErrorType) :: errval
  integer :: error

  integer(HSIZE_T), dimension(1) :: dimsr= (/4/)
  integer :: trank = 1
  integer(HID_T) :: type_id      ! Attribute Dataspace identifier
  integer(HID_T) :: tgt_id,tgt_sid       ! Attribute identifier
  integer(HSIZE_T), dimension(1) :: adims = (/1/) ! Attribute dimension
!-----------------------------------------------------------------------
! See refobjexample.f90
!-----------------------------------------------------------------------
  !
  ! Create dataspace and dataset to store references to the objects
  !
  call h5screate_simple_f(trank, dimsr, tgt_sid, error)
  call h5dcreate_f(inid,tname,H5T_STD_REF_OBJ,tgt_sid,tgt_id,error)
  !
  ! Create a datatype and store in the file
  !
  call h5tcopy_f(H5T_NATIVE_REAL, type_id, error)
  call h5tcommit_f(inid, "MyType", type_id, error)

  errval%errBool = .false.
  return
  end subroutine make_reference
!-----------------------------------------------------------------------
! subprogram 20. make_mesh_group
! Defines a mesh group that points to other variables that define
!  the actual mesh
!-----------------------------------------------------------------------
  subroutine make_mesh_group(gInId,gridId,h5in,meshName,&
              meshKind,axis0,axis1,axis2,transform,trName,errval)
  integer(HID_T), intent(in) :: gInId
  integer(HID_T), intent(inoUT) :: gridId
  TYPE(hdf5InOpts), intent(in) :: h5in
  character*(*), intent(in) :: meshname,axis0,axis1,axis2
  character*(*), intent(in) :: meshKind,transform,trName
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! Open the group
!-----------------------------------------------------------------------
  call make_group(gInId, meshName, gridId,"",errval)
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call write_attribute(gridId,'vsType',"mesh",errval)
  call write_attribute(gridId,'vsKind',meshKind,errval)
  call write_attribute(gridId,'vsAxis0',axis0,errval)
  call write_attribute(gridId,'vsAxis1',axis1,errval)
  call write_attribute(gridId,'vsAxis2',axis2,errval)
  if(len_trim(transform)>0) then
     call write_attribute(gridId,'vsTransform',transform,errval)
     call write_attribute(gridId,'vsTransformedMesh',trName,errval)
  endif
  if(len_trim(h5in%vsCentering)>0) then
    call write_attribute(gridId,'vsCentering',h5in%vsCentering,&
                              errval)
  endif
!-----------------------------------------------------------------------
! vsMD: Multidomain cabilities
!-----------------------------------------------------------------------
  if(len_trim(h5in%vsMD)>0) then
     call write_attribute(gridId,"vsMD",h5in%vsMD,errval)
  endif
!-----------------------------------------------------------------------
  return
  end subroutine make_mesh_group
!-----------------------------------------------------------------------
! subprogram 20. make_time_group
! Make a group that contains the time data.  See:
!    https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
  subroutine make_time_group(gInId,grName,h5in,h5err)
  integer(HID_T), intent(in) :: gInId
  character*(*), intent(in) :: grName
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: h5err
  integer,parameter :: FAIL=-1
  integer(HID_T) :: timeId
  integer :: error
!-----------------------------------------------------------------------
! Open the group
!-----------------------------------------------------------------------
  call make_group(gInId, grName, timeId,"",h5err)
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call write_attribute(timeId,'vsType',"time",h5err)
  call write_attribute(timeId,'vsStep',h5in%vsStep,h5err)
  call write_attribute(timeId,'vsTime',h5in%vsTime,h5err)
  if(h5in%unitConvert) then
     call write_attribute(timeId,'vsUnitConvert',h5in%vsUnitCnv,h5err)
  endif
  call close_group(grName,timeId,h5err)
  return
  end subroutine make_time_group

!-----------------------------------------------------------------------
! subprogram 10. write_attribute_ch_sc
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_ch_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname,attribute
  TYPE(hdf5ErrorType) :: errval

  integer :: arank=1
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
  integer(HID_T) :: attr_id       ! Attribute identifier
  integer(HSIZE_T), dimension(1) :: adims = (/1/) ! Attribute dimension
  integer(SIZE_T) :: attrlen    ! Length of the attribute string
  integer(HSIZE_T), dimension(1) :: data_dims
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
  if (len_trim(attribute)==0) then
     errval%errBool = .false.
     return
  endif
!-----------------------------------------------------------------------
! See attrexample.f90
!-----------------------------------------------------------------------
  data_dims(1) = 1
  attrlen=len_trim(attribute)

  ! Create the data space for the time attribute.
  call h5screate_f(H5S_SCALAR_F, aspace_id, error)
  !call h5screate_simple_f(arank, adims, aspace_id, error)

  ! Create datatype for the attribute.
  call h5tcopy_f(H5T_NATIVE_character, atype_id, error)
  call h5tset_size_f(atype_id, attrlen, error)

  ! Create dataset attribute for the group
  call h5acreate_f(inid, aname, atype_id, aspace_id,attr_id, error)

  if (error==FAIL) then
     errval%errorMsg = 'Cannot create attribute '//aname//attribute
     errval%errBool = .true.
     return
  else
    ! Write the attribute data.
    call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
    ! Close the attribute.
    call h5aclose_f(attr_id, error)
  endif

  ! Close the dataspace.
  call h5sclose_f(aspace_id,error)

  errval%errBool = .false.
  return
  end subroutine write_attribute_ch_sc

!-----------------------------------------------------------------------
! subprogram 11. write_attribute_ch_vec
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_ch_vec(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  character*(*), dimension(:), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  integer :: arank=1
  integer(i4) :: i
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
  integer(HID_T) :: attr_id       ! Attribute identifier
  integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
  integer(SIZE_T) :: attrlen    ! Length of the attribute string
  integer(HSIZE_T), dimension(1) :: data_dims
  integer,parameter :: FAIL=-1
  integer :: error
  integer(SIZE_T) :: ati
!-----------------------------------------------------------------------
! If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
  if (len_trim(attribute(1))==0) then
     errval%errBool = .false.
     return
  endif
!-----------------------------------------------------------------------
! See attrexample.f90
!-----------------------------------------------------------------------
  data_dims(1) = SIZE(attribute)
  adims(:) = (/data_dims(1)/)
  attrlen=0
  do i=1,data_dims(1)
     ! attrlen=MAX(attrlen,LEN(attribute(i)))
     ati = LEN(attribute(i))
     attrlen=MAX(attrlen,ati)
  enddo

  ! Create the data space for the time attribute.
  call h5screate_simple_f(arank, adims, aspace_id, error)

  ! Create datatype for the attribute.
  call h5tcopy_f(H5T_NATIVE_character, atype_id, error)
  call h5tset_size_f(atype_id, attrlen, error)

  ! Create dataset attribute for the group
  call h5acreate_f(inid, aname, atype_id, aspace_id,attr_id, error)

  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Can not create attribute '//aname
     errval%errBool = .true.
     return
  else
    ! Write the attribute data.
    call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
    ! Close the attribute.
    call h5aclose_f(attr_id, error)
  endif

  ! Close the dataspace.
  call h5sclose_f(aspace_id,error)

  errval%errBool = .false.
  return
  end subroutine write_attribute_ch_vec

!-----------------------------------------------------------------------
! subprogram 12. write_attribute_int_sc
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_int_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i4), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  integer :: arank=1
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
  integer(HID_T) :: attr_id       ! Attribute identifier
  integer(HSIZE_T), dimension(1) :: adims = (/1/) ! Attribute dimension
  integer(HSIZE_T), dimension(1) :: data_dims
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! See attrexample.f90
!-----------------------------------------------------------------------
  data_dims(1) = 1

  ! Create the data space for the time attribute.
  call h5screate_simple_f(arank, adims, aspace_id, error)

  ! Create dataset attribute for the group
  call h5tcopy_f(H5T_NATIVE_integer, atype_id, error)
  call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)

  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Can not create attribute '//aname
     errval%errBool = .true.
     return
  else
    ! Write the attribute data.
    call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
    ! Close the attribute.
    call h5aclose_f(attr_id, error)
  endif

  ! Close the dataspace.
  call h5sclose_f(aspace_id,error)

  errval%errBool = .false.
  return
  end subroutine write_attribute_int_sc

!-----------------------------------------------------------------------
! subprogram 13. write_attribute_int_vec
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_int_vec(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i4), dimension(:), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  integer :: arank=1
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
  integer(HID_T) :: attr_id       ! Attribute identifier
  integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
  integer(HSIZE_T), dimension(1) :: data_dims
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! See attrexample.f90
!-----------------------------------------------------------------------
  data_dims(1) = SIZE(attribute)
  adims(:) = (/data_dims(1)/)

  ! Create the data space for the attribute.
  call h5screate_simple_f(arank, adims, aspace_id, error)

  ! Create dataset attribute for the group
  call h5tcopy_f(H5T_NATIVE_integer, atype_id, error)
  call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)


  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Can not create attribute '//aname
     errval%errBool = .true.
     return
  else
    ! Write the attribute data.
    call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
    ! Close the attribute.
    call h5aclose_f(attr_id, error)
  endif

  ! Close the dataspace.
  call h5sclose_f(aspace_id,error)

  errval%errBool = .false.
  return
  end subroutine write_attribute_int_vec

!-----------------------------------------------------------------------
! subprogram 13. write_attribute_rl_sc
!-----------------------------------------------------------------------
  subroutine write_attribute_rl_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), intent(in) :: attribute
  TYPE(hdf5ErrorType), intent(out) :: errval

  integer :: arank=1
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
  integer(HID_T) :: attr_id       ! Attribute identifier
  integer(HSIZE_T), dimension(1) :: adims = (/1/) ! Attribute dimension
  integer(HSIZE_T), dimension(1) :: data_dims
  integer,parameter :: FAIL=-1
  integer :: error
!-----------------------------------------------------------------------
! See attrexample.f90
!-----------------------------------------------------------------------
  data_dims(1) = 1

  ! Create the data space for the time attribute.
  call h5screate_f(H5S_SCALAR_F, aspace_id, error)

  ! Create dataset attribute for the group

  call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
  call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)

  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Can not create attribute '//aname
     errval%errBool = .true.
     return
  else
    ! Write the attribute data.
!-PRE        if (h5in%typeConvert) then
!-PRE          call h5awrite_f(attr_id, atype_id, real(attribute,r4), data_dims, error)
!-PRE        else
      call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
!-PRE        endif
    ! Close the attribute.
    call h5aclose_f(attr_id, error)
  endif

  ! Close the dataspace.
  call h5sclose_f(aspace_id,error)

  errval%errBool = .false.
  return
  end subroutine write_attribute_rl_sc

!-----------------------------------------------------------------------
! subprogram 14. write_attribute_rl_vec
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_rl_vec(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), dimension(:), intent(in) :: attribute
  type(hdf5ErrorType) :: errval
  integer :: arank=1
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
  integer(HID_T) :: attr_id       ! Attribute identifier
  integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
  integer(HSIZE_T), dimension(1) :: data_dims
  integer,parameter :: FAIL=-1
  integer :: error
!-------------------------------------------------------------------
! See attrexample.f90
!-------------------------------------------------------------------
  data_dims(1) = SIZE(attribute)
  adims(:) = (/data_dims(1)/)

  ! Create the data space for the attribute.
  call h5screate_simple_f(arank, adims, aspace_id, error)

  ! Create dataset attribute for the group
  call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
  call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)

  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Can not create attribute '//aname
     errval%errBool = .true.
     return
  else
    ! Write the attribute data.
!-PRE        if(h5in%typeConvert) then
!-PRE          call h5awrite_f(attr_id, atype_id, real(attribute,r4), data_dims, error)
!-PRE        else
      call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
!-PRE        endif
    ! Close the attribute.
    call h5aclose_f(attr_id, error)
  endif

  ! Close the dataspace.
  call h5sclose_f(aspace_id,error)

  errval%errBool = .false.
  return
  end subroutine write_attribute_rl_vec

!-----------------------------------------------------------------------
! subprogram 20. dump_h5in_attributes
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_h5in_attributes(dset_id,h5in,h5err)
  integer(HID_T), intent(in) :: dset_id
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: h5err
  integer,parameter :: FAIL=-1
  integer :: error
!-------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  if(len_trim(h5in%mesh)>0) then
   if(h5in%mesh(1:5)=="mesh-") then
     if(h5in%debug) WRITE(*,*) 'Writing vsType attributes'
     call write_attribute(dset_id,'vsType',"mesh",h5err)
     if(h5in%debug) WRITE(*,*) 'Writing vsMesh attributes',h5in%mesh(6:)
     call write_attribute(dset_id,'vsKind',h5in%mesh(6:),h5err)
   else
     if(h5in%debug) WRITE(*,*) 'Writing vsType attributes'
     call write_attribute(dset_id,'vsType',"variable",h5err)
     if(h5in%debug) WRITE(*,*) 'Writing vsMesh attributes',h5in%mesh
     call write_attribute(dset_id,'vsMesh',h5in%mesh,h5err)
     if(len_trim(h5in%vsCentering)>0) then
       if(h5in%debug) WRITE(*,*) 'Writing vsCentering attributes'
       call write_attribute(dset_id,'vsCentering',h5in%vsCentering, &
                              h5err)
     endif
   endif
  endif
!-----------------------------------------------------------------------
! vsMD: Multidomain cabilities
!-----------------------------------------------------------------------
  if(len_trim(h5in%vsMD)>0) then
     if(h5in%debug) WRITE(*,*) 'Writing vsMD attributes',h5in%vsMD
     call write_attribute(dset_id,"vsMD",h5in%vsMD,h5err)
  endif
!-----------------------------------------------------------------------
! Label the units
!-----------------------------------------------------------------------
  if(len_trim(h5in%units)>0) then
     if(h5in%debug) WRITE(*,*) 'Writing units attributes',h5in%units
     call write_attribute(dset_id,"units",h5in%units,h5err)
  endif
!-----------------------------------------------------------------------
! If we have the ability to write the conversion factors than do so
!-----------------------------------------------------------------------
  if(h5in%unitConvert) then
     call write_attribute(dset_id,'vsUnitConvert',h5in%vsUnitCnv,h5err)
  endif
!-----------------------------------------------------------------------
  if(h5in%debug) WRITE(*,*) "Returning from h5in_attributes"
  return
  end subroutine dump_h5in_attributes
!-----------------------------------------------------------------------
! subprogram 20. dump_h5_int
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_h5_int(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i4), intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: error
  integer(HID_T) :: dspace_id, dset_id
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(1) :: dims=0
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_f(H5S_SCALAR_F, dspace_id, error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname, &
                  H5T_NATIVE_integer,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
  call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,value,dims,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Data set write failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_h5_int
!-------------------------------------------------------------------
! subprogram 20. dump_h5
! Write an hdf5 array + references to independent vars
!-------------------------------------------------------------------
  subroutine dump_h5_rldbl(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: error
  integer(HID_T) :: dspace_id, dset_id
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(1) :: dims=0
!-------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-------------------------------------------------------------------
! Create the data space.
!-------------------------------------------------------------------
  call h5screate_f(H5S_SCALAR_F, dspace_id, error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
       if (error==FAIL) then
          errval%errorMsg = 'ERROR: Creating plist failed for '//aname
          errval%errBool = .true.
          return
       endif
#endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
  if(h5in%typeConvert) then
    call h5dwrite_f(dset_id,h5in%wrd_type,real(value,r4),dims,error)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,value,dims,error)
  endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Data set write failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_h5_rldbl
!-----------------------------------------------------------------------
! subprogram 20. dump_h5_1d
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_h5_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: error
  integer :: rank
  integer(HID_T) :: dspace_id, dset_id
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(1) :: dims
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 1;            dims(:) = (/SIZE(array,1)/)
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
!#ifdef __MPI
!   call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error, &
!                 xfer_prp = plist_id)
!#else
  if(h5in%typeConvert) then
   call h5dwrite_f(dset_id,h5in%wrd_type,real(array,r4),dims,error)
  else
   call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
  endif
!#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Data set write failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_h5_1d
!-----------------------------------------------------------------------
! subprogram 21. dump_h5_2d
! Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
  subroutine dump_h5_2d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r8), dimension(:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error

  integer(HID_T) dspace_id, rank, dset_id, plist_id
  integer(HSIZE_T) :: dims(2)
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 2
  if(h5in%doTranspose) then
     dims(2) = SIZE(array,1);  dims(1) = SIZE(array,2)
  else
     dims(1) = SIZE(array,1);  dims(2) = SIZE(array,2)
  endif
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
  if(h5in%typeConvert) then
   if(h5in%doTranspose) then
     call h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(real(array,r4)),dims, &
                   error,xfer_prp = plist_id)
   else
     call h5dwrite_f(dset_id,h5in%wrd_type,real(array,r4),dims, &
                   error,xfer_prp = plist_id)
   endif
  else
   if(h5in%doTranspose) then
     call h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                   error,xfer_prp = plist_id)
   else
     call h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                   error,xfer_prp = plist_id)
   endif
  endif
#else
  if(h5in%typeConvert) then
   if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(real(array,r4)),dims, &
                   error)
   else
    call h5dwrite_f(dset_id,h5in%wrd_type,real(array,r4),dims,error)
   endif
  else
   if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                   error)
   else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
   endif
  endif
#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Writing data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_h5_2d

!-----------------------------------------------------------------------
! subprogram 22. dump_h5_3d
!-----------------------------------------------------------------------
  subroutine dump_h5_3d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r8), dimension(:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error, i,j
  real(r8), dimension(:,:,:), allocatable :: tmparray

  integer(HID_T) dspace_id, rank, dset_id, plist_id
  integer(HSIZE_T) :: dims(3)
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 3
  if(h5in%doTranspose) then
   dims(3)=SIZE(array,1);dims(2)=SIZE(array,2);dims(1)=SIZE(array,3)
   allocate(tmparray(dims(1),dims(2),dims(3)))
   do i=1,dims(1); do j=1,dims(2)
      tmparray(i,j,:)=array(:,j,i)
   enddo; enddo
  else
   dims(1)=SIZE(array,1);dims(2)=SIZE(array,2);dims(3)=SIZE(array,3)
  endif
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  if(h5in%debug) WRITE(*,*) 'Calling h5screate_simple_f', dims
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  if(h5in%debug) WRITE(*,*) 'Calling h5dcreate_simple_f', h5in%wrd_type
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
  if(h5in%typeConvert) then
   if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,real(tmparray,r4),dims, &
                   error,xfer_prp = plist_id)
    deallocate(tmparray)
   else
    call h5dwrite_f(dset_id,h5in%wrd_type,real(array,r4),dims, &
                   error,xfer_prp = plist_id)
   endif
  else
   if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                   error,xfer_prp = plist_id)
    deallocate(tmparray)
   else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                   error,xfer_prp = plist_id)
   endif
  endif
#else
  if(h5in%typeConvert) then
   if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,real(tmparray,r4),dims,error)
    deallocate(tmparray)
   else
    call h5dwrite_f(dset_id,h5in%wrd_type,real(array,r4),dims,error)
   endif
  else
   if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
    deallocate(tmparray)
   else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
   endif
  endif
#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Writing data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  if(allocated(tmparray)) deallocate(tmparray)
  return
  end subroutine dump_h5_3d
!-----------------------------------------------------------------------
! subprogram 23. dump_h5_4d
!-----------------------------------------------------------------------
  subroutine dump_h5_4d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r8), dimension(:,:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error, i,j,k
  real(r8), dimension(:,:,:,:), allocatable :: tmparray

  integer(HID_T) dspace_id, rank, dset_id, plist_id
  integer(HSIZE_T) :: dims(4)
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 4
  if(h5in%doTranspose) then
   dims(4)=SIZE(array,1);  dims(2)=SIZE(array,3)
   dims(3)=SIZE(array,2);  dims(1)=SIZE(array,4)
   allocate(tmparray(dims(1),dims(2),dims(3),dims(4)))
   do i=1,dims(1); do j=1,dims(2); do k=1,dims(3)
      tmparray(i,j,k,:)=array(:,k,j,i)
  enddo; enddo; enddo
  else
   dims(1)=SIZE(array,1);  dims(2)=SIZE(array,2)
   dims(3)=SIZE(array,3);  dims(4)=SIZE(array,4)
  endif
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! reate property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
  if(h5in%typeConvert) then
    if(h5in%doTranspose) then
      call h5dwrite_f(dset_id,h5in%wrd_type,real(tmparray,r4),dims, &
                     error,xfer_prp = plist_id)
      deallocate(tmparray)
    else
      call h5dwrite_f(dset_id,h5in%wrd_type,real(array,r4),dims, &
                     error,xfer_prp = plist_id)
    endif
  else
    if(h5in%doTranspose) then
      call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                     error,xfer_prp = plist_id)
      deallocate(tmparray)
    else
      call h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                     error,xfer_prp = plist_id)
    endif
  endif
#else
  if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
    deallocate(tmparray)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
  endif
#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Writing data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  if(allocated(tmparray)) deallocate(tmparray)
  return
  end subroutine dump_h5_4d
!-----------------------------------------------------------------------
! subprogram 20. dump_rl_1d
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_rl_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r4), dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: error
  integer :: rank
  integer(HID_T) :: dspace_id, dset_id
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(1) :: dims
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 1;            dims(:) = (/SIZE(array,1)/)
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#if 0
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
!#ifdef __MPI
!   call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error, &
!                 xfer_prp = plist_id)
!#else
   call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
!#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Data set write failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_rl_1d
!-----------------------------------------------------------------------
! subprogram 21. dump_rl_2d
! Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
  subroutine dump_rl_2d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r4), dimension(:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error

  integer(HID_T) dspace_id, rank, dset_id, plist_id
  integer(HSIZE_T) :: dims(2)
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 2
  if(h5in%doTranspose) then
     dims(2) = SIZE(array,1);  dims(1) = SIZE(array,2)
  else
     dims(1) = SIZE(array,1);  dims(2) = SIZE(array,2)
  endif
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!    Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
  if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                   error,xfer_prp = plist_id)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                   error,xfer_prp = plist_id)
  endif
#else
  if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,TRANSPOSE(array),dims, &
                   error)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
  endif
#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Writing data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_rl_2d

!-----------------------------------------------------------------------
! subprogram 22. dump_rl_3d
!-----------------------------------------------------------------------
  subroutine dump_rl_3d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r4), dimension(:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error, i,j
  real(r8), dimension(:,:,:), allocatable :: tmparray

  integer(HID_T) dspace_id, rank, dset_id, plist_id
  integer(HSIZE_T) :: dims(3)
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 3
  if(h5in%doTranspose) then
   dims(3)=SIZE(array,1);dims(2)=SIZE(array,2);dims(1)=SIZE(array,3)
   allocate(tmparray(dims(1),dims(2),dims(3)))
   do i=1,dims(1); do j=1,dims(2)
      tmparray(i,j,:)=array(:,j,i)
   enddo; enddo
  else
   dims(1)=SIZE(array,1);dims(2)=SIZE(array,2);dims(3)=SIZE(array,3)
  endif
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
  if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                   error,xfer_prp = plist_id)
    deallocate(tmparray)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                   error,xfer_prp = plist_id)
  endif
#else
  if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
    deallocate(tmparray)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
  endif
#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Writing data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_rl_3d
!-----------------------------------------------------------------------
! subprogram 23. dump_rl_4d
!-----------------------------------------------------------------------
  subroutine dump_rl_4d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r4), dimension(:,:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  integer,parameter :: FAIL=-1
  integer :: error, i,j,k
  real(r8), dimension(:,:,:,:), allocatable :: tmparray

  integer(HID_T) dspace_id, rank, dset_id, plist_id
  integer(HSIZE_T) :: dims(4)
!-----------------------------------------------------------------------
! Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  rank = 4
  if(h5in%doTranspose) then
   dims(4)=SIZE(array,1);  dims(3)=SIZE(array,2)
   dims(3)=SIZE(array,2);  dims(1)=SIZE(array,4)
   allocate(tmparray(dims(1),dims(2),dims(3),dims(4)))
   do i=1,dims(1); do j=1,dims(2); do k=1,dims(3)
      tmparray(i,j,k,:)=array(:,k,j,i)
  enddo; enddo; enddo
  else
   dims(1)=SIZE(array,1);  dims(2)=SIZE(array,2)
   dims(3)=SIZE(array,3);  dims(4)=SIZE(array,4)
  endif
!-----------------------------------------------------------------------
! Create the data space.
!-----------------------------------------------------------------------
  call h5screate_simple_f(rank,dims,dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Create the data set.
! Note: wrd_type is data type being written into file (r4 or r8)
!-----------------------------------------------------------------------
  call h5dcreate_f(inid,aname,h5in%wrd_type,dspace_id,dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Create data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
!Create property list for collective dataset write
!-----------------------------------------------------------------------
#ifdef __MPI
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
   if (error==FAIL) then
      errval%errorMsg = 'ERROR: Creating plist failed for '//aname
      errval%errBool = .true.
      return
   endif
#endif
!-----------------------------------------------------------------------
! Write stored data to "name" data set.
!-----------------------------------------------------------------------
#ifdef __MPI
  if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims, &
                   error,xfer_prp = plist_id)
    deallocate(tmparray)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims, &
                   error,xfer_prp = plist_id)
  endif
#else
  if(h5in%doTranspose) then
    call h5dwrite_f(dset_id,h5in%wrd_type,tmparray,dims,error)
    deallocate(tmparray)
  else
    call h5dwrite_f(dset_id,h5in%wrd_type,array,dims,error)
  endif
#endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Writing data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Add the VisSchema attributes
!-----------------------------------------------------------------------
  call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine dump_rl_4d
!-----------------------------------------------------------------------
! subprogram 30. add_h5_int
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a 1D array
!-----------------------------------------------------------------------
  subroutine add_h5_int(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  integer(i4), intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: rank = 1
  integer :: error
  integer(HID_T) :: dspace_id, dset_id, filespace
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(1) :: dims,maxdims,chunk_dims,extdims,offset
  integer(HSIZE_T), dimension(1) :: olddims,oldmaxdims
  integer(HID_T) :: cparms        !dataset creatation property identifier 
  LOGICAL :: dset_exists
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! First determine whether the dataset exists because different
! whether data exists or not.
!-----------------------------------------------------------------------
  call h5lexists_f(inid, aname, dset_exists, error)

!-----------------------------------------------------------------------
! Do the case of creating the data
!-----------------------------------------------------------------------
  if (.NOT. dset_exists) then
    !
    ! Create the data space with unlimited dimensions.
    maxdims = (/H5S_UNLIMITED_f/)
    dims = (/1/)
    call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
     if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data space failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
    chunk_dims = (/1/)
    call h5pset_chunk_f(cparms, rank, chunk_dims, error)

    ! Create a new dataset within the file using cparms creation properties.
    call h5dcreate_f(inid,aname,H5T_NATIVE_INTEGER,dspace_id,dset_id,error,cparms)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data set failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Write stored data to "name" data set.
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,value,dims,error)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Data set write failed for '//aname
       errval%errBool = .true.
       return
    endif
    call h5pclose_f(cparms, error) !Close the property list.
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Close property list failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Add the VisSchema attributes on the first step
    call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Case for appending the dataset
!-----------------------------------------------------------------------
  else
    call h5dopen_f(inid,aname, dset_id, error)
    
    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
    
    ! Extend the dataset. This call assures that dataset has the space
    dims = (/1/)
    extdims=dims
    extdims(1)=olddims(1)+dims(1)
    call h5dextend_f(dset_id, extdims, error)

    ! Define memory space
    call h5screate_simple_f(rank, dims, dspace_id, error)

    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    offset=0
    offset(1) = olddims(1)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
    
    ! Write the data to the hyperslab.
    call H5Dwrite_f(dset_id,H5T_NATIVE_INTEGER,value,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    call h5sclose_f(filespace, error)
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine add_h5_int
!-----------------------------------------------------------------------
! subprogram 30. add_h5_dbl
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a 1D array
!-----------------------------------------------------------------------
  subroutine add_h5_dbl(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: rank = 1
  integer :: error
  integer(HID_T) :: dspace_id=0, dset_id=0, filespace=0
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(1) :: dims=0,maxdims,extdims,offset
  integer(HSIZE_T), dimension(1) :: olddims,oldmaxdims
  integer(HID_T) :: cparms        !dataset creatation property identifier 
  LOGICAL :: dset_exists
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! First determine whether the dataset exists because different
! whether data exists or not.
!-----------------------------------------------------------------------
  call h5lexists_f(inid, aname, dset_exists, error)

!-----------------------------------------------------------------------
! Do the case of creating the data
!-----------------------------------------------------------------------
  if (.NOT. dset_exists) then
    !
    ! Create the data space with unlimited dimensions.
    !
    maxdims = (/H5S_UNLIMITED_f/)
    dims = (/1/)
    call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
     if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data space failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
    call h5pset_chunk_f(cparms, rank, dims, error)

    ! Create a new dataset within the file using cparms creation properties.
    call h5dcreate_f(inid,aname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error,cparms)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data set failed for '//aname
       errval%errBool = .true.
       return
    endif
    call h5pclose_f(cparms, error) !Close the property list.
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Close property list failed for '//aname
       errval%errBool = .true.
       return
    endif

    ! Write stored data to "name" data set.
    !call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,value,dims,error)
    call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,value,dims,error)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Data set write failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Add the VisSchema attributes on the first time step
    call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Case for appending the dataset
!-----------------------------------------------------------------------
  else
    call h5dopen_f(inid,aname, dset_id, error)
    
    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
    
    ! Extend the dataset. This call assures that dataset has the space
    dims = (/1/)
    extdims=dims
    extdims(1)=olddims(1)+dims(1)
    call h5dextend_f(dset_id, extdims, error)

    ! Define memory space
    call h5screate_simple_f(rank, dims, dspace_id, error)

    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    offset=0
    offset(1) = olddims(1)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
    
    ! Write the data to the hyperslab.
    call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,value,dims,error, &
                    file_space_id=filespace,mem_space_id=dspace_id)
    call h5sclose_f(filespace, error)
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine add_h5_dbl
!-----------------------------------------------------------------------
! subprogram 30. add_h5_int_1d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_int_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  integer(i4), dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: rank = 2, asize
  integer :: error
  integer(HID_T) :: dspace_id, dset_id, filespace
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(2) :: dims,maxdims,chunk_dims,extdims,offset
  integer(HSIZE_T), dimension(2) :: olddims,oldmaxdims
  integer(HID_T) :: cparms        !dataset creatation property identifier 
  LOGICAL :: dset_exists
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! First determine whether the dataset exists because different
! whether data exists or not.
!-----------------------------------------------------------------------
  call h5lexists_f(inid, aname, dset_exists, error)
  asize=SIZE(array)

!-----------------------------------------------------------------------
! Do the case of creating the data
!-----------------------------------------------------------------------
  if (.NOT. dset_exists) then
    !
    ! Create the data space with unlimited dimensions.
    maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
    ! For convenience, put the time step (extendible set, as the first index
    dims = (/1, asize/)
    call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
     if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data space failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
    chunk_dims = dims
    call h5pset_chunk_f(cparms, rank, chunk_dims, error)

    ! Create a new dataset within the file using cparms creation properties.
    call h5dcreate_f(inid,aname,H5T_NATIVE_INTEGER,dspace_id,dset_id,error,cparms)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data set failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Write stored data to "name" data set.
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,array,dims,error)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Data set write failed for '//aname
       errval%errBool = .true.
       return
    endif
    call h5pclose_f(cparms, error) !Close the property list.
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Close property list failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Add the VisSchema attributes on the first step
    call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Case for appending the dataset
!-----------------------------------------------------------------------
  else
    call h5dopen_f(inid,aname, dset_id, error)
    
    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
    
    ! Extend the dataset. This call assures that dataset has the space
    dims = (/1, asize/)
    extdims=dims
    extdims(1)=olddims(1)+dims(1)
    extdims(2)=asize
    call h5dextend_f(dset_id, extdims, error)

    ! Define memory space
    call h5screate_simple_f(rank, dims, dspace_id, error)

    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    offset=0
    offset(1) = olddims(1)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
    
    ! Write the data to the hyperslab.
    call H5Dwrite_f(dset_id,H5T_NATIVE_INTEGER,array,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    call h5sclose_f(filespace, error)
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine add_h5_int_1d
!-----------------------------------------------------------------------
! subprogram 30. add_h5_1d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: rank = 2, asize
  integer :: error
  integer(HID_T) :: dspace_id, dset_id, filespace
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HSIZE_T), dimension(2) :: dims,maxdims,chunk_dims,extdims,offset
  integer(HSIZE_T), dimension(2) :: olddims,oldmaxdims
  integer(HID_T) :: cparms        !dataset creatation property identifier 
  LOGICAL :: dset_exists
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! First determine whether the dataset exists because different
! whether data exists or not.
!-----------------------------------------------------------------------
  call h5lexists_f(inid, aname, dset_exists, error)
  asize=SIZE(array)

!-----------------------------------------------------------------------
! Do the case of creating the data
!-----------------------------------------------------------------------
  if (.NOT. dset_exists) then
    !
    ! Create the data space with unlimited dimensions.
    maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
    ! For convenience, put the time step (extendible set, as the first index
    dims = (/1, asize/)
    call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
     if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data space failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
    chunk_dims = dims
    call h5pset_chunk_f(cparms, rank, chunk_dims, error)

    ! Create a new dataset within the file using cparms creation properties.
    call h5dcreate_f(inid,aname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error,cparms)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data set failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Write stored data to "name" data set.
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Data set write failed for '//aname
       errval%errBool = .true.
       return
    endif
    call h5pclose_f(cparms, error) !Close the property list.
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Close property list failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Add the VisSchema attributes on the first step
    call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Case for appending the dataset
!-----------------------------------------------------------------------
  else
    call h5dopen_f(inid,aname, dset_id, error)
    
    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
    
    ! Extend the dataset. This call assures that dataset has the space
    dims = (/1, asize/)
    extdims=dims
    extdims(1)=olddims(1)+dims(1)
    extdims(2)=asize
    call h5dextend_f(dset_id, extdims, error)

    ! Define memory space
    call h5screate_simple_f(rank, dims, dspace_id, error)

    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    offset=0
    offset(1) = olddims(1)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
    
    ! Write the data to the hyperslab.
    call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    call h5sclose_f(filespace, error)
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine add_h5_1d
!-----------------------------------------------------------------------
! subprogram 30. add_h5_2d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_2d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:,:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: rank = 3
  integer :: error
  integer(HID_T) :: dspace_id, dset_id, filespace
  integer(HID_T) :: plist_id       ! Property list identifier
  integer, dimension(2)          :: asize
  integer(HSIZE_T), dimension(3) :: dims,maxdims,chunk_dims,extdims,offset
  integer(HSIZE_T), dimension(3) :: olddims,oldmaxdims
  integer(HID_T) :: cparms        !dataset creatation property identifier 
  LOGICAL :: dset_exists
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! First determine whether the dataset exists because different
! whether data exists or not.
!-----------------------------------------------------------------------
  call h5lexists_f(inid, aname, dset_exists, error)
  asize(1)=SIZE(array,1)
  asize(2)=SIZE(array,2)
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname

  if(h5in%doTranspose) then
    dims = (/1, asize(2), asize(1)/)
  else
    dims = (/1, asize(1), asize(2)/)
  endif

!-----------------------------------------------------------------------
! Do the case of creating the data
!-----------------------------------------------------------------------
  if (.NOT. dset_exists) then
    !
    ! Create the data space with unlimited dimensions.
    maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
    ! For convenience, put the time step (extendible set, as the first index
    call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
     if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data space failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
    chunk_dims = dims
    call h5pset_chunk_f(cparms, rank, chunk_dims, error)

    ! Create a new dataset within the file using cparms creation properties.
    call h5dcreate_f(inid,aname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error,cparms)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data set failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Write stored data to "name" data set.
    if(h5in%doTranspose) then
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,TRANSPOSE(array),dims,error)
    else
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error)
    endif

    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Data set write failed for '//aname
       errval%errBool = .true.
       return
    endif
    call h5pclose_f(cparms, error) !Close the property list.
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Close property list failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Add the VisSchema attributes on the first step
    call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Case for appending the dataset
!-----------------------------------------------------------------------
  else
    call h5dopen_f(inid,aname, dset_id, error)
    
    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
    
    ! Extend the dataset. This call assures that dataset has the space
    extdims=dims
    extdims(1)=olddims(1)+dims(1)
    !extdims(2)=asize(1)
    extdims(2)=dims(2)
    !extdims(3)=asize(2)
    extdims(3)=dims(3)

    call h5dextend_f(dset_id, extdims, error)

    ! Define memory space
    call h5screate_simple_f(rank, dims, dspace_id, error)

    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    offset=0
    offset(1) = olddims(1)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
    
    ! Write the data to the hyperslab.
    if(h5in%doTranspose) then
      call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,TRANSPOSE(array),dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    else
      call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    endif

    call h5sclose_f(filespace, error)
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine add_h5_2d

!-----------------------------------------------------------------------
! subprogram 30a. add_h5_3d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_3d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:,:,:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: rank = 4
  integer :: error
  integer(HID_T) :: dspace_id, dset_id, filespace
  integer(HID_T) :: plist_id       ! Property list identifier
  integer, dimension(3)          :: asize
  integer(HSIZE_T), dimension(4) :: dims,maxdims,chunk_dims,extdims,offset
  integer(HSIZE_T), dimension(4) :: olddims,oldmaxdims
  integer(HID_T) :: cparms        !dataset creatation property identifier 
  LOGICAL :: dset_exists
  integer ::i,j
  real(r8), dimension(:,:,:), allocatable :: tmparray
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! First determine whether the dataset exists because different
! whether data exists or not.
!-----------------------------------------------------------------------
  call h5lexists_f(inid, aname, dset_exists, error)
  asize(1)=SIZE(array,1)
  asize(2)=SIZE(array,2)
  asize(3)=SIZE(array,3)
  ! For convenience, put the time step (extendible set, as the first index
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  if(h5in%doTranspose) then
    dims = (/1, asize(3), asize(2),asize(1)/)
    allocate(tmparray(dims(2),dims(3),dims(4)))
    do i=1,dims(2); do j=1,dims(3)
      tmparray(i,j,:)=array(:,j,i)
    enddo; enddo
  else
    dims = (/1, asize(1), asize(2),asize(3)/)
  endif

!-----------------------------------------------------------------------
! Do the case of creating the data
!-----------------------------------------------------------------------
  if (.NOT. dset_exists) then
    !
    ! Create the data space with unlimited dimensions.
    maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f/)

    call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
     if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data space failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
    chunk_dims = dims
    call h5pset_chunk_f(cparms, rank, chunk_dims, error)

    ! Create a new dataset within the file using cparms creation properties.
    call h5dcreate_f(inid,aname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error,cparms)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data set failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Write stored data to "name" data set.
    if(h5in%doTranspose) then
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,tmparray,dims,error)
    else 
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error)
    endif

    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Data set write failed for '//aname
       errval%errBool = .true.
       return
    endif
    call h5pclose_f(cparms, error) !Close the property list.
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Close property list failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Add the VisSchema attributes on the first step
    call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Case for appending the dataset
!-----------------------------------------------------------------------
  else
    call h5dopen_f(inid,aname, dset_id, error)
    
    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
    
    ! Extend the dataset. This call assures that dataset has the space
    extdims=dims
    extdims(1)=olddims(1)+dims(1)
!    extdims(2)=asize(1)
!    extdims(3)=asize(2)
!    extdims(4)=asize(3)
    extdims(2)=dims(2)
    extdims(3)=dims(3)
    extdims(4)=dims(4)

    call h5dextend_f(dset_id, extdims, error)

    ! Define memory space
    call h5screate_simple_f(rank, dims, dspace_id, error)

    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    offset=0
    offset(1) = olddims(1)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
    
    ! Write the data to the hyperslab.
    if(h5in%doTranspose) then
      call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,tmparray,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    else 
      call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    endif
    call h5sclose_f(filespace, error)
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  if(allocated(tmparray)) deallocate(tmparray)
  return
  end subroutine add_h5_3d
!-----------------------------------------------------------------------
! subprogram 30a. add_h5_4d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_4d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:,:,:,:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer,parameter :: FAIL=-1
  integer :: rank = 5
  integer :: error
  integer(HID_T) :: dspace_id, dset_id, filespace
  integer(HID_T) :: plist_id       ! Property list identifier
  integer, dimension(4)          :: asize
  integer(HSIZE_T), dimension(5) :: dims,maxdims,chunk_dims,extdims,offset
  integer(HSIZE_T), dimension(5) :: olddims,oldmaxdims
  integer(HID_T) :: cparms        !dataset creatation property identifier 
  LOGICAL :: dset_exists
  integer :: i,j,k
  real(r8), dimension(:,:,:,:), allocatable :: tmparray
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
! First determine whether the dataset exists because different
! whether data exists or not.
!-----------------------------------------------------------------------
  call h5lexists_f(inid, aname, dset_exists, error)
  asize(1)=SIZE(array,1)
  asize(2)=SIZE(array,2)
  asize(3)=SIZE(array,3)
  asize(4)=SIZE(array,4)

  if(h5in%verbose) WRITE(*,*) 'Writing ', aname
  if(h5in%doTranspose) then
   dims = (/1, asize(4), asize(3),asize(2),asize(1)/)
   allocate(tmparray(dims(2),dims(3),dims(4),dims(5)))
   do i=1,dims(2); do j=1,dims(3); do k=1,dims(4)
      tmparray(i,j,k,:)=array(:,k,j,i)
  enddo; enddo; enddo
  else
   dims = (/1, asize(1), asize(2),asize(3),asize(4)/)
  endif


!-----------------------------------------------------------------------
! Do the case of creating the data
!-----------------------------------------------------------------------
  if (.NOT. dset_exists) then
    !
    ! Create the data space with unlimited dimensions.
    maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
    ! For convenience, put the time step (extendible set, as the first index
    call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
     if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data space failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
    chunk_dims = dims
    call h5pset_chunk_f(cparms, rank, chunk_dims, error)

    ! Create a new dataset within the file using cparms creation properties.
    call h5dcreate_f(inid,aname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error,cparms)
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Create data set failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Write stored data to "name" data set.
    if(h5in%doTranspose) then
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,tmparray,dims,error)
    else
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error)
    endif
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Data set write failed for '//aname
       errval%errBool = .true.
       return
    endif
    call h5pclose_f(cparms, error) !Close the property list.
    if (error==FAIL) then
       errval%errorMsg = 'ERROR: Close property list failed for '//aname
       errval%errBool = .true.
       return
    endif
    ! Add the VisSchema attributes on the first step
    call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
! Case for appending the dataset
!-----------------------------------------------------------------------
  else
    call h5dopen_f(inid,aname, dset_id, error)
    
    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
    
    ! Extend the dataset. This call assures that dataset has the space
    extdims=dims
    extdims(1)=olddims(1)+dims(1)
!     extdims(2)=asize(1)
!    extdims(3)=asize(2)
!    extdims(4)=asize(3)
!    extdims(5)=asize(4)
    extdims(2)=dims(2)
    extdims(3)=dims(3)
    extdims(4)=dims(4)
    extdims(5)=dims(5)
    call h5dextend_f(dset_id, extdims, error)

    ! Define memory space
    call h5screate_simple_f(rank, dims, dspace_id, error)

    ! Open filespace in existing dataset and get existing shape and size of data
    call h5dget_space_f(dset_id, filespace, error)

    offset=0
    offset(1) = olddims(1)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
    
    ! Write the data to the hyperslab.
    if(h5in%doTranspose) then
      call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,tmparray,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    else
      call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,error,  &
                    file_space_id=filespace,mem_space_id=dspace_id)
    endif
    call h5sclose_f(filespace, error)
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
  call h5sclose_f(dspace_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data space failed for '//aname
     errval%errBool = .true.
     return
  endif
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  if(allocated(tmparray)) deallocate(tmparray)
  return
  end subroutine add_h5_4d
!-----------------------------------------------------------------------
! subprogram 40. read_dims
! Read the dimensions of dataset associated with aname: 1d array
!-----------------------------------------------------------------------
  subroutine read_dims(dset_id,aname,dims,errval)
  integer(HID_T), intent(in) :: dset_id
  character(*), intent(in) :: aname
  integer(HSIZE_T), dimension(:), intent(inoUT) :: dims
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer,parameter :: FAIL=-1
  integer :: error

  integer(HSIZE_T), dimension(:), allocatable :: maxdims
  integer(HID_T) dspace_id
!-----------------------------------------------------------------------
! Get dataset's dataspace handle.
!-----------------------------------------------------------------------
  allocate(maxdims(SIZE(dims)))
  call h5dget_space_f(dset_id, dspace_id, error)
!-----------------------------------------------------------------------
! Get dataspace's dimensinons.
!-----------------------------------------------------------------------
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine read_dims

!-----------------------------------------------------------------------
! subprogram 41. read_h5_1d
! Read simple data set: 1d array
!-----------------------------------------------------------------------
  subroutine read_h5_1d(fid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  TYPE(hdf5InOpts), intent(in) :: h5in
  real(r8), dimension(:), intent(inoUT) :: array
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer(HSIZE_T), dimension(1) :: dims, fdims
  integer,parameter :: FAIL=-1
  integer :: error

  integer(HID_T) dset_id
!-----------------------------------------------------------------------
! Open the dataset specified by aname
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) ' Reading 1d array: ', aname
  call h5dopen_f(fid, aname, dset_id, error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Find data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Read data set
!-----------------------------------------------------------------------
  dims(1)=SIZE(array)
  call read_dims(dset_id,aname,fdims,errval)
  call check_dims(dims,fdims, errval)
  if (errval%errBool) return
  call h5dread_f(dset_id,h5in%wrd_type,array,dims,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Reading data set failed for '//aname
     errval%errBool = .true.
     call h5dclose_f(dset_id,error)
     return
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine read_h5_1d

!-----------------------------------------------------------------------
! subprogram 42. read_h5_2d
! Read simple data set: 2d array
!-----------------------------------------------------------------------
  subroutine read_h5_2d(fid,aname,array,h5in,errval)
  integer, intent(in) :: fid
  character*(*), intent(in) :: aname
  real(r8), dimension(:,:), intent(inoUT) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer(HSIZE_T), dimension(2) :: dims, fdims
  integer,parameter :: FAIL=-1
  integer :: error, i
  real(r8), dimension(:,:), allocatable :: tmparray

  integer(HID_T) dset_id
!-----------------------------------------------------------------------
! Open the dataset specified by aname
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) ' Reading 2d array: ', aname
  call h5dopen_f(fid, aname, dset_id, error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Find data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Check dims
!-----------------------------------------------------------------------
  if(.NOT. h5in%doTranspose) then
     dims(1)=SIZE(array,1); dims(2)=SIZE(array,2)
  else
     dims(1)=SIZE(array,2); dims(2)=SIZE(array,1)
  endif
  call read_dims(dset_id,aname,fdims,errval)
  call check_dims(dims,fdims, errval)
  if (errval%errBool) then
    call h5dclose_f(dset_id,error)
    return
  endif
!-----------------------------------------------------------------------
! Read data set
!-----------------------------------------------------------------------
  if(.NOT. h5in%doTranspose) then
     call h5dread_f(dset_id,h5in%wrd_type,array,dims,error)
  else
     allocate(tmparray(dims(1),dims(2)))
     call h5dread_f(dset_id,h5in%wrd_type,tmparray,dims,error)
     do i=1,dims(1);      array(:,i)=tmparray(i,:);     enddo
     deallocate(tmparray)
  endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Read data set failed for '//aname
     errval%errBool = .true.
     call h5dclose_f(dset_id,error)
     return
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine read_h5_2d

!-----------------------------------------------------------------------
! subprogram 43. read_h5_3d
! Read simple data set: 3d array
!-----------------------------------------------------------------------
  subroutine read_h5_3d(fid,aname,array,h5in,errval)
  integer, intent(in) :: fid
  character*(*), intent(in) :: aname
  real(r8), dimension(:,:,:), intent(inoUT) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer(HSIZE_T), dimension(3) :: dims, fdims
  integer,parameter :: FAIL=-1
  integer :: error, i,j
  real(r8), dimension(:,:,:), allocatable :: tmparray

  integer dset_id
!-----------------------------------------------------------------------
! Open the dataset specified by aname
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) ' Reading 3d array: ', aname
  call h5dopen_f(fid, TRIM(aname), dset_id, error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Find data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Check dims
!-----------------------------------------------------------------------
  if(.NOT. h5in%doTranspose) then
     dims(1)=SIZE(array,1)
     dims(2)=SIZE(array,2)
     dims(3)=SIZE(array,3)
  else
     dims(1)=SIZE(array,3)
     dims(2)=SIZE(array,2)
     dims(3)=SIZE(array,1)
  endif
  call read_dims(dset_id,aname,fdims,errval)
  !call check_dims(dims,fdims, errval)
  if (errval%errBool) then
    call h5dclose_f(dset_id,error)
    return
  endif
!-----------------------------------------------------------------------
! Read data set
!-----------------------------------------------------------------------
  if(.NOT. h5in%doTranspose) then
     call h5dread_f(dset_id,h5in%wrd_type,array,dims,error)
  else
     allocate(tmparray(dims(1),dims(2),dims(3)))
     call h5dread_f(dset_id,h5in%wrd_type,tmparray,dims,error)
     do i=1,dims(1); do j=1,dims(2)
      array(:,j,i)=tmparray(i,j,:)
     enddo; enddo
     deallocate(tmparray)
  endif
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Read data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset
!-----------------------------------------------------------------------
  call h5dclose_f(dset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine read_h5_3d
!-----------------------------------------------------------------------
! subprogram 41. read_attribute_rl_sc
! Read real scalar attribute
!-----------------------------------------------------------------------
  subroutine read_attribute_rl_sc(fid,aname,val,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  TYPE(hdf5InOpts), intent(in) :: h5in
  double precision, intent(inoUT) :: val
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  integer(HSIZE_T), dimension(1) :: dims, fdims
  integer,parameter :: FAIL=-1
  integer :: error

  integer(HID_T) aset_id
!-----------------------------------------------------------------------
! Open the dataset specified by aname
!-----------------------------------------------------------------------
  if(h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
  call h5aopen_name_f(fid, aname, aset_id, error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Find data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
! Read attribute
!-----------------------------------------------------------------------
  dims(1)=1
  call h5aread_f(aset_id,h5in%wrd_type,val,dims,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Reading data set failed for '//aname
     errval%errBool = .true.
     call h5aclose_f(aset_id,error)
     return
  endif
!-----------------------------------------------------------------------
! Terminate access to the dataset
!-----------------------------------------------------------------------
  call h5aclose_f(aset_id,error)
  if (error==FAIL) then
     errval%errorMsg = 'ERROR: Close data set failed for '//aname
     errval%errBool = .true.
     return
  endif
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine read_attribute_rl_sc
!-----------------------------------------------------------------------
! close module.
!-------------------------------------------------------------------
  END MODULE hdf5_api
