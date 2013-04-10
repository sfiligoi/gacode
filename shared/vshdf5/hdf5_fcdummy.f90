!-----------------------------------------------------------------------
! dummy routines to resolve the linker symbols
!------------------------------------------------------------------
! module hdf_api
!-----------------------------------------------------------------------
  module hdf5_api
  implicit none
  character(5), parameter, private :: h5fortranapiversion="1.0"
  integer, parameter, private :: i4=selected_int_kind(9)
  integer, parameter, private :: i8=selected_int_kind(18)
  integer, parameter, private :: r4=selected_real_kind(6,37)
  integer, parameter, private :: r8=selected_real_kind(13,307)
  integer,parameter :: HID_T=i4
  integer,parameter :: HIDT=i4
  integer,parameter :: HSIZET=i4
  integer,parameter :: HSIZE_T=i4
  integer,parameter :: SIZET=i4
  integer,parameter :: SIZE_T=i4
!-----------------------------------------------------------------------
  type hdf5inopts
     !integer(hid_t) :: wrd_type
     integer :: write_kind_real
     integer :: write_kind_int
     logical :: dotranspose        
     logical :: verbose               ! whether to write verbose output
     logical :: debug                 ! write even more verbose output for debugging
     logical :: pio                   ! whether parallel i/o is used
     logical :: wrvstime              ! whether to write time to attribute
     logical :: typeconvert           ! whether to demote the type
     logical :: unitconvert           ! whether to write vsunitcnv
     integer(selected_int_kind(9)) :: comm     ! Communicator associated w/ open file
     integer(selected_int_kind(9)) :: info     ! mpi_info: set to MPI_INFO_NULL if unsure
     character(len=30) :: mesh             ! See above
     character(len=30) :: units            ! Units
     character(len=30) :: vsAxisLabels     ! Axis labels
     character(len=30) :: vsCentering      ! How to center variables on mesh
     character(len=30) :: vsMD             ! Multidomain variable
     character(len=30) :: vsTimeGroup      ! Time group label
     character(len=30) :: vsIndexOrder     ! Data ordering
     character(len=1000) :: vsLabels        ! Labels for (ic) (nqty)
     real(r8) :: vsTime                    ! Time
     integer(HID_T) :: vsStep=0            ! Step # associated with time
     real(r8) :: vsUnitCnv=1.              ! conversion factor
     character(len=30), dimension(3) :: vsAxis   ! For rectilinear meshes
     INTEGER(HID_T) :: h5_kind_type_r4
     INTEGER(HID_T) :: h5_kind_type_r8
     INTEGER(HID_T) :: h5_kind_type_i4
     INTEGER(HID_T) :: h5_kind_type_i8
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
                   ,write_attribute_rl_sc,write_attribute_rls_sc &
                   ,write_attribute_rl_vec &
                   ,write_attribute_int_sc,write_attribute_int_vec &
                   ,write_attribute_intl_sc,write_attribute_intl_vec
  end interface
  interface dump_h5
    module procedure dump_int,dump_intl,dump_int_1d,dump_intl_1d, &
                     dump_h5_rldbl,dump_h5_rls, &
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
    module procedure read_h5_intl,read_h5_intl_1d,read_h5_1d,read_h5_2d,read_h5_3d
  end interface
  interface read_attribute
    module procedure read_attribute_intl_sc,read_attribute_rl_sc,read_attribute_intl_vec
  end interface


  contains


!-----------------------------------------------------------------------
! Commonly used hdf5 routines
!-----------------------------------------------------------------------
  subroutine hdf5_close()
  return
  end subroutine hdf5_close

!-----------------------------------------------------------------------
! subprogram 0. check_dims 
! Write mismatched dimension errors
!-----------------------------------------------------------------------
  subroutine check_dims(dims, fdims, errval)
  integer(hsize_t), dimension(:), intent(in) :: dims, fdims
  type(hdf5errortype), intent(inout) :: errval
  return
  end subroutine check_dims

!-----------------------------------------------------------------------
! function 1. h5accessMethod
! Return something the correct hdf5 access method given a more
! memoral names
!-----------------------------------------------------------------------
  function h5accessmethod(access_method)
  integer(hid_t) :: h5accessmethod
  character(*), intent(in) :: access_method
  return
  end function h5accessmethod

!-----------------------------------------------------------------------
! subprogram 0. vshdf5_fcinit
! Open fortran hdf5 and set open/close parameters.
!-----------------------------------------------------------------------
  subroutine vshdf5_fcinit()
  return
  end subroutine vshdf5_fcinit

!-----------------------------------------------------------------------
! subprogram 3 vshdf5_inith5vars
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
  h5in%vsCentering=" "
  h5in%doTranspose=.false.
  h5in%verbose=.false.
  h5in%debug=.false.
  h5in%pIO=.false.
  h5in%wrVsTime=.false.
  h5in%unitConvert=.false.
  h5in%mesh =  " "
  h5in%vsAxisLabels = " " 
  h5in%units =  " "
  h5in%vsCentering =  " "
  h5in%vsMD =  " "
  ! Data-ordering options:
  ! [ix][iy][iz][ic] compMinorC [iz][iy][ix][ic] compMinorF
  ! [ic][ix][iy][iz] compMajorC [ic][iz][iy][ix] compMajorF
  h5in%vsIndexOrder = " "
  h5in%vsLabels = " "
  h5err%errBool = .false.
  h5err%errorMsg =  " "

  h5in%h5_kind_type_r4 = r4
  h5in%h5_kind_type_r8 = r8
  h5in%h5_kind_type_i4 = i4
  h5in%h5_kind_type_i8 = i8

  h5in%write_kind_real=h5in%h5_kind_type_r8
  h5in%write_kind_int=h5in%h5_kind_type_i4

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
  return
  end subroutine open_h5file
!-----------------------------------------------------------------------
! subprogram 5. open_oldh5file
! Open file for writing and write file attributes
! Create the group for the independent variables at this stage
!-----------------------------------------------------------------------
  subroutine open_oldh5file(fname,fileid,rootgid,h5in,h5err)
  character(*), intent(in) :: fname
  integer(hid_t), intent(out) :: fileid,rootgid
  type(hdf5errortype), intent(inout) :: h5err
  type(hdf5inopts), intent(inout) :: h5in
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
  return
  end subroutine open_newh5file

!-----------------------------------------------------------------------
! subprogram 7. close_h5file
!     Close the file associated with fileId.
!-----------------------------------------------------------------------
  subroutine close_h5file(fileId,root_id,h5err)
  integer(HID_T), intent(in)         :: fileId
  integer(HID_T), intent(in)                :: root_id
  TYPE(hdf5ErrorType), intent(inoUT) :: h5err
  return
  end subroutine close_h5file

!-----------------------------------------------------------------------
! subprogram 8. open_group
!     Open a group in a safe way
!-----------------------------------------------------------------------
  subroutine open_group(inid,gname,gid,errval)
  character(*), intent(in) :: gname
  integer(HID_T), intent(in) :: inid
  integer(HID_T), intent(out) :: gid
  type(hdf5ErrorType), intent(inout) :: errval
  return
  end subroutine open_group

!-----------------------------------------------------------------------
! subprogram 9. make_group
! Create a group in a safe way
!-----------------------------------------------------------------------
  subroutine make_group(inid,gname,gid,h5in,errval)
  character(*), intent(in) :: gname
  integer(HID_T), intent(in) :: inid
  integer(HID_T), intent(out) :: gid
  type(hdf5InOpts), intent(in) :: h5in
  type(hdf5ErrorType), intent(inout) :: errval
  return
  end subroutine make_group

!-----------------------------------------------------------------------
! subprogram 10. close_group
! Close a group in a safe way
!-----------------------------------------------------------------------
  subroutine close_group(gname,inid,errval)
  character(*), intent(in) :: gname
  integer(HID_T), intent(in) :: inid
  type(hdf5ErrorType), intent(inout) :: errval
  errval%errBool = .false.
  return
  end subroutine close_group

!-----------------------------------------------------------------------
! subprogram 11. test_group
! See if group exists
!-----------------------------------------------------------------------
  subroutine test_group(inid,gname,group_exists,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: gname
  LOGICAL, intent(out) :: group_exists
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  return
  end subroutine test_group

!-----------------------------------------------------------------------
! subprogram 12. get_nmember
! Get the number of members of a group
!-----------------------------------------------------------------------
  subroutine get_nmembers(inid,gname,nmembers,errval)
  character(*), intent(in) :: gname
  integer(HID_T), intent(in) :: inid
  integer(HID_T), intent(out) :: nmembers
  type(hdf5errortype), intent(inout) :: errval
!-----------------------------------------------------------------------
  errval%errBool = .false.
  return
  end subroutine get_nmembers

!-----------------------------------------------------------------------
! subprogram 13. make reference
! Simplify referencing of one object (source) to another (target)
!-----------------------------------------------------------------------
  subroutine make_reference(inid,tname,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: tname
  type(hdf5errortype), intent(inout) :: errval
  errval%errBool = .false.
  return
  end subroutine make_reference
!-----------------------------------------------------------------------
! subprogram 20. make_mesh_group
! Defines a mesh group that points to other variables that define
!  the actual mesh
!-----------------------------------------------------------------------
!  subroutine make_mesh_group(gInId,gridId,h5in,meshName,&
!              meshKind,axis0,axis1,axis2,transform,trName,errval)
  subroutine make_mesh_group(gInId,meshName,gridId,h5in,&
              errval)
  integer(HID_T), intent(in) :: gInId
  character*(*), intent(in) :: meshName!,axis0,axis1,axis2
  integer(HID_T), intent(inoUT) :: gridId
  TYPE(hdf5InOpts), intent(inout) :: h5in
  !character*(*), intent(in) :: meshKind,transform,trName
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
!-----------------------------------------------------------------------
  return
  end subroutine make_mesh_group
!-----------------------------------------------------------------------
! subprogram 21. make_time_group
! Make a group that contains the time data.  See:
!    https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
  subroutine make_time_group(gInId,h5in,h5err)
  integer(HID_T), intent(in) :: gInId
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: h5err
  return
  end subroutine make_time_group

!-----------------------------------------------------------------------
! subprogram 22. make_vec_group
! Make a group that defines a vector.  See:
!    https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
  subroutine make_vec_group(gInId,grName,veclabel,h5in,h5err)
  integer(HID_T), INTENT(IN) :: gInId
  character*(*), INTENT(IN) :: grName,veclabel
  type(hdf5InOpts), INTENT(IN) :: h5in
  type(hdf5ErrorType), INTENT(INOUT) :: h5err
  return
  end subroutine make_vec_group
!-----------------------------------------------------------------------
! subprogram 23. make_limits_group
! Make a group that contains the visualization region data.  See:
!    https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
  SUBROUTINE make_limits_group(gInId,grName,vsKind,lowerBound,      &
    upperBound,h5in,h5err)
  integer(HID_T), intent(in) :: gInId
  character*(*), intent(in) :: grName,vsKind
  real(r8), dimension(:), intent(in) :: lowerBound, upperBound      
  type(hdf5InOpts), intent(in) :: h5in
  type(hdf5ErrorType), intent(inoUT) :: h5err
  return
  END SUbroutine make_limits_group
!-----------------------------------------------------------------------
! subprogram 23b. make_external_link
! 
! Creates an external link, a soft link to an object in a different file
!-----------------------------------------------------------------------

  SUBROUTINE make_external_link(file_name, obj_name, link_loc_id, link_name, &
             h5in,h5err)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: file_name  
                       ! Name of the file containing the target object. Neither 
                       ! the file nor the target object is required to exist. 
                       ! May be the file the link is being created in.
  CHARACTER(LEN=*), INTENT(IN) :: obj_name  
                       ! Name of the target object, which need not already exist.
  INTEGER(HID_T), INTENT(IN) :: link_loc_id 
                       ! The file or group identifier for the new link.
  CHARACTER(LEN=*), INTENT(IN) :: link_name 
                       ! The name of the new link.
  INTEGER :: hdferr        
                       ! Error code: 
                       ! 0 on success and -1 on failure
  type(hdf5InOpts), intent(in) :: h5in
  type(hdf5ErrorType), intent(inoUT) :: h5err

  return
  END SUBROUTINE make_external_link
!-----------------------------------------------------------------------
! subprogram 14. write_attribute_ch_sc
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_ch_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname,attribute
  TYPE(hdf5ErrorType) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_ch_sc

!-----------------------------------------------------------------------
! subprogram 15. write_attribute_ch_vec
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_ch_vec(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  character*(*), dimension(:), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_ch_vec

!-----------------------------------------------------------------------
! subprogram 16. write_attribute_int_sc
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_int_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i4), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_int_sc

!-----------------------------------------------------------------------
! subprogram 17. write_attribute_int_vec
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_int_vec(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i4), dimension(:), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_int_vec

!-----------------------------------------------------------------------
! subprogram 12. write_attribute_intl_sc
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_intl_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i8), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_intl_sc

!-----------------------------------------------------------------------
! subprogram 13. write_attribute_int_vec
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_intl_vec(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i8), dimension(:), intent(in) :: attribute
  TYPE(hdf5ErrorType) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_intl_vec

!-----------------------------------------------------------------------
! subprogram 18. write_attribute_rl_sc
!-----------------------------------------------------------------------
  subroutine write_attribute_rl_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), intent(in) :: attribute
  TYPE(hdf5ErrorType), intent(out) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_rl_sc

!-----------------------------------------------------------------------
! subprogram 13. write_attribute_rls_sc
!-----------------------------------------------------------------------
  subroutine write_attribute_rls_sc(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r4), intent(in) :: attribute
  TYPE(hdf5ErrorType), intent(out) :: errval

  errval%errBool = .false.
  return
  end subroutine write_attribute_rls_sc

!-----------------------------------------------------------------------
! subprogram 19. write_attribute_rl_vec
! Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
  subroutine write_attribute_rl_vec(inid,aname,attribute,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), dimension(:), intent(in) :: attribute
  type(hdf5ErrorType) :: errval
  errval%errBool = .false.
  return
  end subroutine write_attribute_rl_vec

!-----------------------------------------------------------------------
! subprogram 25. dump_h5in_attributes
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_h5in_attributes(dset_id,h5in,h5err)
  integer(HID_T), intent(in) :: dset_id
  type(hdf5InOpts), intent(in) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  return
  end subroutine dump_h5in_attributes
!-----------------------------------------------------------------------
! subprogram 20. dump_int
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_int(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i4), intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_int
!-----------------------------------------------------------------------
! subprogram 20. dump_int_1d
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_int_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i4), dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_int_1d
!-----------------------------------------------------------------------
! subprogram 20. dump_intl
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_intl(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i8), intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_intl
!-----------------------------------------------------------------------
! subprogram 20. dump_intl_1d
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_intl_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  integer(i8), dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_intl_1d

!-----------------------------------------------------------------------
! subprogram 27. dump_h5_rldbl
! Write an hdf5 array + references to independent vars
!-------------------------------------------------------------------
  subroutine dump_h5_rldbl(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_h5_rldbl

!-------------------------------------------------------------------
! subprogram 20. dump_h5_rls
! Write an hdf5 array + references to independent vars
!-------------------------------------------------------------------
  subroutine dump_h5_rls(inid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r4), intent(in) :: value
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_h5_rls

!-----------------------------------------------------------------------
! subprogram 28. dump_h5_1d
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_h5_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r8), dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_h5_1d
!-----------------------------------------------------------------------
! subprogram 29. dump_h5_2d
! Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
  subroutine dump_h5_2d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r8), dimension(:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  errval%errBool = .false.
  return
  end subroutine dump_h5_2d

!-----------------------------------------------------------------------
! subprogram 30. dump_h5_3d
!-----------------------------------------------------------------------
  subroutine dump_h5_3d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r8), dimension(:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  type(hdf5ErrorType), intent(inout) :: errval
  return
  end subroutine dump_h5_3d
!-----------------------------------------------------------------------
! subprogram 31. dump_h5_4d
!-----------------------------------------------------------------------
  subroutine dump_h5_4d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r8), dimension(:,:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType) :: errval
  errval%errBool = .false.
  return
  end subroutine dump_h5_4d
!-----------------------------------------------------------------------
! subprogram 32. dump_rl_1d
! Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
  subroutine dump_rl_1d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character*(*), intent(in) :: aname
  real(r4), dimension(:), intent(in) :: array
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine dump_rl_1d
!-----------------------------------------------------------------------
! subprogram 33. dump_rl_2d
! Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
  subroutine dump_rl_2d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r4), dimension(:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  type(hdf5errortype), intent(inout) :: errval
  errval%errBool = .false.
  return
  end subroutine dump_rl_2d

!-----------------------------------------------------------------------
! subprogram 34. dump_rl_3d
!-----------------------------------------------------------------------
  subroutine dump_rl_3d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r4), dimension(:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  type(hdf5errortype), intent(inout) :: errval
  errval%errBool = .false.
  return
  end subroutine dump_rl_3d
!-----------------------------------------------------------------------
! subprogram 35. dump_rl_4d
!-----------------------------------------------------------------------
  subroutine dump_rl_4d(inid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  real(r4), dimension(:,:,:,:), intent(in) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  type(hdf5errortype), intent(inout) :: errval
  errval%errBool = .false.
  return
  end subroutine dump_rl_4d
!-----------------------------------------------------------------------
! subprogram 30. add_h5_int
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a 1D array
!-----------------------------------------------------------------------
  subroutine add_h5_int(inid,aname,value,h5in,data_step,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  integer(i4), intent(in) :: value
  integer, intent(in)   :: data_step ! length of array to preserve
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine add_h5_int
!-----------------------------------------------------------------------
! subprogram 30. add_h5_dbl
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a 1D array
!-----------------------------------------------------------------------
  subroutine add_h5_dbl(inid,aname,value,h5in,data_step,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, intent(in) :: value
  integer, intent(in)   :: data_step ! length of array to preserve
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine add_h5_dbl
!-----------------------------------------------------------------------
! subprogram 30. add_h5_int_1d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_int_1d(inid,aname,array,h5in,data_step,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  integer(i4), dimension(:), intent(in) :: array
  integer, intent(in)   :: data_step ! length of array to preserve
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine add_h5_int_1d
!-----------------------------------------------------------------------
! subprogram 30. add_h5_1d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_1d(inid,aname,array,h5in,data_step,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:), intent(in) :: array
  integer, intent(in)   :: data_step ! length of array to preserve
  type(hdf5ErrorType), intent(inout) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine add_h5_1d
!-----------------------------------------------------------------------
! subprogram 30. add_h5_2d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_2d(inid,aname,array,h5in,data_step,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:,:), intent(in) :: array
  integer, intent(in)   :: data_step ! length of array to preserve
  type(hdf5ErrorType), intent(inout) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine add_h5_2d

!-----------------------------------------------------------------------
! subprogram 30a. add_h5_3d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_3d(inid,aname,array,h5in,data_step,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:,:,:), intent(in) :: array
  integer, intent(in)   :: data_step ! length of array to preserve
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine add_h5_3d
!-----------------------------------------------------------------------
! subprogram 30a. add_h5_4d
! This adds data to an unlimited data space.  This is used for
! writing things like time data
! You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
  subroutine add_h5_4d(inid,aname,array,h5in,data_step,errval)
  integer(HID_T), intent(in) :: inid
  character(*), intent(in) :: aname
  double precision, dimension(:,:,:,:), intent(in) :: array
  integer, intent(in)   :: data_step ! length of array to preserve
  TYPE(hdf5ErrorType) :: errval
  TYPE(hdf5InOpts), intent(in) :: h5in
  errval%errBool = .false.
  return
  end subroutine add_h5_4d
!-----------------------------------------------------------------------
! subprogram 36. read_dims
! Read the dimensions of dataset associated with aname: 1d array
!-----------------------------------------------------------------------
  subroutine read_dims(dset_id,dims,errval)
  integer(HID_T), intent(in) :: dset_id
  integer(HSIZE_T), dimension(:), intent(inoUT) :: dims
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_dims

!-----------------------------------------------------------------------
! subprogram 41. read_h5_intl
! Read simple data set: 1d array
!-----------------------------------------------------------------------
  subroutine read_h5_intl(fid,aname,value,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer(i8), intent(inoUT) :: value
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_h5_intl

!-----------------------------------------------------------------------
! subprogram 41. read_h5_intl_1d
! Read simple data set: 1d array
!-----------------------------------------------------------------------
  subroutine read_h5_intl_1d(fid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer(i8), dimension(:), intent(inoUT) :: array
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_h5_intl_1d

!-----------------------------------------------------------------------
! subprogram 37. read_h5_1d
! Read simple data set: 1d array
!-----------------------------------------------------------------------
  subroutine read_h5_1d(fid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  TYPE(hdf5InOpts), intent(in) :: h5in
  real(r8), dimension(:), intent(inoUT) :: array
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_h5_1d

!-----------------------------------------------------------------------
! subprogram 42. read_h5_2d
! Read simple data set: 2d array
!-----------------------------------------------------------------------
  subroutine read_h5_2d(fid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  real(r8), dimension(:,:), intent(inoUT) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_h5_2d

!-----------------------------------------------------------------------
! subprogram 39. read_h5_3d
! Read simple data set: 3d array
!-----------------------------------------------------------------------
  subroutine read_h5_3d(fid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  real(r8), dimension(:,:,:), intent(inoUT) :: array
  TYPE(hdf5InOpts), intent(in) :: h5in
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_h5_3d

!-----------------------------------------------------------------------
! subprogram 41. read_attribute_intl_sc
! Read real scalar attribute
!-----------------------------------------------------------------------
  subroutine read_attribute_intl_sc(fid,aname,val,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer(i8), intent(inoUT) :: val
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_attribute_intl_sc

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
  errval%errBool = .false.
  return
  end subroutine read_attribute_rl_sc

!-----------------------------------------------------------------------
! subprogram 41. read_attribute_intl_vec
! Read integer vector attribute
!-----------------------------------------------------------------------
  subroutine read_attribute_intl_vec(fid,aname,array,h5in,errval)
  integer(HID_T), intent(in) :: fid
  character*(*), intent(in) :: aname
  TYPE(hdf5InOpts), intent(in) :: h5in
  integer(i8), dimension(:), intent(inoUT) :: array
  TYPE(hdf5ErrorType), intent(inoUT) :: errval
  errval%errBool = .false.
  return
  end subroutine read_attribute_intl_vec
!-----------------------------------------------------------------------
! close module.
!-----------------------------------------------------------------------
  end module hdf5_api
