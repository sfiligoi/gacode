!------------------------------------------------------
! gyro_write_timedata_wedge_hdf5.f90
!
! PURPOSE:
!  HDF5 output of wedge data.
!-----------------------------------------------------

subroutine gyro_write_timedata_wedge_hdf5

  use gyro_globals
  use hdf5_api

  !---------------------------------------------------
  implicit none
  include 'mpif.h'
  !
  integer, parameter :: hr4=SELECTED_REAL_KIND(6,37)
  !
  real :: pi=3.141592653589793
  !
  character(60) :: description
  character(64) :: step_name
  character(128) :: dumpfile
  character(128) :: meshfile
  integer(HID_T) :: fidwedge,gidwedge
  integer :: n_wedge
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label


  !---------------------------------------------------
  ! Grid
  !
  n_wedge = n_theta_plot*n_theta_mult

  if (i_proc == 0) then
     if (n_wedge <= 1) then
        write(*,*) "Wedge caluculations need n_theta_plot*n_theta_mult > 1."
        return 
     endif
  endif

  !---------------------------------------------------
  ! SEK: Should I do this every time?
  !---------------------------------------------------
  if (i_proc == 0) then
     call vshdf5_inith5vars(h5in, h5err)
     h5in%comm=MPI_COMM_SELF
     h5in%info=MPI_INFO_NULL
     h5in%wrd_type=H5T_NATIVE_REAL
     !h5in%wrd_type=h5in%h5_kind_type_r4
     h5in%typeConvert=.true.
     !h5in%wrd_type=H5T_NATIVE_DOUBLE
     h5in%doTranspose=.true.
     h5in%verbose=.true.
     h5in%debug=.false.
     h5in%wrVsTime=.true.
     h5in%vsTime=t_current
     h5in%vsStep=step

     !---------------------------------------------------
     ! Timestep data:
     !
     number_label=NINT(t_current/dt)
     if (number_label>999999) THEN
        write(step_name,fmt='(i7.7)') number_label
     else if (number_label>99999) THEN
        write(step_name,fmt='(i6.6)') number_label
     else
        write(step_name,fmt='(i5.5)') number_label
     endif

     dumpfile=TRIM(path)//"gyrowedge"//TRIM(step_name)//".h5"
     description="GYRO real space file at single phi plane.  On wedge mesh"
     call open_newh5file(dumpfile,fidwedge,description,gidwedge,h5in,h5err)

     ! make link to grid file

     call make_external_link(TRIM(meshfile),"wedgeMesh", &
        gidwedge,"wedgeMesh", h5in,h5err)


  endif
  !---------------------------------------------------

  !--------------------------------------------------
  ! Output of field-like quantities:
  !
  !

  if (plot_n_flag == 1) then

     ! DENSITY
     h5in%units=" "
     h5in%mesh="/cartMesh"
     call write_distributed_complex_sorf_h5("density",&
          gidwedge,gidwedge,&
          n_wedge*n_x*n_kinetic,&
          n_wedge,n_x,n_kinetic,&
          moments_plot_wedge(:,:,:,1),&
          .true.,&
          .true.,&
          h5in,h5err)
  endif

  if (plot_e_flag == 1) then
     ! ENERGY
     h5in%units="energy units"
     h5in%mesh="/cartMesh"
     call write_distributed_complex_sorf_h5("energy",&
          gidwedge,gidwedge,&
          n_wedge*n_x*n_kinetic,&
          n_wedge,n_x,n_kinetic,&
          moments_plot_wedge(:,:,:,2),&
          .true.,&
          .true.,&
          h5in,h5err)
  endif

  if (plot_v_flag == 1) then
     ! PARALLEL VELOCITY
     h5in%units="vpar units"
     h5in%mesh="/cartMesh"
     call write_distributed_complex_sorf_h5("v_par",&
          gidwedge,gidwedge,&
          n_wedge*n_x*n_kinetic,&
          n_wedge,n_x,n_kinetic,&
          moments_plot_wedge(:,:,:,3),&
          .true.,&
          .true.,&
          h5in,h5err)
  endif

  !--------------------------------------------------

  if (i_proc == 0) then
     call close_h5file(fidwedge,gidwedge,h5err)
  endif

  return

end subroutine gyro_write_timedata_wedge_hdf5

