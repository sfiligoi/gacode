!------------------------------------------------------
! gyro_write_timedata_wedge_hdf5.f90
!
! PURPOSE:
!  HDF5 output of wedge data.
!-----------------------------------------------------

subroutine gyro_write_timedata_wedge_hdf5

  use gyro_globals
  !use hdf5
  use hdf5_api
  use gyro_vshdf5_mod

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

     call hdf5_write_wedge_coords
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
     call write_distributed_complex_h5("density",&
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
     call write_distributed_complex_h5("energy",&
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
     call write_distributed_complex_h5("v_par",&
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

contains
  subroutine hdf5_write_wedge_coords
    use GEO_interface
    !------------------------------------------
    !  Write the coordinates out
    !  We want to have same coordinate system as:
    !    allocate(phi_plot(n_theta_plot,n_x,n_field+eparallel_plot_flag))
    !  This should be generalized to include the other GEO options
    !------------------------------------------
    real, dimension(:,:), allocatable :: Rf,Zf
    real, dimension(:,:,:), allocatable :: bufferwedgeMesh
    real :: rmajc,zmagc,kappac,deltac,zetac,xdc,rhoc
    real :: dRdr,dZdr,dkappadr,ddeltadr,dzetadr 
    real, dimension (:), allocatable :: theta,r_c
    real, dimension (:,:), allocatable :: dRdrho,DRdtheta
    real, dimension (:,:), allocatable :: dZdrho,DZdtheta
    real, dimension (:,:), allocatable :: rtJacobian
    real :: zeta_wedge
    integer :: iphi,ix,j,ncoarse,nwedge

    ncoarse = n_theta_plot
    nwedge = n_theta_plot*n_theta_mult
    allocate(Rf(1:nwedge,n_x), Zf(1:nwedge,n_x))
    allocate(theta(1:nwedge), r_c(1:n_x))

    allocate(dRdrho(1:nwedge,1:n_x))
    allocate(dZdrho(1:nwedge,1:n_x))
    allocate(dRdtheta(1:nwedge,1:n_x))
    allocate(dZdtheta(1:nwedge,1:n_x))
    allocate(rtJacobian(1:nwedge,1:n_x))

    !----------------------------------------
    ! Calculate the R,Z coordinates.  See write_geometry_arrays.f90
    !---------------------------------------- 

    do ix=1,n_x
       if (flat_profile_flag == 0) then
          r_c(ix)=r_s(ix)
       else
          r_c(ix)=r(ix)
       endif
    enddo
    do j=1,nwedge
       ! This needs to match up with what's in gyro_set_blend_arrays.f90
       theta(j)=theta_wedge_offset+real(j-1)*theta_wedge_angle/       &
            real(n_theta_plot*n_theta_mult-1)
       !theta = -pi+REAL(j)*pi*2./REAL(nwedge)
    enddo

    do ix=1,n_x
       rhoc     = r_c(ix)
       rmajc  = rmaj_s(ix)
       dRdr   = drmaj_s(ix)
       zmagc  = zmag_s(ix)
       dZdr   = dzmag_s(ix)
       kappac = kappa_s(ix)
       dkappadr = (kappac/rhoc)*s_kappa_s(ix)
       deltac   = delta_s(ix)
       ddeltadr = (1./rhoc)*s_delta_s(ix)
       xdc    = asin(deltac)
       zetac  = zeta_s(ix)
       dzetadr = (1./rhoc)*s_zeta_s(ix)

       do j=1,nwedge
          if(radial_profile_method==1) then
             Rf(j,ix)=rmajc+rhoc*cos(theta(j))
             Zf(j,ix)=zmagc+rhoc*sin(theta(j))
          else
             Rf(j,ix)=rmajc+rhoc*cos(theta(j)+xdc*sin(theta(j)))
             Zf(j,ix)=zmagc+kappac*rhoc*sin(theta(j)+zetac*sin(2.*theta(j)))
          endif
          !pieces of jacobian 
          dRdrho(j,ix) = COS(theta(j)+xdc*SIN(theta(j))) + &
               dRdr - (1./SQRT(1-deltac**2))* &
               rhoc*SIN(theta(j))*SIN(theta(j)+xdc*SIN(theta(j)))*ddeltadr

          dZdrho(j,ix) = SIN(theta(j)+SIN(2.*theta(j))*zetac*kappac) + &
               dZdr + rhoc*COS(theta(j)+SIN(2.*theta(j))*zetac)*SIN(2.*theta(j))* &
               kappac*dzetadr + &
               rhoc*SIN(theta(j)+SIN(2.*theta(j))*zetac)*dkappadr

          dRdtheta(j,ix) = -1.*rhoc*(1.+xdc*COS(theta(j)))* &
               SIN(theta(j)+xdc)*SIN(theta(j))

          dZdtheta(j,ix)= rhoc*COS(theta(j)+SIN(2.*theta(j))*zetac) * &
               (1.+2.*COS(2*theta(j))*zetac*kappac)

          rtJacobian(j,ix) = dRdrho(j,ix)*dZdtheta(j,ix) - &
               dZdrho(j,ix)*dRdtheta(j,ix)

       enddo
    enddo

    !-------------------------------------------------------
    ! Set up the alpha grid
    ! These are saved in a module so no need to recalculate
    !-------------------------------------------------------

    if (.not. allocated(alpha_phi_wedge) ) then
       allocate(alpha_phi_wedge(nwedge,n_x,n_torangle_wedge))
       do iphi=1,n_torangle_wedge
          !Don't store zeta_wedge b/c analysis is on a plane by plane basis
          zeta_wedge=REAL(iphi-1)/REAL(n_torangle_wedge)*2.*pi
          alpha_phi_wedge(:,:,iphi)=torangle_offset+zeta_wedge+nu_wedge(:,:)
       end do
    endif

    !----------------------------------------
    ! Dump the wedge meshes
    !---------------------------------------- 

    h5in%units=""
    call dump_h5(gidwedge,'Rgyro',Rf,h5in,h5err)
    call dump_h5(gidwedge,'Zgyro',Zf,h5in,h5err)
    call dump_h5(gidwedge,'torangle_offset',torangle_offset,h5in,h5err)
    call dump_h5(gidwedge,'alpha',alpha_phi_wedge,h5in,h5err)
    h5in%units="m"
    call dump_h5(gidwedge,'R',Rf*a_meters,h5in,h5err)
    call dump_h5(gidwedge,'Z',Zf*a_meters,h5in,h5err)
    call dump_h5(gidwedge,'r_min',r_c*a_meters,h5in,h5err)
    call dump_h5(gidwedge,'rtJacobian',rtJacobian*a_meters,h5in,h5err)
    h5in%units="radians"
    call dump_h5(gidwedge,'theta',theta,h5in,h5err)

    ! For ease of use, have a single data set that has R,Z. 
    allocate(bufferwedgeMesh(nwedge,n_x,2))
    bufferwedgeMesh(:,:,1) = Rf(:,:)*a_meters
    bufferwedgeMesh(:,:,2) = Zf(:,:)*a_meters
    h5in%units="m"
    h5in%mesh="mesh-structured"
    call dump_h5(gidwedge,'cartMesh',bufferwedgeMesh(:,:,:),h5in,h5err)
    h5in%mesh=""
    deallocate(bufferwedgeMesh)

    !----------------------------------------
    ! 
    !---------------------------------------- 
    deallocate(Rf, Zf)

    deallocate(theta, dRdrho,dZdrho,dRdtheta,dZdtheta)

    return
  end subroutine hdf5_write_wedge_coords

end subroutine gyro_write_timedata_wedge_hdf5

