  !------------------------------------------------------
  subroutine hdf5_write_coords
  use gyro_globals
  use hdf5_api
  use GEO_interface
    !------------------------------------------
    !  Write the coordinates out
    !  We want to have same coordinate system as:
    !    allocate(phi_plot(n_theta_plot,n_x,n_field+eparallel_plot_flag))
    !  This should be generalized to include the other GEO options
    !------------------------------------------
    real, dimension(:,:), allocatable :: Rc,Zc
    real, dimension(:,:,:,:), allocatable :: buffer
    real, dimension(:,:,:), allocatable :: bufferMesh
    real :: theta,rmajc,zmagc,kappac,deltac,zetac,r_c,xdc
    integer :: iphi,ix,j,ncoarse,nphi

  integer(HID_T) :: dumpGid,dumpFid,gid3D,fid3D
  integer(HID_T) :: dumpTGid,dumpTFid
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label
  logical :: write_threed
  logical :: h5_rewind=.false.


    ncoarse = n_theta_plot
    allocate(Rc(0:ncoarse,n_x), Zc(0:ncoarse,n_x))

    !----------------------------------------
    ! Calculate the R,Z coordinates.  See write_geometry_arrays.f90
    ! The theta grid needs to correspond to the the theta_plot
    ! array which sets the interpolation arrays in 
    ! gyro_set_blend_arrays.  The theta_plot array is defined as:
    !  do j=1,n_theta_plot
    !          theta_plot(j) = -pi+(j-1)*pi_2/n_theta_plot
    !  enddo
    ! such that theta E [0,2 pi) in gyro_banana_operators.f90
    ! For the 3D arrays, we want the periodic point repeated for
    ! nice plots; i.e., theta E [0,2 pi], but we plot the raw
    ! mode data on theta E [0, 2 pi).  Can be a bit confusing.
    !---------------------------------------- 

    do ix=1,n_x
       r_c=r(ix)
       rmajc = rmaj_s(ix)
       zmagc = zmag_s(ix)
       kappac = kappa_s(ix)
       deltac = delta_s(ix)
       xdc    = asin(deltac)
       zetac  = zeta_s(ix)
       ! Note:  This needs to be generalized for all geometries
       do j=0,ncoarse
          theta = -pi+REAL(j)*pi*2./REAL(ncoarse)
          if(radial_profile_method==1) then
             Rc(j,ix)=rmajc+r_c*cos(theta)
             Zc(j,ix)=zmagc+r_c*sin(theta)
          else
             Rc(j,ix)=rmajc+r_c*cos(theta+xdc*sin(theta))
             Zc(j,ix)=zmagc+kappac*r_c*sin(theta+zetac*sin(2.*theta))
          endif
       enddo
    enddo

    !----------------------------------------
    ! Dump the coarse meshes
    !---------------------------------------- 

    h5in%units=""
    call dump_h5(dumpGid,'Rgyro',Rc,h5in,h5err)
    call dump_h5(dumpGid,'Zgyro',Zc,h5in,h5err)
    h5in%units="m"
    call dump_h5(dumpGid,'R',Rc*a_meters,h5in,h5err)
    call dump_h5(dumpGid,'Z',Zc*a_meters,h5in,h5err)

    ! Here we do not repeat the points since this is the grid
    ! that will be used for the mode plots on thete E [0,2 pi)
    allocate(bufferMesh(0:ncoarse,n_x,2))
    bufferMesh(:,:,1)= Rc*a_meters
    bufferMesh(:,:,2)= Zc*a_meters
    h5in%units="m"
    h5in%mesh="mesh-structured"
    call dump_h5(dumpGid,'cartMesh',bufferMesh(:,:,:),h5in,h5err)
    h5in%mesh=" "
    deallocate(bufferMesh)

    !----------------------------------------
    ! Dump the coarse mesh(es) in 3D
    !---------------------------------------- 
    if (write_threed) then
    !------------------------------------------------
    ! Set up the toroidal grid.  Only used for coarse grid
    !-------------------------------------------------
      allocate(zeta_phi(n_torangle_3d))
      do iphi=1,n_torangle_3d
         zeta_phi(iphi)=REAL(iphi-1)/REAL(n_torangle_3d-1)*2.*pi
      end do

     !-------------------------------------------------------
     ! Set up the alpha grid
     ! These are set up in a module so no need to recalculate
     !-------------------------------------------------------
     if (.not. allocated(alpha_phi) ) then 
        nphi=1
        if (n_torangle_3d > 0 ) nphi=n_torangle_3d 
        allocate(alpha_phi(0:ncoarse,n_x,nphi))
        do iphi=1,n_torangle_3d
           alpha_phi(:,:,iphi)=zeta_phi(iphi)+nu_coarse(:,:)
        end do
     endif

     h5in%units="m"
     call dump_h5(gid3d,'R',Rc*a_meters,h5in,h5err)
     call dump_h5(gid3d,'Z',Zc*a_meters,h5in,h5err)
     h5in%units="radians"
     call dump_h5(gid3d,'torAngle',zeta_phi,h5in,h5err)
     call dump_h5(gid3d,'torangle_offset',torangle_offset,h5in,h5err)
     call dump_h5(gid3d,'alpha',alpha_phi,h5in,h5err)

     allocate(buffer(0:ncoarse,n_x,n_torangle_3d,3))
      buffer = -9999.999
     do iphi=1,n_torangle_3d
        buffer(:,:,iphi,1)= Rc(:,:)*COS(zeta_phi(iphi))
        buffer(:,:,iphi,2)=-Rc(:,:)*SIN(zeta_phi(iphi))
        buffer(:,:,iphi,3)= Zc(:,:)
     enddo
!    do iphi=1,n_torangle_3d
!      do j=0,ncoarse
!        do ix=1,n_x
!          buffer(j,ix,iphi,1)= Rc(j,ix)*COS(zeta_phi(iphi))
!!	 write(*,*) "j=",j,"ix=",ix,"iphi=",iphi &
!!		,"Rc=",Rc(j,ix),"zeta=",zeta_phi(iphi)
!!	write(*,*) " buffer(j,ix,iphi,1) = ", buffer(j,ix,iphi,1)
!	buffer(j,ix,iphi,2)=-Rc(j,ix)*SIN(zeta_phi(iphi))
!          buffer(j,ix,iphi,3)= Zc(j,ix)
!        enddo
!      enddo
!     enddo


     h5in%units="m"; h5in%mesh="mesh-structured"
 !    call dump_h5_4d(gid3d,'cartMesh',buffer*a_meters,h5in,h5err)
! This is what I want for the external linking.
!     call dump_h5(gid3d,'cartMesh',figGrid,gidGrid,'threeDMesh',buffer*a_meters,h5in,h5err)
     call dump_h5(gid3d,'cartMesh',buffer*a_meters,h5in,h5err)
     deallocate(buffer)
    endif

    !----------------------------------------
    ! Dump the wedge mesh(es)
    !---------------------------------------- 
    if (allocated(Rc)) deallocate(Rc)
    if (allocated(Zc)) deallocate(Zc)
    if (allocated(zeta_phi)) deallocate(zeta_phi)

  end subroutine hdf5_write_coords
