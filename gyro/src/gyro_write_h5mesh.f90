  !------------------------------------------------------
  subroutine hdf5_write_coords
  use mpi
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
    integer :: iphi,ncoarse,nphi,nwedge

    character(60) :: description
    character(128) :: dumpfile

    integer(HID_T) :: gridFileID,gridGroupID
    integer(HID_T) :: poloidalGridID, threeDGridID, wedgeGridID 
    type(hdf5InOpts) :: h5in
    type(hdf5ErrorType) :: h5err
    integer :: number_label
    logical :: write_threed
    logical :: h5_rewind=.false.
    
    real :: pi=3.141592653589793

  !wedge w/ duplicates **CLEAN UP**
    real, dimension(:,:), allocatable :: Rf,Zf
    real, dimension(:,:,:), allocatable :: bufferwedgeMesh
    real :: rhoc
    real :: dRdr,dZdr,dkappadr,ddeltadr,dzetadr 
    real, dimension (:), allocatable :: w_theta,w_r_c,zeta_phi
    real, dimension (:,:), allocatable :: dRdrho,DRdtheta
    real, dimension (:,:), allocatable :: dZdrho,DZdtheta
    real, dimension (:,:), allocatable :: rtJacobian
    real :: zeta_wedge

!---------------------------------------------------
      ! Determine if the 3D files need to be written 
      if (n_torangle_3d > 1 ) then
         write_threed = .true.
      else
         write_threed = .false.
      endif

    !---------------------------------------
    !0) init h5in and set things 

         call vshdf5_inith5vars(h5in, h5err)
         h5in%comm=MPI_COMM_SELF
         h5in%info=MPI_INFO_NULL
         h5in%wrd_type=H5T_NATIVE_REAL
         h5in%typeConvert=.true.
         h5in%doTranspose=.true.
         h5in%verbose=.true.
         h5in%debug=.false.
         h5in%wrVsTime=.true.
         h5in%vsTime=t_current
         h5in%vsStep=step


         if (debug_flag==1) h5in%debug=.true.
 
    !---------------------------------------
    ! 1) open grid file 
         dumpfile=TRIM(path)//"gyroMesh.h5" 
         description="Store GYRO physical grids."
         call open_newh5file(dumpfile,gridFileID,description, &
            gridGroupID,h5in,h5err)
         if (h5err%errBool) call catch_error(h5err%errorMsg)


      !1.1) make groups
         call make_mesh_group(gridGroupID,"poloidalMesh",poloidalGridID,h5in,&
              h5err)
         call make_mesh_group(gridGroupID,"threeDMesh",threeDGridID,h5in,&
              h5err)
         call make_mesh_group(gridGroupID,"wedgeMesh",wedgeGridID,h5in,&
              h5err)


 

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
    call dump_h5(poloidalGridID,'Rgyro',Rc,h5in,h5err)
    call dump_h5(poloidalGridID,'Zgyro',Zc,h5in,h5err)
    h5in%units="m"
    call dump_h5(poloidalGridID,'R',Rc*a_meters,h5in,h5err)
    call dump_h5(poloidalGridID,'Z',Zc*a_meters,h5in,h5err)

    ! Here we do not repeat the points since this is the grid
    ! that will be used for the mode plots on thete E [0,2 pi)
    allocate(bufferMesh(0:ncoarse,n_x,2))
    bufferMesh(:,:,1)= Rc*a_meters
    bufferMesh(:,:,2)= Zc*a_meters
    h5in%units="m"
    h5in%mesh="mesh-structured"
    call dump_h5(poloidalGridID,'cartMesh',bufferMesh(:,:,:),h5in,h5err)
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
     call dump_h5(threeDGridID,'R',Rc*a_meters,h5in,h5err)
     call dump_h5(threeDGridID,'Z',Zc*a_meters,h5in,h5err)
     h5in%units="radians"
     call dump_h5(threeDGridID,'torAngle',zeta_phi,h5in,h5err)
     call dump_h5(threeDGridID,'torangle_offset',torangle_offset,h5in,h5err)
     call dump_h5(threeDGridID,'alpha',alpha_phi,h5in,h5err)

     allocate(buffer(0:ncoarse,n_x,n_torangle_3d,3))
      buffer = -9999.999
     do iphi=1,n_torangle_3d
        buffer(:,:,iphi,1)= Rc(:,:)*COS(zeta_phi(iphi))
        buffer(:,:,iphi,2)=-Rc(:,:)*SIN(zeta_phi(iphi))
        buffer(:,:,iphi,3)= Zc(:,:)
     enddo

     h5in%units="m"; h5in%mesh="mesh-structured"
     call dump_h5(threeDGridID,'cartMesh',buffer*a_meters,h5in,h5err)
     deallocate(buffer)
    endif

    if (allocated(Rc)) deallocate(Rc)
    if (allocated(Zc)) deallocate(Zc)
    if (allocated(zeta_phi)) deallocate(zeta_phi)


!----------------------dump wedge grid ------------------

  wedge: if (time_skip_wedge > 0) then
    
      ncoarse = n_theta_plot
      nwedge = n_theta_plot*n_theta_mult
      allocate(Rf(1:nwedge,n_x), Zf(1:nwedge,n_x))
      allocate(w_theta(1:nwedge), w_r_c(1:n_x))

      allocate(dRdrho(1:nwedge,1:n_x))
      allocate(dZdrho(1:nwedge,1:n_x))
      allocate(dRdtheta(1:nwedge,1:n_x))
      allocate(dZdtheta(1:nwedge,1:n_x))
      allocate(rtJacobian(1:nwedge,1:n_x))



          !----------------------------------------
      ! Calculate the R,Z coordinates.  
      !---------------------------------------- 

      do ix=1,n_x
         if (flat_profile_flag == 0) then
            w_r_c(ix)=r_s(ix)
         else
            w_r_c(ix)=r(ix)
         endif
      enddo
      do j=1,nwedge
         ! This needs to match up with what's in gyro_set_blend_arrays.f90
         w_theta(j)=theta_wedge_offset+real(j-1)*theta_wedge_angle/       &
              real(n_theta_plot*n_theta_mult-1)
         !theta = -pi+REAL(j)*pi*2./REAL(nwedge)
      enddo

      do ix=1,n_x
         rhoc     = w_r_c(ix)
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
               Rf(j,ix)=rmajc+rhoc*cos(w_theta(j))
               Zf(j,ix)=zmagc+rhoc*sin(w_theta(j))
            else
               Rf(j,ix)=rmajc+rhoc*cos(w_theta(j)+xdc*sin(w_theta(j)))
               Zf(j,ix)=zmagc+kappac*rhoc*sin(w_theta(j)+zetac*sin(2.*w_theta(j)))
            endif
            !pieces of jacobian 
            dRdrho(j,ix) = COS(w_theta(j)+xdc*SIN(w_theta(j))) + &
                 dRdr - (1./SQRT(1-deltac**2))* &
                 rhoc*SIN(w_theta(j))*SIN(w_theta(j)+xdc*SIN(w_theta(j)))*ddeltadr

            dZdrho(j,ix) = SIN(w_theta(j)+SIN(2.*w_theta(j))*zetac*kappac) + &
                 dZdr + rhoc*COS(w_theta(j)+SIN(2.*w_theta(j))*zetac)*SIN(2.*w_theta(j))* &
                 kappac*dzetadr + &
                 rhoc*SIN(w_theta(j)+SIN(2.*w_theta(j))*zetac)*dkappadr

            dRdtheta(j,ix) = -1.*rhoc*(1.+xdc*COS(w_theta(j)))* &
                 SIN(w_theta(j)+xdc)*SIN(w_theta(j))

            dZdtheta(j,ix)= rhoc*COS(w_theta(j)+SIN(2.*w_theta(j))*zetac) * &
                 (1.+2.*COS(2*w_theta(j))*zetac*kappac)

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
      call dump_h5(wedgeGridID,'Rgyro',Rf,h5in,h5err)
      call dump_h5(wedgeGridID,'Zgyro',Zf,h5in,h5err)
      call dump_h5(wedgeGridID,'torangle_offset',torangle_offset,h5in,h5err)
      call dump_h5(wedgeGridID,'alpha',alpha_phi_wedge,h5in,h5err)
      h5in%units="m"
      call dump_h5(wedgeGridID,'R',Rf*a_meters,h5in,h5err)
      call dump_h5(wedgeGridID,'Z',Zf*a_meters,h5in,h5err)
      call dump_h5(wedgeGridID,'r_min',w_r_c*a_meters,h5in,h5err)
      call dump_h5(wedgeGridID,'rtJacobian',rtJacobian*a_meters,h5in,h5err)
      h5in%units="radians"
      call dump_h5(wedgeGridID,'w_theta',w_theta,h5in,h5err)

      ! For ease of use, have a single data set that has R,Z. 
      allocate(bufferwedgeMesh(nwedge,n_x,2))
      bufferwedgeMesh(:,:,1) = Rf(:,:)*a_meters
      bufferwedgeMesh(:,:,2) = Zf(:,:)*a_meters
      h5in%units="m"
      h5in%mesh="mesh-structured"
      call dump_h5(wedgeGridID,'cartMesh',bufferwedgeMesh(:,:,:),h5in,h5err)
      h5in%mesh=""
      deallocate(bufferwedgeMesh)

      !----------------------------------------
      ! 
      !---------------------------------------- 
      deallocate(Rf, Zf)
      deallocate(w_theta, dRdrho,dZdrho,dRdtheta,dZdtheta)

    endif wedge

         call close_group("poloidalMesh",poloidalGridID,h5err)
         call close_group("threeDMesh",threeDGridID,h5err)
         call close_group("wedgeMesh",wedgeGridID,h5err)


         call close_h5file(gridFileID,gridGroupID,h5err)






  end subroutine hdf5_write_coords

