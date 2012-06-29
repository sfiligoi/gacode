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
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
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

