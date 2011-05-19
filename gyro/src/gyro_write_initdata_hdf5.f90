 !------------------------------------------
 !  Data that does not change with time.  
 !  It is equivalent to:
 !    profile_vugyro.out
 !------------------------------------------

subroutine gyro_write_initdata_hdf5(datafile,action)

  use gyro_globals
  use hdf5_api
  use math_constants
  use GEO_interface
  !------------------------------------------
  implicit none
  include 'mpif.h'
  !
  integer, intent(in) :: action
  character (len=*), intent(in) :: datafile
  !
  integer :: data_loop
  integer :: i_dummy
  real :: dummy
  integer(HID_T) :: fid, rootid
  character(90) :: description
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: n_wedge
  integer :: io_mode
  real :: theta
  real :: dr, buff
  real :: kt
  real, allocatable :: buffer(:,:,:)
  !double precision :: buff
  !------------------------------------------
  ! Do the initialization here.  Might need
  ! better logic here based on action.
  !------------------------------------------
  call vshdf5_fcinit()
  if(i_proc/=0) return
  call vshdf5_inith5vars(h5in, h5err)
  h5in%comm=MPI_COMM_SELF
  h5in%info=MPI_INFO_NULL
  h5in%typeConvert=.true.
  h5in%wrd_type=H5T_NATIVE_REAL
  !h5in%wrd_type=H5T_NATIVE_DOUBLE
  h5in%doTranspose=.true.
  !h5in%vsTime=intime
  h5in%wrVsTime=.false.
  h5in%verbose=.true.

  !---------------------------------------------------------------------
  ! Write the variables to an hdf5 file
  ! These variables are essentially the write_profile_vugyro.f90 
  !---------------------------------------------------------------------
  description=" "
  call open_newh5file(datafile,fid,description,rootid,h5in,h5err)

  h5in%mesh=" "; h5in%units=" "
  call dump_h5(rootid,"n_x",n_x,h5in,h5err)
  call dump_h5(rootid,"n_theta_section",n_theta_section,h5in,h5err)
  call dump_h5(rootid,"n_pass",n_pass,h5in,h5err)
  call dump_h5(rootid,"n_trap",n_trap,h5in,h5err)
  call dump_h5(rootid,"n_energy",n_energy,h5in,h5err)
  call dump_h5(rootid,"n_theta_plot",n_theta_plot,h5in,h5err)
  call dump_h5(rootid,"n0",n0,h5in,h5err)
  call dump_h5(rootid,"n_n",n_n,h5in,h5err)
  call dump_h5(rootid,"d_n",d_n,h5in,h5err)
  call dump_h5(rootid,"nonlinear_flag",nonlinear_flag,h5in,h5err)
  call dump_h5(rootid,"electron_method",electron_method,h5in,h5err)
  call dump_h5(rootid,"n_field",n_field,h5in,h5err)
  call dump_h5(rootid,"n_ion",n_ion,h5in,h5err)
  call dump_h5(rootid,"n_kinetic",n_kinetic,h5in,h5err)
  call dump_h5(rootid,"n_spec",n_spec,h5in,h5err)
  call dump_h5(rootid,"field_r0_flag",field_r0_flag,h5in,h5err)
  call dump_h5(rootid,"field_r0_grid",field_r0_grid,h5in,h5err)
  call dump_h5(rootid,"boundary_method",boundary_method,h5in,h5err)
  call dump_h5(rootid,"r",r,h5in,h5err)
  call dump_h5(rootid,"q",q,h5in,h5err)
  call dump_h5(rootid,"Rmaj_s",rmaj_s,h5in,h5err)
  call dump_h5(rootid,"r_s",r_s,h5in,h5err)
  call dump_h5(rootid,"q_s",q_s,h5in,h5err)
  call dump_h5(rootid,"dlntdr_s",dlntdr_s,h5in,h5err)
  call dump_h5(rootid,"dlnndr_s",dlnndr_s,h5in,h5err)
  call dump_h5(rootid,"tem_s",tem_s,h5in,h5err)
  call dump_h5(rootid,"den_s",den_s,h5in,h5err)
  call dump_h5(rootid,"aspect_s",rmaj_s/r_s,h5in,h5err)
  call dump_h5(rootid,"delta_s",delta_s,h5in,h5err)
  call dump_h5(rootid,"zeta_s",zeta_s,h5in,h5err)
  call dump_h5(rootid,"kappa_s",kappa_s,h5in,h5err)
  call dump_h5(rootid,"drmaj_s",drmaj_s,h5in,h5err)
  call dump_h5(rootid,"shat_s",shat_s,h5in,h5err)
  call dump_h5(rootid,"s_delta_s",s_delta_s,h5in,h5err)
  call dump_h5(rootid,"s_zeta_s",s_zeta_s,h5in,h5err)
  call dump_h5(rootid,"s_kappa_s",s_kappa_s,h5in,h5err)
  call dump_h5(rootid,"zmag_s",zmag_s,h5in,h5err)
  call dump_h5(rootid,"dzmag_s",dzmag_s,h5in,h5err)
  call dump_h5(rootid,"beta_unit_s",beta_unit_s,h5in,h5err)
  call dump_h5(rootid,"gamma_e_s",gamma_e_s,h5in,h5err)
  call dump_h5(rootid,"gamma_p_s",gamma_p_s,h5in,h5err)
  call dump_h5(rootid,"mach_s",mach_s,h5in,h5err)
  call dump_h5(rootid,"b_unit_s",b_unit_s,h5in,h5err)
  call dump_h5(rootid,"dr_eodr",dr_eodr,h5in,h5err)
  call dump_h5(rootid,"z_eff_s",z_eff_s,h5in,h5err)
  call dump_h5(rootid,"nu_s",nu_s,h5in,h5err)
  call dump_h5(rootid,"w0_s",w0_s,h5in,h5err)
  call dump_h5(rootid,"box_multiplier",box_multiplier,h5in,h5err)

  call dump_h5(rootid,"lambda", lambda(ir_norm,:),h5in,h5err)
  call dump_h5(rootid,"energy", energy,h5in,h5err)
  call dump_h5(rootid,"lambda_tp", lambda_tp(ir_norm),h5in,h5err)
  call dump_h5(rootid,"krho_collect", krho_collect(:),h5in,h5err)
  call dump_h5(rootid,"rhos_norm", rhos_norm,h5in,h5err)
  call dump_h5(rootid,"zcharge", z(:),h5in,h5err)
  call dump_h5(rootid,"n_moment", n_moment ,h5in,h5err)

  !---------------------------------------------------------------------
  ! These variables are essentially gyro_write_units.f90
  !---------------------------------------------------------------------
  ! kT in MJ (note the conversion 1.6022e-22 MJ/keV)
  kt = 1.6022e-22*tem_norm

  h5in%units="kg"
  call dump_h5(rootid,'m_ref', 2.*kg_proton ,h5in,h5err)
  h5in%units="Tesla"
  call dump_h5(rootid,'b_unit',b_unit_norm,h5in,h5err)
  h5in%units="m"
  call dump_h5(rootid,'a',a_meters,h5in,h5err)
  h5in%units="1/s"
  call dump_h5(rootid,'csda_norm',csda_norm,h5in,h5err)
  h5in%units="m/s"
  call dump_h5(rootid,'csda_norm_D',csda_norm*a_meters,h5in,h5err)
  h5in%units="keV"
  call dump_h5(rootid,'Te',tem_norm,h5in,h5err)
  h5in%units="10^19/m^3"
  call dump_h5(rootid,'ne',den_norm,h5in,h5err)
  h5in%units="m"
  call dump_h5(rootid,'rho_sD',rhos_norm*a_meters,h5in,h5err)
  h5in%units="m^2/s"
  call dump_h5(rootid,'chi_gBD',csda_norm*(rhos_norm*a_meters)**2,h5in,h5err)
  buff=1.e19*den_norm*(csda_norm*a_meters)*rhos_norm**2/0.624e22
  h5in%units="MW/keV/m^2"
  call dump_h5(rootid,'Gamma_gBD',buff,h5in,h5err)
  h5in%units="MW/m^2"
  buff=1.e19*den_norm*(csda_norm*a_meters)*kt*rhos_norm**2
  call dump_h5(rootid,'Q_gBD',buff,h5in,h5err)
  h5in%units="Nm/m^2"
  buff=1.e19*den_norm*a_meters*kt*rhos_norm**2*1e6
  call dump_h5(rootid,'Pi_gBD',buff,h5in,h5err)
  h5in%units="MW/m^3"
  buff=1.e19*den_norm*csda_norm*kt*rhos_norm**2
  call dump_h5(rootid,'S_gBD',buff,h5in,h5err)

  !---------------------------------------------------------------------
  ! Data from geometry_arrays.out
  !---------------------------------------------------------------------

  n_wedge = n_theta_plot*n_theta_mult

  allocate(buffer(14,1:n_x,n_wedge))
  do i=1,n_x

     if (flat_profile_flag == 0) then

        ! All profiles are global and radial variation is consistent

        GEO_rmin_in      = r_s(i)
        GEO_rmaj_in      = rmaj_s(i)
        GEO_drmaj_in     = drmaj_s(i)
        GEO_zmag_in      = zmag_s(i)
        GEO_dzmag_in     = dzmag_s(i)
        GEO_q_in         = q_s(i)
        GEO_s_in         = shat_s(i)
        GEO_kappa_in     = kappa_s(i)
        GEO_s_kappa_in   = s_kappa_s(i)
        GEO_delta_in     = delta_s(i)
        GEO_s_delta_in   = s_delta_s(i)
        GEO_zeta_in      = zeta_s(i)
        GEO_s_zeta_in    = s_zeta_s(i)
        GEO_beta_star_in = beta_star_s(i)

     else

        ! Profiles are flat and so some parameters need to be linearly extrapolated.

        dr = r(i)-r(ir_norm)

        GEO_rmin_in  = r(i)
        GEO_rmaj_in  = rmaj_s(ir_norm)+drmaj_s(ir_norm)*dr
        GEO_drmaj_in = drmaj_s(ir_norm)
        GEO_zmag_in  = zmag_s(ir_norm)+dzmag_s(ir_norm)*dr
        GEO_dzmag_in = dzmag_s(ir_norm)
        GEO_q_in     = q(i)
        GEO_s_in     = shat_s(ir_norm)
        GEO_kappa_in = kappa_s(ir_norm)+&
             kappa_s(ir_norm)*s_kappa_s(ir_norm)/r(ir_norm)*dr
        GEO_s_kappa_in = s_kappa_s(ir_norm)
        GEO_delta_in   = delta_s(ir_norm)+s_delta_s(ir_norm)/r(ir_norm)*dr
        GEO_s_delta_in = s_delta_s(ir_norm)
        GEO_zeta_in    = zeta_s(ir_norm)+s_zeta_s(ir_norm)/r(ir_norm)*dr
        GEO_s_zeta_in  = s_zeta_s(ir_norm)
        GEO_beta_star_in = beta_star_s(ir_norm)

     endif

     GEO_fourier_in(:,:) = a_fourier_geo_s(:,0:n_fourier_geo,i)
     call GEO_do()

     do j=1,n_wedge

        theta = -pi+(j-1)*pi_2/n_wedge

        ! Test for special case
        if (n_wedge == 1) theta = 0.0

        call GEO_interp(theta)

        buffer(1 ,i,j)=GEO_nu
        buffer(2 ,i,j)=GEO_gsin
        buffer(3 ,i,j)=GEO_gcos1
        buffer(4 ,i,j)=GEO_gcos2
        buffer(5 ,i,j)=GEO_usin
        buffer(6 ,i,j)=GEO_ucos
        buffer(7 ,i,j)=GEO_b
        buffer(8 ,i,j)=GEO_g_theta
        buffer(9 ,i,j)=GEO_grad_r
        buffer(10,i,j)=GEO_gq
        buffer(11,i,j)=GEO_captheta

     enddo ! j

     ! Set up nu for plotting and synthetic diagnostic.
     ! Note we need theta=-pi and pi so include 0 index
     do j=0,n_theta_plot
        theta = -pi+real(j)*pi_2/n_theta_plot
        if (n_theta_plot == 1) theta = 0.0 ! Test for special case
        call GEO_interp(theta)
        nu_coarse(j,i)=GEO_nu
     enddo
     do j=1,n_wedge
        theta=theta_wedge_offset+real(j-1)*theta_wedge_angle/       &
             real(n_theta_plot*n_theta_mult-1)
        if (n_wedge == 1) theta = 0.0 ! Test for special case
        call GEO_interp(theta)
        nu_wedge(j,i)=GEO_nu
     enddo

  enddo ! i


  call dump_h5(rootid,"nu",      buffer(1,:,:),h5in,h5err)
  call dump_h5(rootid,"gsin",    buffer(2,:,:),h5in,h5err)
  call dump_h5(rootid,"gcos1",   buffer(3,:,:),h5in,h5err)
  call dump_h5(rootid,"gcos2",   buffer(4,:,:),h5in,h5err)
  call dump_h5(rootid,"usin",    buffer(5,:,:),h5in,h5err)
  call dump_h5(rootid,"ucos",    buffer(6,:,:),h5in,h5err)
  call dump_h5(rootid,"b",       buffer(7,:,:),h5in,h5err)
  call dump_h5(rootid,"g_theta", buffer(8,:,:),h5in,h5err)
  call dump_h5(rootid,"grad_r",  buffer(9,:,:),h5in,h5err)
  call dump_h5(rootid,"gq",      buffer(10,:,:),h5in,h5err)
  call dump_h5(rootid,"captheta",buffer(11,:,:),h5in,h5err)
  deallocate(buffer)

  call close_h5file(fid,rootid,h5err)

end subroutine gyro_write_initdata_hdf5
