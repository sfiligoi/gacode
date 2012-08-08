!-----------------------------------------------------------
! gyro_write_initdata.f90
!
! PURPOSE:
!  Consolidated code to write all initial data.  
!  file1: profile data
!  file2: units
!  file3: geometry data
!
!-----------------------------------------------------------

subroutine gyro_write_initdata(datafile1,datafile2,datafile3,io,h5datafile)

  use gyro_globals
  use gyro_profile_exp
  use math_constants
  use GEO_interface
#ifdef HAVE_HDF5
  use hdf5_api
#endif

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile1
  character (len=*), intent(in) :: datafile2
  character (len=*), intent(in) :: datafile3
  character (len=*), intent(in) :: h5datafile 
  !
  integer :: n_wedge
  real :: kt
  real :: theta
  real :: dr

#ifdef HAVE_HDF5
   include 'mpif.h'
  integer(HID_T) :: fid, rootid
  character(90) :: description
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  real :: buff
  real, allocatable :: buffer(:,:,:)
#endif

  !---------------------------------------------------

  if (output_flag == 0) return

  if (i_proc /= 0) return

  if(io_method < 3) then
     open(unit=io,file=datafile1,status='replace')

     ! Essential scalars
     write(io,20) n_x
     write(io,20) n_theta_section
     write(io,20) n_pass
     write(io,20) n_trap
     write(io,20) n_energy
     write(io,20) n_theta_plot
     write(io,20) n0
     write(io,20) n_n
     write(io,20) d_n
     write(io,20) n_explicit_damp
     write(io,20) nonlinear_flag
     write(io,20) electron_method
     write(io,20) n_field
     write(io,20) n_ion
     write(io,20) n_kinetic
     write(io,20) n_spec
     write(io,20) field_r0_flag
     write(io,20) field_r0_grid
     write(io,20) n_grid_exp
     write(io,20) boundary_method

     ! Basic profile data (** CORRECT IN VUGYRO **)
     write(io,10) r(:)
     write(io,10) q(:)
     write(io,10) r_s(:)
     write(io,10) q_s(:)
     write(io,10) dlntdr_s(:,:)
     write(io,10) dlnndr_s(:,:)
     write(io,10) tem_s(:,:)
     write(io,10) den_s(:,:)
     write(io,10) rmaj_s(:)/r_s(:)
     write(io,10) delta_s(:)
     write(io,10) zeta_s(:)
     write(io,10) kappa_s(:)
     write(io,10) drmaj_s(:)
     write(io,10) shat_s(:)
     write(io,10) s_delta_s(:)
     write(io,10) s_zeta_s(:)
     write(io,10) s_kappa_s(:)
     write(io,10) zmag_s(:)
     write(io,10) dzmag_s(:)
     write(io,10) beta_unit_s(:)
     write(io,10) gamma_e_s(:)
     write(io,10) gamma_p_s(:) 
     write(io,10) mach_s(:) 
     write(io,10) b_unit_s(:)
     write(io,10) dr_eodr(:)
     write(io,10) z_eff_s(:)
     write(io,10) nu_s(n_spec,:)
     write(io,10) w0_s(:)

     write(io,10) box_multiplier

     write(io,10) lambda(ir_norm,:)
     write(io,10) energy(:,1)
     write(io,10) lambda_tp(ir_norm)
     write(io,10) krho_collect(:)
     write(io,10) rhos_norm
     write(io,10) z(:)

     ! Recent additions
     write(io,20) n_theta_plot*n_theta_mult

     ! Number of flux moments
     write(io,20) n_moment

     close(io)

     !================================================================

     open(unit=io,file=datafile2,status='replace')

     ! kT in MJ (note the conversion 1.6022e-22 MJ/keV)
     kt = 1.6022e-22*tem_norm

     write(io,26) 2.0*kg_proton,'m_ref (kg)'
     write(io,26) b_unit_norm,'b_unit (Tesla)'
     write(io,26) a_meters,'a (m)'
     write(io,26) csda_norm,'csD/a (1/s)'
     write(io,26) csda_norm*a_meters,'csD (m/s)'
     write(io,26) tem_norm,'Te (keV)'
     write(io,26) den_norm,'ne (10^19/m^3)'
     write(io,26) rhos_norm*a_meters,'rho_sD (m)'
     write(io,26) csda_norm*(rhos_norm*a_meters)**2,&
          'chi_gBD (m^2/s)'

     write(io,26) 1e19*den_norm*(csda_norm*a_meters)*rhos_norm**2/0.624e22,&
          'Gamma_gBD (0.624e22/m^2/s) = (MW/keV/m^2)'

     write(io,26) 1e19*den_norm*(csda_norm*a_meters)*kt*rhos_norm**2,&
          'Q_gBD (MJ/m^2/s) = (MW/m^2)'

     write(io,26) 1e19*den_norm*a_meters*kt*rhos_norm**2*1e6,&
          'Pi_gBD (J/m^2) = (Nm/m^2)'

     write(io,26) 1e19*den_norm*csda_norm*kt*rhos_norm**2,&
          'S_gBD (MJ/m^3/s) = (MW/m^3)'

     close(io)

    endif !io_method <3
     !================================================================

    if(io_method < 3) then
      open(unit=io,file=datafile3,status='replace')
    endif !io_method <3


     n_wedge = n_theta_plot*n_theta_mult
              
#ifdef HAVE_HDF5
  if(io_method > 1) allocate(buffer(14,1:n_x,n_wedge))
#endif

     do i=1,n_x

        if (flat_profile_flag == 0) then

           ! All profiles are global and radial variation 
           ! is consistent

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

           ! Profiles are flat and so some parameters need
           ! to be linearly extrapolated.

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

#ifdef HAVE_HDF5
            if(io_method > 1) then
                
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
           endif
#endif

           write(io,10) GEO_nu
           write(io,10) GEO_gsin
           write(io,10) GEO_gcos1
           write(io,10) GEO_gcos2
           write(io,10) GEO_usin
           write(io,10) GEO_ucos
           write(io,10) GEO_b
           write(io,10) GEO_g_theta
           write(io,10) GEO_grad_r
           write(io,10) GEO_gq
           write(io,10) GEO_captheta

        enddo ! j
#ifdef HAVE_HDF5
     if(io_method > 1 ) then
       ! Set up nu for plotting and synthetic diagnostic.
       ! Note we need theta=-pi and pi so include 0 index
       do j=0,n_theta_plot
          theta = -pi+real(j)*pi_2/n_theta_plot
          if (n_theta_plot == 1) theta = 0.0 ! Test for special case
          call GEO_interp(theta)
          nu_coarse(j,i)=GEO_nu
       enddo
         do j=1,n_wedge
            if (n_wedge == 1) then
              theta = 0.0 ! Test for special case
            else
            theta=theta_wedge_offset+real(j-1)*theta_wedge_angle/       &
                 real(n_theta_plot*n_theta_mult-1)
            endif
            call GEO_interp(theta)
            nu_wedge(j,i)=GEO_nu
         enddo
      endif
#endif
     enddo ! i

    if(io_method <3) then
     close(io)
    endif



!hdf5 output handling
#ifdef HAVE_HDF5
  if(io_method >1) then
  !------------------------------------------
  ! Do the initialization here.  
  !------------------------------------------
    call vshdf5_fcinit()
    call vshdf5_inith5vars(h5in, h5err)
    h5in%comm=MPI_COMM_SELF
    h5in%info=MPI_INFO_NULL
    h5in%typeConvert=.true.
    h5in%wrd_type=H5T_NATIVE_REAL
    h5in%doTranspose=.true.
    h5in%wrVsTime=.false.
    h5in%verbose=.true.
    !---------------------------------------------------------------------
    ! Write the variables to an hdf5 file
    ! These variables are essentially the write_profile_vugyro.f90 
    !---------------------------------------------------------------------
    description=" "
    call open_newh5file(h5datafile,fid,description,rootid,h5in,h5err)
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
    call dump_h5(rootid,"n_grid_exp",n_grid_exp,h5in,h5err)
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
    call dump_h5(rootid,"n_explicit_damp",n_explicit_damp,h5in,h5err)

    call dump_h5(rootid,"lambda", lambda(ir_norm,:),h5in,h5err)
    call dump_h5(rootid,"energy", energy,h5in,h5err)
    call dump_h5(rootid,"lambda_tp", lambda_tp(ir_norm),h5in,h5err)
    call dump_h5(rootid,"krho_collect", krho_collect(:),h5in,h5err)
    call dump_h5(rootid,"rhos_norm", rhos_norm,h5in,h5err)
    call dump_h5(rootid,"zcharge", z(:),h5in,h5err)
    call dump_h5(rootid,"n_moment", n_moment ,h5in,h5err)


    call dump_h5(rootid,"time_skip", time_skip,h5in,h5err)
    call dump_h5(rootid,"time_step", dt ,h5in,h5err)
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
    call dump_h5(rootid,"nu",      buffer(1,:,:),h5in,h5err)
    call dump_h5(rootid,"gsin",    buffer(2,:,:),h5in,h5err)
    call dump_h5(rootid,"gcos1",   buffer(3,:,:),h5in,h5err)
    call dump_h5(rootid,"gcos2",   buffer(4,:,:),h5in,h5err)
    call dump_h5(rootid,"b",       buffer(7,:,:),h5in,h5err)
    call dump_h5(rootid,"g_theta", buffer(8,:,:),h5in,h5err)
    call dump_h5(rootid,"grad_r",  buffer(9,:,:),h5in,h5err)
    call dump_h5(rootid,"gq",      buffer(10,:,:),h5in,h5err)
    call dump_h5(rootid,"captheta",buffer(11,:,:),h5in,h5err)

    call close_h5file(fid,rootid,h5err)
    deallocate(buffer)
  endif !io_method >1 
#endif 
  if (debug_flag == 1 ) then
     print *,'[gyro_write_initdata done]'
  endif

10 format(1pe15.8)
20 format(i4)
26 format(t2,1pe15.8,2x,a) 

!  endif !i_proc == 0 

end subroutine gyro_write_initdata
