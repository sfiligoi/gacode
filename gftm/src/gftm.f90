!---------------------------------------------------------
! gftm_run.f90
!
! PURPOSE:
!  Manage standalone gftm run by calling
!
!    call gftm_read_input()
!    call gftm_run()
!
!  and then adding extra I/O. 
!---------------------------------------------------------

program gftm

  use gftm_pkg
  use gftm_interface
  use gftm_global
  use gftm_dimensions

  implicit none

  integer :: i
  integer :: n
  character (len=4) :: tag(6)=(/'ion1','ion2','ion3','ion4','ion5','ion6'/)
  real :: prec

  call gftm_read_input()
  call gftm_run()
  call gftm_dump_global()
  if(units_in.eq.'GENE')then
     print 30,'GENE reference units used'
     print 30,'Conversion to gftm units:'
     print 30,'Bunit/Bref = ',1.0/Bref_out
     print 30,'Te/Tref = ',taus_in(1)
     print 30,'mi/mref = ',mass_in(2)
     print 30,'a/Lref = ',1.0
     print 30,'cs/cref = ',SQRT(taus_in(1)/mass_in(2))
     print 30,'rhos/rhoref = ',SQRT(mass_in(2)*taus_in(1))*Bref_out
  endif
  if(alpha_zf_in.lt.0.0)print 30,' kx_geo0_out = ',kx_geo0_out, &
      ' SAT_geo0_out = ',SAT_geo0_out

! write interchange stability criteria with ELITE conventions
  Print 30,'  D(R) = ',-interchange_DR,'  D(I) = ',0.25-interchange_DM
!  write species info
  print 40,'  kinetic species = ',ns_in,'  non-kinetic species = ',nstotal_in - ns_in

  if (gftm_use_transport_model_in) then

     ! Output to screen

     print 20,'Gam/Gam_GB','    Q/Q_GB','Q_low/Q_GB','  Pi/Pi_GB', '    S/S_GB'
     print 10,'elec',&
          gftm_elec_pflux_out,&
          gftm_elec_eflux_out,&
          gftm_elec_eflux_low_out,&
          gftm_elec_mflux_out,&
          gftm_elec_expwd_out

     prec = abs(gftm_elec_pflux_out)+&
          abs(gftm_elec_eflux_out)+&
          abs(gftm_elec_eflux_low_out)+&
          abs(gftm_elec_mflux_out)

     do i=1,gftm_ns_in-1
        print 10,tag(i),&
             gftm_ion_pflux_out(i),&
             gftm_ion_eflux_out(i),&
             gftm_ion_eflux_low_out(i),&
             gftm_ion_mflux_out(i),&
             gftm_ion_expwd_out(i)

        prec = prec+&
             abs(gftm_ion_pflux_out(i))+&
             abs(gftm_ion_eflux_out(i))+&
             abs(gftm_ion_eflux_low_out(i))+&
             abs(gftm_ion_mflux_out(i))
     enddo

     ! Output to file

     n = gftm_ns_in-1
     open(unit=1,file='out.gftm.gbflux',status='replace')
     write(1,'(32(1pe11.4,1x))') gftm_elec_pflux_out,gftm_ion_pflux_out(1:n),&
          gftm_elec_eflux_out,gftm_ion_eflux_out(1:n),&
          gftm_elec_mflux_out,gftm_ion_mflux_out(1:n),&
          gftm_elec_expwd_out,gftm_ion_expwd_out(1:n)
     close(1)

     open(unit=1,file='out.gftm.grid',status='replace')
     write(1,'(i2)') gftm_ns_in,gftm_nxgrid_in
     close(1)

     write(*,*)"nbasis = ",nbasis,"  nu = ",nu,"  ne = ",ne

     ! write ky spectrum to file out.gftm.ky_spectrum
     CALL write_gftm_ky_spectrum

     ! write flux spectrum summed over nmodes to file out.gftm.sum_flux_spectrum
     ! this can be compared directly with the flux spectrum from CGYRO 

     CALL write_gftm_sum_flux_spectrum

     ! write density fluctuation amplitude spectrum to file out.gftm.density_spectrum
     CALL write_gftm_density_spectrum

     ! write temperature fluctuation amplitude spectrum to file out.gftm.temperature_spectrum
     CALL write_gftm_temperature_spectrum

     ! write intensity fluctuation amplitude spectrum per mode to file out.gftm.intensity_spectrum
     CALL write_gftm_intensity_spectrum

     ! write field fluctuation amplitude per mode spectrum to file out.gftm.field_spectrum
     CALL write_gftm_field_spectrum

     ! write eigenvalue spectrum to file out.gftm.eigenvalue_spectrum
     CALL write_gftm_eigenvalue_spectrum

     ! write ne-te crossphase spectrum to file out.gftm.nete_crossphase_spectrum
     CALL write_gftm_nete_crossphase_spectrum

     ! write ns-ts crossphase spectrum to file out.gftm.nsts_crossphase_spectrum
     CALL write_gftm_nsts_crossphase_spectrum

     ! write QL flux (weight) spectrum per mode to file out.gftm.QL_weight_spectrum
     CALL write_gftm_QL_flux_spectrum

     ! write kx/ky-spectral shift spectrum to file out.gftm.spectral_shift
     CALL write_gftm_spectral_shift_spectrum

     ! write ave_p0 spectrum to file out.gftm.ave_p0_spectrum
     CALL write_gftm_ave_p0_spectrum

     ! write intensity fluctuation amplitude spectrum per mode to file out.gftm.scalar_saturation_parameters
     CALL write_gftm_scalar_saturation_parameters

     ! write Gaussian width spectrum to file out.gftm.width_spectrum
     CALL write_gftm_width_spectrum

  else

     print 10,'     ky:',gftm_ky_in
     print 10,'     ft:',get_ft()
     print 10,'Guassian width = ',get_gaussian_width()
     ! Collect linear eigenvalues
     do i=1,gftm_nmodes_in
        print 10,'(wr,wi):',gftm_eigenvalue_out(i)
     enddo

     prec = sum(abs(gftm_eigenvalue_out))

    ! write single ky-eigenmode wavefunction to file
    CALL write_wavefunction_out('out.gftm.wavefunction')

  endif

  open(unit=1,file=trim(gftm_path_in)//'out.gftm.prec')
  write(1,*) prec
  close(1)

10 format(a,10(1x,1pe11.4))
20 format(t7,a,t19,a,t31,a,t43,a,t55,a)
30 format(a,1pe11.4,a,1pe11.4)
40 format(a,I10,a,I10)

end program gftm
