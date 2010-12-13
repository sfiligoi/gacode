!-----------------------------------------------------
! write_units.f90
!
! PURPOSE:
!  Write normalizing parameters.
!-----------------------------------------------------

subroutine gyro_write_units(datafile,io)

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  real :: kt
  !---------------------------------------------------

  select case (output_flag)

  case (1)

     open(unit=io,file=datafile,status='replace')

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

  end select

  if (debug_flag == 1) then
     print *,'[gyro_write_units called]'
  endif

26 format(t2,1pe15.8,2x,a) 

end subroutine gyro_write_units
