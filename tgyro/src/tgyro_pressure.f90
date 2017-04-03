!-------------------------------------------------------------
! tgyro_pressure.f90
!
! PURPOSE:
!  Localize the calculation of total pressure, as well as 
!  effective temperature needed by TOQ/EPED.
!--------------------------------------------------------------

subroutine tgyro_pressure

  use tgyro_globals
  use tgyro_ped
  
  implicit none
  integer :: i_ion
 
  ! First compute thermal pressure
  ptot_exp = exp_ne*exp_te
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        ptot_exp = ptot_exp + exp_ni(i_ion,:)*exp_ti(i_ion,:)
     endif
  enddo

  ! Mean thermal temperature on axis.  Need this for input 
  ! in tgyro_eped_nn
  !
  !            p_thermal
  ! Te_TOQ = ------------
  !            n_thermal

  te_toq = ptot_exp(1)/(exp_ne(1)+sum(exp_ni(therm_vec,1)))

  ! Now add fast ions
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 0) then
        ptot_exp = ptot_exp + exp_ni(i_ion,:)*exp_ti(i_ion,:)
     endif
  enddo

  ! Convert to Pa: n[1/cm^3]*(kT[ev])/10  
  ptot_exp = ptot_exp*k/10.0

end subroutine tgyro_pressure
