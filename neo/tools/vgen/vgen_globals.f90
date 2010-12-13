module vgen_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err
  
  ! vgen inputs
  integer :: er_method
  integer :: erspecies_indx
  integer :: vel_method
  
  real, parameter :: pi=3.1415926535897932
  real, parameter :: mass_deuterium = 3.3452   ! (x 10-27 kg)
  real, parameter :: temp_norm_fac   = 1602.2
  real, parameter :: charge_norm_fac = 1.6022
  
  real :: dens_norm, temp_norm, mass_norm, vth_norm
  
end module vgen_globals
