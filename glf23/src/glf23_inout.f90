!-----------------------------------------------------------------
!  input routines
!-----------------------------------------------------------------

subroutine put_glf23_model_parameters(use_tm, adiabatic, alpha_p, &
     alpha_quench,version,lprint)

  use glf23_gf

  implicit none
  logical,intent(in) :: use_tm,adiabatic
  real,intent(in) :: alpha_p,alpha_quench
  integer,intent(in) :: version,lprint

  ! transfer values
  use_transport_model_gf = use_tm
  use_adiabatic_electrons_gf = adiabatic
  alpha_p_mult_gf = alpha_p
  alpha_e_mult_gf = alpha_quench

  if (version < 1 .or. version > 3) then
     write(*,*)"version specified does not exist",version,"reverting to version 2"
     version_gf = 2
  else  
     version_gf = version
  endif
  lprint_gf = lprint

end subroutine put_glf23_model_parameters

!-----------------------------------------------------------------------

subroutine put_glf23_species(nsp,zsp,msp)

  use glf23_gf

  implicit none
  integer,intent(in):: nsp
  real,intent(in) :: zsp(3),msp(3)
  integer :: is

  ! check for valid input
  if (nsp < 2 .or. nsp > 3) call glf23_error(1,"input number of species invlaid")
  do is=1,nsp
     if (zsp(is) == 0.0) call glf23_error(1,"input zs_in is = 0")
     if (msp(is) <= 0.0) call glf23_error(1,"input mass_in is <= 0")
  enddo

  ! transfer values
  ns_gf = nsp
  amassgas_gf = msp(2)
  if (nsp > 2) then
     amassimp_gf = msp(3)
     zimp_gf = zsp(3)
  endif

end subroutine put_glf23_species

!-----------------------------------------------------------------

subroutine put_glf23_kys(kys)

  use glf23_gf

  implicit none
  real, intent(in) :: kys

  if (kys <= 0.0) call glf23_error(1,"input kys is <= 0")

  ! transfer values
  xky0_gf = kys

end subroutine put_glf23_kys

!-----------------------------------------------------------------

subroutine put_glf23_gradients(rln,rlt,vpar_shear,vexb_shear)

  use glf23_gf

  implicit none
  real, intent(in) :: rln(3),rlt(3),vpar_shear(3)
  real, intent(in) :: vexb_shear
  integer :: is

  ! transfer values
  rlte_gf = rlt(1)
  rlti_gf = rlt(2)
  rlne_gf = rln(1)
  rlni_gf = rln(2)
  rlnimp_gf = rln(3)
  gamma_p_gf = vpar_shear(2)
  gamma_e_gf = vexb_shear   

end subroutine put_glf23_gradients

!-----------------------------------------------------------------

subroutine put_glf23_averages(tsp,asp,betae,xnue)
  
  use glf23_gf
  
  implicit none
  real,intent(in) :: tsp(3),asp(3)
  real,intent(in) :: betae,xnue
  integer :: is

  do is=1,3
     if(tsp(is) < 0.0) call glf23_error(1,"input taus_in is < 0")
     if(asp(is) < 0.0) call glf23_error(1,"input as_in is < 0")
  enddo
  if (betae < 0.0) call glf23_error(1,"input betae_in is < 0")
  if (xnue  < 0.0) call glf23_error(1,"input xnue_in is < 0")

  ! transfer values
  taui_gf = tsp(2)
  apwt_gf = asp(2)
  aiwt_gf = asp(3)
  betae_gf = betae

  if (betae_gf == 0.0) betae_gf = 1e-12

  xnu_gf = xnue
        
end subroutine put_glf23_averages

!-----------------------------------------------------------------

subroutine put_glf23_geometry(rmin,rmaj,q,shat,alpha,xwell,theta0)
  
  use glf23_gf
  
  implicit none
  real,intent(in) :: rmin,rmaj,q,shat,alpha,theta0,xwell

  ! transfer values
  rmin_gf = rmin
  rmaj_gf = rmaj
  q_gf = abs(q)
  shat_gf = shat
  alpha_gf = alpha
  xwell_gf = xwell
  theta0_gf = theta0
  
  ! validity checks
  if (rmin_gf >= rmaj_gf) rmin_gf = 0.999*rmaj_gf  
  
end subroutine put_glf23_geometry

!-----------------------------------------------------------------
!  output routines
!-----------------------------------------------------------------

real function get_glf23_growthrate(index1)
  
  use glf23_gf
  
  implicit none
  integer,intent(in) :: index1
  integer :: i3
  
  i3 = size(gamma_gf)
  if (index1 > i3)then
     write(*,*)"requested growthrate index out of bounds",i3
     get_glf23_growthrate = 0.0
  else
     get_glf23_growthrate = gamma_gf(index1)
  endif
  
end function get_glf23_growthrate

!-----------------------------------------------------------------

real function get_glf23_frequency(index1)
  
  use glf23_gf
  
  implicit none
  integer,intent(in) :: index1
  integer :: i3
  
  i3 = size(freq_gf)
  if (index1 > i3)then
     write(*,*)"requested frequency index is of bounds",i3
     get_glf23_frequency = 0.0
  else
     get_glf23_frequency = freq_gf(index1)
  endif
  
end function get_glf23_frequency

!-----------------------------------------------------------------

real function get_glf23_particle_flux(i1,i2)
  
  use glf23_gf
  
  implicit none
  integer,intent(in) :: i1,i2
  integer :: i3
  
  ! note that only electrostatic part of the transport is computed in glf23
  get_glf23_particle_flux = 0.0
  
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested particle flux index is of bounds",i3
  elseif(i2 == 1)then
     if(i1 == 1)then
        get_glf23_particle_flux = diff_gf+zimp_gf*diff_im_gf
     elseif(i1 == 2)then
        get_glf23_particle_flux = diff_gf
     elseif(i1 == 3)then
        get_glf23_particle_flux = diff_im_gf
     endif
  endif

end function get_glf23_particle_flux

!-----------------------------------------------------------------

real function get_glf23_energy_flux(i1,i2)
  !
  use glf23_gf
  !
  implicit none
  integer,intent(in) :: i1,i2
  integer :: i3
  !
  ! note that only electrostatic part of the transport is computed in glf23
  !
  get_glf23_energy_flux = 0.0
  !
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested energy flux index is of bounds",i3
  elseif(i2 == 1)then
     if(i1 == 1)then
        get_glf23_energy_flux = chie_gf
     elseif(i1 == 2)then
        get_glf23_energy_flux = chii_gf
     elseif(i1 == 3)then 
        get_glf23_energy_flux = chi_im_gf
     endif
  endif
  !
end function get_glf23_energy_flux

!-----------------------------------------------------------------

real function get_glf23_stress_par(i1,i2)

  use glf23_gf

  implicit none
  integer,intent(in) :: i1,i2
  integer :: i3

  ! note that only electrostatic part of the transport is computed in glf23
  get_glf23_stress_par = 0.0

  i3=6
  if (i1*i2 > i3)then
     write(*,*)"requested stress_par index is of bounds",i3
  else if(i2 == 1)then
     if (i1 == 1)then
        get_glf23_stress_par = 0.0
     else if (i1 == 2)then
        get_glf23_stress_par = eta_par_gf
     else if (i1 == 3)then 
        get_glf23_stress_par = 0.0
     endif
  endif

end function get_glf23_stress_par

!-----------------------------------------------------------------

real function get_glf23_stress_tor(i1,i2)
  
  use glf23_gf
  
  implicit none
  integer,intent(in) :: i1,i2
  integer :: i3
 
  ! note that only electrostatic part of the transport is computed in glf23
  get_glf23_stress_tor = 0.0
  
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested stress_tor index is of bounds",i3
  elseif(i2 == 1)then
     if(i1 == 1)then
        get_glf23_stress_tor = 0.0
     elseif(i1 == 2)then
        get_glf23_stress_tor = eta_phi_gf
     elseif(i1 == 3)then 
        get_glf23_stress_tor = 0.0
     endif
  endif
  
end function get_glf23_stress_tor

!-----------------------------------------------------------------

real function get_glf23_exchange(i1,i2)
  
  use glf23_gf
  
  implicit none
  integer, intent(in) :: i1,i2
  integer :: i3
  
  ! note that only electrostatic part of the transport is computed in glf23  
  get_glf23_exchange = 0.0
  
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested exchange index is of bounds",i3
  elseif(i2 == 1)then
     if(i1 == 1)then
        get_glf23_exchange = exch_gf
     elseif(i1 == 2)then
        get_glf23_exchange = -exch_gf
     elseif(i1 == 3)then 
        get_glf23_exchange = 0.0
     endif
  endif
  
end function get_glf23_exchange

!---------------------------------------------

subroutine glf23_error(status,message)

  use glf23_gf
  implicit none

  character (len=*) :: message
  integer :: status

  write(*,*) trim(message)

  if (status == 1) stop

end subroutine glf23_error
