!-----------------------------------------------------------------
!  input routines
!-----------------------------------------------------------------
SUBROUTINE put_model_parameters(use_tm, adiabatic, alpha_p, &
           alpha_quench,version,lprint)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  LOGICAL,INTENT(IN) :: use_tm,adiabatic
  REAL,INTENT(IN) :: alpha_p,alpha_quench
  INTEGER,INTENT(IN) :: version,lprint
  !
  ! transfer values
  ! 
  use_transport_model_gf = use_tm
  use_adiabatic_electrons_gf = adiabatic
  alpha_p_mult_gf = alpha_p
  alpha_e_mult_gf = alpha_quench
  version_gf = version
  if(version.eq.1)then
     iglf = 0   !original glf23
  elseif(version.eq.2)then
     iglf = 1   !retuned version v1.61
  elseif(version.eq.3)then
     iglf = 98  !renormed + real geometry version 
  else
    write(*,*)"version specified does not exist",version,"reverting to version 3"
    iglf = 3
  endif
  lprint_gf = lprint
  !
END SUBROUTINE put_model_parameters
!
!-----------------------------------------------------------------
!
SUBROUTINE put_species(nsp,zsp,msp)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN):: nsp
  REAL,INTENT(IN) :: zsp(3),msp(3)
  INTEGER :: is
  !
! check for valid input
  if(nsp.lt.2.or.nsp.gt.3)call glf23_error(1,"input number of species invlaid")
  do is=1,nsp
    if(glf23_isnan(zsp(is)))call glf23_error(1,"input zs_in is NAN")
    if(glf23_isinf(zsp(is)))call glf23_error(1,"input zs_in is INF")
    if(glf23_isnan(msp(is)))call glf23_error(1,"input mass_in is NAN")
    if(glf23_isinf(msp(is)))call glf23_error(1,"input mass_in is INF")
    if(zsp(is).eq.0.0)call glf23_error(1,"input zs_in is = 0")
    if(msp(is).le.0.0)call glf23_error(1,"input mass_in is <= 0")
  enddo
  ! transfer values
  ns_gf = nsp
  amassgas_gf = msp(2)
  if(nsp.gt.2)then
    amassimp_gf = msp(3)
    zimp_gf = zsp(3)
  endif
  !
END SUBROUTINE put_species
!
!-----------------------------------------------------------------
!
SUBROUTINE put_kys(kys)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: kys

  if(glf23_isnan(kys))call glf23_error(1,"input ky_in is NAN")
  if(glf23_isinf(kys))call glf23_error(1,"input ky_in is INF")
  if(kys.le.0)call glf23_error(1,"input kys is <= 0")
  !
  ! transfer values
  !
  xky0_gf = kys
  !
END SUBROUTINE put_kys
!
!-----------------------------------------------------------------
!
SUBROUTINE put_gradients(rln,rlt,vpar_shear,vexb_shear)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: rln(3),rlt(3),vpar_shear(3)
  REAL,INTENT(IN) :: vexb_shear
  INTEGER :: is

  do is=1,nsm
   if(glf23_isnan(rln(is)))call glf23_error(1,"input rlns_in is NAN")
   if(glf23_isinf(rln(is)))call glf23_error(1,"input rlns_in is INF")
   if(glf23_isnan(rlt(is)))call glf23_error(1,"input rlts_in is NAN")
   if(glf23_isnan(rlt(is)))call glf23_error(1,"input rlts_in is INF")
   if(glf23_isnan(vpar_shear(is)))call glf23_error(1,"input vpar_shear_in is NAN")
   if(glf23_isnan(vpar_shear(is)))call glf23_error(1,"input vpar_shear_in is INF")
  enddo
  if(glf23_isnan(vexb_shear))call glf23_error(1,"input vexb_shear_in is NAN")
  if(glf23_isinf(vexb_shear))call glf23_error(1,"input vexb_shear_in is INF")

  !
  ! transfer values
  !
    rlte_gf = rlt(1)
    rlti_gf = rlt(2)
    rlne_gf = rln(1)
    rlni_gf = rln(2)
    rlnimp_gf = rln(3)
    gamma_p_gf = vpar_shear(2)
    gamma_e_gf = vexb_shear   
  !    
END SUBROUTINE put_gradients
!
!-----------------------------------------------------------------
!
SUBROUTINE put_averages(tsp,asp,betae,xnue)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: tsp(3),asp(3)
  REAL,INTENT(IN) :: betae,xnue
  INTEGER :: is

  do is=1,3
    if(glf23_isnan(tsp(is)))call glf23_error(1,"input taus_in is NAN")
    if(glf23_isinf(tsp(is)))call glf23_error(1,"input taus_in is INF")
    if(tsp(is).lt.0.0)call glf23_error(1,"input taus_in is < 0")
    if(glf23_isnan(asp(is)))call glf23_error(1,"input as_in is NAN")
    if(glf23_isinf(asp(is)))call glf23_error(1,"input as_in is INF")
    if(asp(is).lt.0.0)call glf23_error(1,"input as_in is < 0")
  enddo
    if(glf23_isnan(betae))call glf23_error(1,"input betae_in is NAN")
    if(glf23_isinf(betae))call glf23_error(1,"input betae_in is INF")
    if(betae.lt.0.0)call glf23_error(1,"input betae_in is < 0")
    if(glf23_isnan(xnue))call glf23_error(1,"input xnue_in is NAN")
    if(glf23_isinf(xnue))call glf23_error(1,"input xnue_in is INF")
    if(xnue.lt.0.0)call glf23_error(1,"input xnue_in is < 0")

  !
  ! transfer values
   taui_gf = tsp(2)
   apwt_gf = asp(2)
   aiwt_gf = asp(3)
   betae_gf = betae
   if(betae_gf.eq.0.0)betae_gf = 1.0E-12
   xnu_gf = xnue
  !      
END SUBROUTINE put_averages
!
!-----------------------------------------------------------------
!
SUBROUTINE put_s_alpha_geometry(rmin,rmaj,q,shat,alpha,xwell,theta0)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: rmin,rmaj,q,shat,alpha,theta0,xwell

  if(glf23_isnan(rmin))call glf23_error(1,"input rmin_sa is NAN")
  if(glf23_isinf(rmin))call glf23_error(1,"input rmin_sa is INF")
  if(glf23_isnan(rmaj))call glf23_error(1,"input rmaj_sa is NAN")
  if(glf23_isinf(rmaj))call glf23_error(1,"input rmaj_sa is INF")
  if(glf23_isnan(q))call glf23_error(1,"input q_sa is NAN")
  if(glf23_isinf(q))call glf23_error(1,"input q_sa is INF")
  if(glf23_isnan(shat))call glf23_error(1,"input shat_sa is NAN")
  if(glf23_isinf(shat))call glf23_error(1,"input shat_sa is INF")
  if(glf23_isnan(alpha))call glf23_error(1,"input alpha_sa is NAN")
  if(glf23_isinf(alpha))call glf23_error(1,"input alpha_sa is INF")
  if(glf23_isnan(xwell))call glf23_error(1,"input xwell_sa is NAN")
  if(glf23_isinf(xwell))call glf23_error(1,"input xwell_sa is INF")
  if(glf23_isnan(theta0))call glf23_error(1,"input theta0_sa is NAN")
  if(glf23_isinf(theta0))call glf23_error(1,"input theta0_sa is INF")
  !
  ! transfer values
  !
  rmin_gf = rmin
  rmaj_gf = rmaj
  q_gf = ABS(q)
  shat_gf = shat
  alpha_gf = alpha
  xwell_gf = xwell
  theta0_gf = theta0
  !
  ! validatiy checks
  !
  if(rmin_gf.ge.rmaj_gf)rmin_gf=0.999*rmaj_gf  
  !
END SUBROUTINE put_s_alpha_geometry
!
!-----------------------------------------------------------------
!  output routines
!-----------------------------------------------------------------
!
REAL FUNCTION get_growthrate(index1)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: index1
  INTEGER :: i3
  !
  i3 = SIZE(gamma_gf)
  if(index1.gt.i3)then
     write(*,*)"requested growthrate index out of bounds",i3
     get_growthrate = 0.0
  else
     get_growthrate = gamma_gf(index1)
  endif
  !
END FUNCTION get_growthrate
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_frequency(index1)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: index1
  INTEGER :: i3
  !
  i3 = SIZE(freq_gf)
  if(index1.gt.i3)then
     write(*,*)"requested frequency index is of bounds",i3
     get_frequency = 0.0
  else
     get_frequency = freq_gf(index1)
  endif
  !
END FUNCTION get_frequency
!
!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_particle_flux(i1,i2)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  ! note that only electrostatic part of the transport is computed in glf23
  !
   get_particle_flux = 0.0
  !
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested particle flux index is of bounds",i3
  elseif(i2.eq.1)then
     if(i1.eq.1)then
       get_particle_flux = rlne_gf*diff_gf
     elseif(i1.eq.2)then
       get_particle_flux = rlne_gf*diff_gf-rlnimp_gf*diff_im_gf
     elseif(i1.eq.3)then
       get_particle_flux = rlnimp_gf*diff_im_gf
     endif
  endif
  !
END FUNCTION get_particle_flux
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_energy_flux(i1,i2)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  ! note that only electrostatic part of the transport is computed in glf23
  !
  get_energy_flux = 0.0
  !
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested energy flux index is of bounds",i3
  elseif(i2.eq.1)then
     if(i1.eq.1)then
       get_energy_flux = rlte_gf*chie_gf
     elseif(i1.eq.2)then
       get_energy_flux = rlti_gf*chii_gf
     elseif(i1.eq.3)then ! impurity energy transport is in chii_gf
       get_energy_flux = 0.0
     endif
  endif
  !
END FUNCTION get_energy_flux
!-----------------------------------------------------------------
!
REAL FUNCTION get_stress_par(i1,i2)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  ! note that only electrostatic part of the transport is computed in glf23
  !
  get_stress_par = 0.0
  !
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested stress_par index is of bounds",i3
  elseif(i2.eq.1)then
     if(i1.eq.1)then
       get_stress_par = 0.0
     elseif(i1.eq.2)then
       get_stress_par = gamma_p_gf*eta_par_gf
     elseif(i1.eq.3)then 
       get_stress_par = 0.0
     endif
  endif
  !
END FUNCTION get_stress_par
!!-----------------------------------------------------------------
!
REAL FUNCTION get_stress_tor(i1,i2)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  ! note that only electrostatic part of the transport is computed in glf23
  !
  get_stress_tor = 0.0
  !
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested stress_tor index is of bounds",i3
  elseif(i2.eq.1)then
     if(i1.eq.1)then
       get_stress_tor = 0.0
     elseif(i1.eq.2)then
       get_stress_tor = gamma_p_gf*eta_phi_gf
     elseif(i1.eq.3)then 
       get_stress_tor = 0.0
     endif
  endif
  !
END FUNCTION get_stress_tor
!-----------------------------------------------------------------
!
REAL FUNCTION get_exchange(i1,i2)
  !
  USE glf23_gf
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  ! note that only electrostatic part of the transport is computed in glf23
  !
  get_exchange = 0.0
  !
  i3=6
  if(i1*i2.gt.i3)then
     write(*,*)"requested exchange index is of bounds",i3
  elseif(i2.eq.1)then
     if(i1.eq.1)then
       get_exchange = exch_gf
     elseif(i1.eq.2)then
       get_exchange = -exch_gf
     elseif(i1.eq.3)then 
       get_exchange = 0.0
     endif
  endif
  !
END FUNCTION get_exchange
!
!---------------------------------------------
!
SUBROUTINE glf23_error(status,message)
!
  use glf23_gf
  implicit none

  character (len=*) :: message
  integer :: status

  write(*,*)trim(message)
!
  if(status.eq.1)STOP
!
END SUBROUTINE glf23_error
