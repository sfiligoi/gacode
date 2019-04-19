!----------------------------------------------------------------
! expro_locsim_profiles.f90
!
! PURPOSE:
!  Read experimental profiles and generate local profile 
!  parameters at rmin = r/a.
!
! INPUTS:
!  path              : path to data
!  comm              : MPI communicator
!  numeq_flag        : Fourier series equilibrium (0=no,1=yes)
!  udsymmetry_flag   : enforce up-down symmetry (0=no,1=yes)
!  quasineutral_flag : enforce quasineutrality (0=no,1=yes)
!  n_species_in      : total species (e+i)
!  z                 : vector of charges (length n_species_in-1)
!  rmin              : r/a
!
! OUTPUTS:
!  btccw             : (+1 or -1) 
!  ipccw             : (+1 or -1)
!  a_meters          : minor radius i metres
!
! OUTPUTS (interface):
!  real :: shift_loc
!  real :: q_loc
!  real :: s_loc
!  real :: kappa_loc
!  real :: delta_loc
!  real :: zeta_loc
!  real :: s_kappa_loc
!  real :: s_delta_loc
!  real :: s_zeta_loc
!  real :: zmag_loc
!  real :: dzmag_loc
!  real :: gamma_e_loc
!  real :: gamma_p_loc
!  real :: mach_loc
!  real :: rmaj_loc
!  real :: rhos_loc [m]
!  real :: z_eff_loc
!  real :: b_unit_loc
!  real :: rho_norm_loc
!  real :: psi_norm_loc
!  real :: psi_a_loc
!  real :: cs_loc
!  real :: beta_star_loc
!
!  real, dimension(9) :: dens_loc
!  real, dimension(9) :: temp_loc
!  real, dimension(9) :: dlnndr_loc
!  real, dimension(9) :: dlntdr_loc
!  real, dimension(9) :: sdlnndr_loc
!  real, dimension(9) :: sdlntdr_loc
!----------------------------------------------------------------

subroutine expro_locsim_profiles(&
     path,&
     comm,&
     numeq_flag,&
     udsymmetry_flag,&
     quasineutral_flag,&
     n_species_in,&
     rmin,&
     btccw,&
     ipccw,&
     a_meters)

  use vpro
  use expro_locsim_interface

  implicit none

  character(len=*), intent(in) :: path 
  integer, intent(in) :: comm
  integer, intent(in) :: numeq_flag
  integer, intent(in) :: udsymmetry_flag
  integer, intent(in) :: quasineutral_flag
  integer, intent(in) :: n_species_in
  double precision, intent(in) :: rmin
  double precision, intent(inout) :: btccw,ipccw,a_meters
  double precision, parameter :: pi=3.14159265358979323846
  integer :: i,j,i_ion

  rmin_loc = rmin
  
  n_species_exp = n_species_in

  !--------------------------------------------------------------
  ! use expro routines to read data:
  !
  expro_ctrl_quasineutral_flag = 1  ! quasi-neutrality density flag
  expro_ctrl_numeq_flag = numeq_flag
  expro_ctrl_n_ion = n_species_exp-1

  !if (comm == -1) then
     call vpro_read
  !else
  !   call expro_palloc(comm,path,1)
  !   call expro_pread
  !endif
  call expro_locsim_alloc(1)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !
  btccw = -expro_signb
  ipccw = -expro_signq*expro_signb

  rmin_exp(:) = expro_rmin(:)

  if (udsymmetry_flag == 1) then
     expro_zmag(:) = 0.0   
     expro_dzmag(:) = 0.0
  endif

  ! Minor radius, a, in meters:
  a_meters = rmin_exp(expro_n_exp)

  rmin_exp(:) = rmin_exp(:)/a_meters

  ! Pack electrons into top of species vector.
  temp_exp(n_species_exp,:)    = expro_te(:)
  dlntdr_exp(n_species_exp,:)  = expro_dlntedr(:)*a_meters 
  sdlntdr_exp(n_species_exp,:) = expro_sdlntedr(:)*a_meters
  dens_exp(n_species_exp,:)    = expro_ne(:)
  dlnndr_exp(n_species_exp,:)  = expro_dlnnedr(:)*a_meters 
  sdlnndr_exp(n_species_exp,:) = expro_sdlnnedr(:)*a_meters

  mass_loc(n_species_exp) = 1.0/1837
  z_loc(n_species_exp) = -1.0
 
  ! Pack ions from the bottom
  do i_ion=1,n_species_exp-1
     ! ion temps should be equal, but not enforced 
     temp_exp(i_ion,:)    = expro_ti(i_ion,:)
     dlntdr_exp(i_ion,:)  = expro_dlntidr(i_ion,:)*a_meters 
     sdlntdr_exp(i_ion,:) = expro_sdlntidr(i_ion,:)*a_meters 

     mass_loc(i_ion) = expro_mass(i_ion)
     z_loc(i_ion) = expro_z(i_ion)
   
     ! First species density is reset by quasi-neutrality
     if (quasineutral_flag == 1 .and. i_ion == 1) then
        dens_exp(i_ion,:)    = expro_ni_new(:)
        dlnndr_exp(i_ion,:)  = expro_dlnnidr_new(:)*a_meters
        sdlnndr_exp(i_ion,:) = expro_sdlnnidr_new(:)*a_meters
     else
        dens_exp(i_ion,:)    = expro_ni(i_ion,:)
        dlnndr_exp(i_ion,:)  = expro_dlnnidr(i_ion,:)*a_meters
        sdlnndr_exp(i_ion,:) = expro_sdlnnidr(i_ion,:)*a_meters
     endif
  enddo

  ! Sanity check for densities
  do i=1,n_species_exp
     if (minval(dens_exp(i,:)) <= 0.0) then
        !call cgyro_error('Nonpositive in exp. density profile')
        return
     endif
  enddo

  ! Rotation
  gamma_e_exp(:) = -expro_w0p(:)*(a_meters*rmin_exp(:))/expro_q(:)
  gamma_p_exp(:) = -expro_w0p(:)*expro_rmaj(:)
  mach_exp(:)    = expro_w0(:)*expro_rmaj(:)

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_exp) ones.
  ! 
  call cub_spline(rmin_exp,expro_rmaj/a_meters,expro_n_exp,rmin,rmaj_loc,1)
  call cub_spline(rmin_exp,expro_q,expro_n_exp,rmin,q_loc,1)
  call cub_spline(rmin_exp,expro_s,expro_n_exp,rmin,s_loc,1)
  call cub_spline(rmin_exp,expro_drmaj,expro_n_exp,rmin,shift_loc,1)
  call cub_spline(rmin_exp,expro_kappa,expro_n_exp,rmin,kappa_loc,1)
  call cub_spline(rmin_exp,expro_skappa,expro_n_exp,rmin,s_kappa_loc,1)
  call cub_spline(rmin_exp,expro_delta,expro_n_exp,rmin,delta_loc,1)
  call cub_spline(rmin_exp,expro_sdelta,expro_n_exp,rmin,s_delta_loc,1)
  call cub_spline(rmin_exp,expro_zeta,expro_n_exp,rmin,zeta_loc,1)
  call cub_spline(rmin_exp,expro_szeta,expro_n_exp,rmin,s_zeta_loc,1)
  call cub_spline(rmin_exp,expro_zmag/a_meters,expro_n_exp,rmin,zmag_loc,1)
  call cub_spline(rmin_exp,expro_dzmag,expro_n_exp,rmin,dzmag_loc,1)
  call cub_spline(rmin_exp,gamma_e_exp,expro_n_exp,rmin,gamma_e_loc,1)
  call cub_spline(rmin_exp,gamma_p_exp,expro_n_exp,rmin,gamma_p_loc,1)
  call cub_spline(rmin_exp,mach_exp,expro_n_exp,rmin,mach_loc,1)
  call cub_spline(rmin_exp,expro_rhos,expro_n_exp,rmin,rhos_loc,1)
  call cub_spline(rmin_exp,expro_cs,expro_n_exp,rmin,cs_loc,1)
  call cub_spline(rmin_exp,expro_z_eff,expro_n_exp,rmin,z_eff_loc,1)
  call cub_spline(rmin_exp,expro_bunit,expro_n_exp,rmin,b_unit_loc,1)
  call cub_spline(rmin_exp,expro_rho,expro_n_exp,rmin,rho_norm_loc,1)
  call cub_spline(rmin_exp,expro_polflux,expro_n_exp,rmin,psi_norm_loc,1)
  psi_norm_loc = psi_norm_loc/expro_polflux(expro_n_exp)
  psi_a_loc = expro_polflux(expro_n_exp)


  beta_star_loc = 0.0  
  do i=1,n_species_exp
     ! Note: mapping is only done for n_species (not n_species_exp)
     call cub_spline(rmin_exp,dens_exp(i,:),expro_n_exp,rmin,dens_loc(i),1)
     call cub_spline(rmin_exp,temp_exp(i,:),expro_n_exp,rmin,temp_loc(i),1)
     call cub_spline(rmin_exp,dlntdr_exp(i,:),expro_n_exp,rmin,dlntdr_loc(i),1)
     call cub_spline(rmin_exp,dlnndr_exp(i,:),expro_n_exp,rmin,dlnndr_loc(i),1)
     call cub_spline(rmin_exp,sdlntdr_exp(i,:),expro_n_exp,rmin,sdlntdr_loc(i),1)
     call cub_spline(rmin_exp,sdlnndr_exp(i,:),expro_n_exp,rmin,sdlnndr_loc(i),1)
     beta_star_loc = beta_star_loc+dens_loc(i)*temp_loc(i)*(dlnndr_loc(i)+dlntdr_loc(i))
  enddo
  ! CGS beta calculation
  betae_loc = 4.027e-3*dens_loc(n_species_exp)*temp_loc(n_species_exp)/b_unit_loc**2

  beta_star_loc = beta_star_loc*betae_loc/(dens_loc(n_species_exp)*temp_loc(n_species_exp))
  
  if (numeq_flag == 1) then

     if (expro_nfourier <= 0) then
        return
     endif

     geo_ny_loc = expro_nfourier
     allocate(geo_yin_exp(8,0:geo_ny_loc,expro_n_exp))
     if(allocated(geo_yin_loc)) deallocate(geo_yin_loc)
     allocate(geo_yin_loc(8,0:geo_ny_loc))
     geo_yin_exp(1:4,:,:) = expro_geo(:,:,:)/a_meters
     geo_yin_exp(5:8,:,:) = expro_dgeo(:,:,:)

     do i=1,8
        do j=0,geo_ny_loc
           call cub_spline(rmin_exp,geo_yin_exp(i,j,:),expro_n_exp,rmin, &
                geo_yin_loc(i,j),1)
        enddo
     enddo

  endif

end subroutine expro_locsim_profiles
