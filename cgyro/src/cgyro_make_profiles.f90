subroutine cgyro_make_profiles

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is,ir,ix 
  integer :: j
  integer :: num_ele

  !-------------------------------------------------------------
  ! Manage electrons
  !
  num_ele = 0
  do is=1,n_species
     if (z(is) == -1) then
        num_ele = num_ele + 1
        is_ele = is
     endif
  enddo

  if (num_ele == 0) then

     ! Adiabatic electrons

     ae_flag = 1
     call cgyro_info('Using adiabatic electron model.')

     dens_ele = ne_ade
     temp_ele = te_ade
     mass_ele = masse_ade

  else if (num_ele == 1) then

     ! GK electrons

     ae_flag = 0
     call cgyro_info('Using gyrokinetic electrons.')

     dens_ele = dens(is_ele)
     temp_ele = temp(is_ele)
     mass_ele = mass(is_ele)

  else

     call cgyro_error('Only one electron species allowed.')
     return

  endif
  !-------------------------------------------------------------


  !-------------------------------------------------------------
  ! Manage simulation type (n=0,linear,nonlinear)
  !
  q = abs(q) 

  if (zf_test_flag == 1) then

     ! Zonal flow (n=0) test

     k_theta = q/rmin
     rho     = ky/k_theta
     length  = box_size/(s*k_theta)

     k_theta = 0

     call cgyro_info('Triggered zonal flow test.')

     if (n_radial /= 1) then
        call cgyro_error('For zonal flow test, set n_radial=1.')
        return
     endif

     n = 0

  else if (n_toroidal == 1) then

     ! Single linear mode (assume n=1, compute rho)

     k_theta = q/rmin
     rho     = ky/k_theta
     length  = box_size/(s*k_theta)

     n = 1

     call cgyro_info('Single-mode linear analysis.')

  else

     ! Multiple modes (n=0,1,2,...,n_toroidal-1)

     k_theta = q/rmin
     rho     = ky/k_theta
     length  = box_size/(s*k_theta)

     ! Now define individual k_thetas

     n = i_group_1

     k_theta = n*k_theta

     call cgyro_info('Multiple toroidal harmonics.')

  endif
  !-------------------------------------------------------------

  !------------------------------------------------------------------------
  ! ExB shear
  !
  if (abs(gamma_e) > 1e-10) then
     call cgyro_info('Triggered ExB shear.') 
     omega_eb = k_theta*length*gamma_e/(2*pi)
  endif
  !------------------------------------------------------------------------

  !-------------------------------------------------------------
  ! Species-dependent quantities
  !
  do is=1,n_species

     ! thermal velocity
     vth(is) = sqrt(temp(is)/mass(is))

     ! collision frequency
     nu(is) = nu_ee_in *(1.0*z(is))**4 &
          * dens(is) / dens_ele &
          * sqrt(mass_ele/mass(is)) * (temp_ele/temp(is))**1.5

  enddo
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Fourier index mapping
  !
  allocate(indx_xi(n_xi))
  do ix=1,n_xi
     indx_xi(ix) = ix-1
  enddo
  allocate(px(n_radial))
  do ir=1,n_radial
     px(ir) = -n_radial/2 + (ir-1)
  enddo
  if (zf_test_flag == 1) px(1) = 1
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! General geometry -- accessible only from interface 
  ! via parameters geo_ny_in and geo_yin_in
  !
  geo_numeq_flag = 0
  geo_ny = 0
  allocate(geo_yin(8,0:geo_ny))
  geo_yin(:,:) = 0.0
  if (equilibrium_model == 3) then
     geo_numeq_flag = 1
     geo_ny = geo_ny_in  
     deallocate(geo_yin)
     allocate(geo_yin(8,0:geo_ny))
     do j=0,geo_ny
        geo_yin(:,j) = geo_yin_in(:,j)
     enddo
  endif
  !-------------------------------------------------------------

end subroutine cgyro_make_profiles
