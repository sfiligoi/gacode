subroutine cgyro_make_profiles

  use cgyro_globals

  implicit none

  integer :: is, num_ele, j
  integer, parameter :: io=20
  
  num_ele = 0
  do is=1,n_species
     if (z(is) == -1) then
        num_ele = num_ele + 1
     endif
  enddo
  if (num_ele == 0) then
     ae_flag = 1
  else if (num_ele == 1) then
     ae_flag = 0
  else
     call cgyro_error('ERROR: (CGYRO) Only one electron species allowed.')
     return
  endif
  
  do is=1,n_species
     nu(is) = nu_1_in *(1.0*z(is))**4/(1.0*z(1))**4 &
          * dens(is) / dens(1) &
          * sqrt(mass(1)/mass(is)) * (temp(1)/temp(is))**1.5
  enddo

  ! Standard local simulation (one point)
  
  q   = abs(q) 
  rho = abs(rho)
  
  if (zf_test_flag == 1) then

     ! Zonal flow (n=0) test
  
     toroidal_num = 0
     k_theta      = 0.0
     k_theta_rho  = ky
     r_length_inv = ky

  else if(toroidal_model == 0) then
     ! k_theta_rho and n are specified; compute rho
     k_theta_rho = abs(k_theta_rho)
     rho = k_theta_rho * rmin / (q * toroidal_num)
     k_theta = k_theta_rho / rho
     r_length_inv =  q * toroidal_num * shat / rmin
  else if(toroidal_model == 1) then
     ! rho and n are specified; compute k_theta
     rho  = abs(rho) 
     k_theta = (q * toroidal_num) / rmin
     k_theta_rho = (q * toroidal_num) / rmin * rho
     k_theta_rho = abs(k_theta_rho)
     r_length_inv =  q * toroidal_num * shat / rmin
  endif

  ! general geometry -- accessible only from interface 
  ! via parameters geo_ny_in and geo_yin_in
  geo_numeq_flag = 0
  geo_ny = 0
  allocate(geo_yin(8,0:geo_ny))
  geo_yin(:,:) = 0.0
  if(equilibrium_model == 3) then
     geo_numeq_flag = 1
     geo_ny = geo_ny_in  
     deallocate(geo_yin)
     allocate(geo_yin(8,0:geo_ny))
     do j=0,geo_ny
        geo_yin(:,j) = geo_yin_in(:,j)
     enddo
  endif
  
  do is=1, n_species
     vth(is) = sqrt(temp(is)/mass(is))
  enddo

  if(adiabatic_ele_model == 0) then
     do is=1, n_species
        if(Z(is) == -1) then
           is_ele = is
           exit
        endif
     enddo
  endif
  if(adiabatic_ele_model == 0) then
     dens_ele = dens(is_ele)
     temp_ele = temp(is_ele)
  else
     dens_ele = ne_ade
     temp_ele = te_ade
  endif

  ! Print the re-mapped equilibrium data
  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//'out.cgyro.equil',status='replace')
     write (io,'(e16.8)',advance='no') rmin
     write (io,'(e16.8)',advance='no') rmaj
     write (io,'(e16.8)',advance='no') q
     write (io,'(e16.8)',advance='no') shat
     write (io,'(e16.8)',advance='no') rho
     write (io,'(e16.8)',advance='no') k_theta_rho
     do is=1,n_species
        write (io,'(e16.8)',advance='no') dens(is)
        write (io,'(e16.8)',advance='no') temp(is)
        write (io,'(e16.8)',advance='no') dlnndr(is)
        write (io,'(e16.8)',advance='no') dlntdr(is)
        write (io,'(e16.8)',advance='no') nu(is)
     enddo
     write (io,*)
     close(io)
  end if

end subroutine cgyro_make_profiles
