subroutine gkcoll_make_profiles

  use gkcoll_globals
  use gkcoll_profile_exp, only : &
       rmaj_p,&
       b_unit_p,& 
       q_exp,&
       n_grid_exp
  use gkcoll_allocate_profile

  implicit none

  integer :: is, num_ele, j
  integer, parameter :: io=20

  real :: cc
  real :: loglam

  call PROFILE_SIM_alloc(1)

  num_ele = 0
  do is=1, n_species
     if(z(is) == -1) then
        num_ele = num_ele + 1
     endif
  enddo
  if(num_ele == 0) then
     adiabatic_ele_model = 1
  else if(num_ele == 1) then
     adiabatic_ele_model = 0
  else
     call gkcoll_error('ERROR: (GKCOLL) Only one electron species allowed')
     return
  end if
  
  if(btccw > 0) then
     sign_bunit = -1.0
  else
     sign_bunit =  1.0
  endif

  if(ipccw > 0) then
     sign_q = -sign_bunit
  else
     sign_q =  sign_bunit
  endif

  select case (profile_model) 

  case (1)

     ! Standard local simulation (one point)

     q    = abs(q) * sign_q

     if(toroidal_model == 2) then
        ! n=0 test
        ! specify rho and r_length_inv
        toroidal_num = 0
        k_theta      = 0.0
        k_theta_rho  = 0.0
        rho  = abs(rho) * sign_bunit
        r_length_inv = 1.0 / (r_length_rho) / rho
     else if(toroidal_model == 0) then
        ! k_theta_rho and n are specified; compute rho
        k_theta_rho = abs(k_theta_rho)
        rho = k_theta_rho * rmin / (q * toroidal_num)
        k_theta = k_theta_rho / rho
        r_length_inv =  q * toroidal_num * shat / rmin
     else if(toroidal_model == 1) then
        ! rho and n are specified; compute k_theta
        rho  = abs(rho) * sign_bunit
        k_theta = (q * toroidal_num) / rmin
        k_theta_rho = (q * toroidal_num) / rmin * rho
        k_theta_rho = abs(k_theta_rho)
        r_length_inv =  q * toroidal_num * shat / rmin
     endif



     ! general geometry -- accessible only from interface 
     ! via parameters geo_ny_in and geo_yin_in
     if(equilibrium_model == 3) then
        geo_numeq_flag = 1
        deallocate(geo_yin)
        geo_ny = geo_ny_in   
        allocate(geo_yin(8,0:geo_ny))
        do j=0,geo_ny
           geo_yin(:,j) = geo_yin_in(:,j)
        enddo
     endif

     do is=1, n_species
        vth(is) = sqrt(temp(is)/mass(is))
     enddo

     ! These normalizations are arbitrary for local profiles
     temp_norm_fac   = 1.0
     charge_norm_fac = 1.0
     a_norm          = 1.0
     dens_norm       = 1.0
     temp_norm       = 1.0
     vth_norm        = 1.0
     b_norm          = 1.0

  case (2)

     ! Standard simulation with experimental profiles

     call gkcoll_experimental_profiles
     if(error_status > 0) return
     call gkcoll_map_experimental_profiles
           
     ! Normalizing quantities
     ! a_norm and b_norm are set in map_experimental_profiles
     temp_norm_fac   = 1.6022*1000
     charge_norm_fac = 1.6022
     dens_norm       = dens(1)
     temp_norm       = temp(1)

     ! Compute vth/a (1/s) using dimensional quantities.  
     ! mass(i) is thus measured in units of deuterium mass.
     do is=1,n_species
        vth(is) = sqrt(temp(is) * temp_norm_fac &
             / (mass(is) * mass_deuterium)) &
             * 1.0e4 / a_norm
     enddo
     vth_norm  = vth(1) * sqrt(mass(1))

     ! Compute rho/a for species 1 using dimensional quantities
     ! mass(i) is thus measured in units of deuterium mass.
     rho = sqrt(temp(1) * temp_norm_fac &
             * mass(1) * mass_deuterium) &
             / (charge_norm_fac * b_norm) &
             * 1.0e-4 / a_norm

     ! rho and n are specified; compute k_theta
     k_theta = (q * toroidal_num) / rmin
     k_theta_rho = (q * toroidal_num) / rmin * rho
     k_theta_rho = abs(k_theta_rho)
     ! compute r_length
      r_length_inv =  q * toroidal_num * shat / rmin

     ! Debye length
     lambda_debye = profile_lambda_debye_scale*7.43* &
          sqrt((1e3*te_ade)/(1e13*ne_ade))/a_norm

     ! Compute collision frequency

     ! Numerical coefficient (relative to M_D)
     cc = sqrt(2.0) * pi * charge_norm_fac**4 &
          * 1.0 / (4.0 * pi * 8.8542)**2 &
          * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
          * 1e9

     do is=1,n_species

        loglam = 24.0 - log(sqrt(ne_ade*1e13)/(te_ade*1000))

        ! Collision rate (1/sec)
        nu(is) = cc * loglam * dens(is) * z(is)**4 &
             / (sqrt(mass(is)) * temp(is)**1.5)
        
        ! Express in local dimensionless GKCOLL units:
        nu(is) = nu(is)/vth_norm

     enddo

     do is=1,n_species
        ! Normalize n and t to value at r for species 1.
        dens(is) = dens(is)/dens_norm
        temp(is) = temp(is)/temp_norm
        if(is == 1) then
           te_ade = te_ade / temp_norm
           ne_ade = ne_ade / dens_norm
        endif
        ! Normalize vth/a
        vth(is) = vth(is)/vth_norm
     enddo
     
     call PROFILE_EXP_alloc(0)
     
  end select

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
     open(unit=io,file=trim(path)//'out.gkcoll.equil',status='replace')
     write (io,'(e16.8,$)') rmin
     write (io,'(e16.8,$)') rmaj
     write (io,'(e16.8,$)') q
     write (io,'(e16.8,$)') shat
     write (io,'(e16.8,$)') rho
     write (io,'(e16.8,$)') k_theta_rho
     write (io,'(e16.8,$)') dens_norm
     write (io,'(e16.8,$)') temp_norm
     write (io,'(e16.8,$)') vth_norm
     write (io,'(e16.8,$)') a_norm
     write (io,'(e16.8,$)') b_norm
     do is=1,n_species
        write (io,'(e16.8,$)') dens(is)
        write (io,'(e16.8,$)') temp(is)
        write (io,'(e16.8,$)') dlnndr(is)
        write (io,'(e16.8,$)') dlntdr(is)
        write (io,'(e16.8,$)') nu(is)
     enddo
     write (io,*)
     close(io)
  end if

end subroutine gkcoll_make_profiles
