subroutine gkcoll_make_profiles

  use gkcoll_globals
  use gkcoll_profile_exp, only : &
       rmaj_p,&
       b_unit_p,& 
       q_exp,&
       n_grid_exp
  use gkcoll_allocate_profile

  implicit none

  integer :: ir, is, ip, num_ele, j
  integer, parameter :: io=20

  real :: cc
  real :: loglam

  call PROFILE_SIM_alloc(1)

  ! radially-local problem only
  r(1) = rmin_in

  num_ele = 0
  do is=1, n_species
     if(z_in(is) == -1) then
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
  
  if(btccw_in > 0) then
     sign_bunit = -1.0
  else
     sign_bunit =  1.0
  endif

  if(ipccw_in > 0) then
     sign_q = -sign_bunit
  else
     sign_q =  sign_bunit
  endif

  select case (profile_model) 

  case (1)

     ! Standard local simulation (one point)

     ir = 1

     rmaj(ir)      = rmaj_in
     q(ir)         = abs(q_in) * sign_q
     rho(ir)       = abs(rho_in) * sign_bunit
     shat(ir)      = shat_in      
     shift(ir)     = shift_in     
     kappa(ir)     = kappa_in     
     s_kappa(ir)   = s_kappa_in   
     delta(ir)     = delta_in    
     s_delta(ir)   = s_delta_in
     zeta(ir)      = zeta_in    
     s_zeta(ir)    = s_zeta_in
     zmag(ir)      = zmag_in    
     s_zmag(ir)    = s_zmag_in

     ! general geometry -- accessible only from interface 
     ! via parameters geo_ny_in and geo_yin_in
     if(equilibrium_model == 3) then
        geo_numeq_flag = 1
        deallocate(geo_yin)
        geo_ny = geo_ny_in   
        allocate(geo_yin(8,0:geo_ny,1))
        do j=0,geo_ny
           geo_yin(:,j,1) = geo_yin_in(:,j)
        enddo
     endif

     te_ade(ir)    = te_ade_in
     ne_ade(ir)    = ne_ade_in

     do is=1,n_species
        z(is)         = z_in(is)
        mass(is)      = mass_in(is)
        dens(is,ir)   = dens_in(is)
        temp(is,ir)   = temp_in(is)
        dlnndr(is,ir) = dlnndr_in(is)
        dlntdr(is,ir) = dlntdr_in(is)
        nu(is,ir)     = nu_in(is)
     enddo

     do is=1, n_species
        vth(is,ir) = sqrt(temp(is,ir)/mass(is))
     enddo

     ! These normalizations are arbitrary for local profiles
     temp_norm_fac = 1.0
     charge_norm_fac = 1.0
     dens_norm(:) = 1.0
     vth_norm(:) = 1.0
     b_unit(:) = 1.0
     a_meters        = 1.0
     temp_norm(:)    = 1.0


  case (2)

     ! Standard simulation with experimental profiles

     do is=1,n_species
        z(is)    = z_in(is)
        mass(is) = mass_in(is)
     enddo

     call gkcoll_experimental_profiles
     if(error_status > 0) return
     call gkcoll_map_experimental_profiles
           
     ! Normalizing quantities
     temp_norm_fac   = 1.6022*1000
     charge_norm_fac = 1.6022
     dens_norm(:) = dens(1,:)
     temp_norm(:) = temp(1,:)

     ! Compute vth/a (1/s) using dimensional quantities.  
     ! mass(i) is thus measured in units of deuterium mass.
     do is=1,n_species
        vth(is,:) = sqrt(temp(is,:) * temp_norm_fac &
             / (mass(is) * mass_deuterium)) &
             * 1.0e4 / a_meters
     enddo
     vth_norm(:)  = vth(1,:) * sqrt(mass(1))

     ! Determine the equilibrium parameters
     select case (profile_equilibrium_model) 
     case (0)
        ! input equilibrium not used -- re-set for s-alpha geometry
        ! use the input profile-averaged rmaj, q, and b_unit
        equilibrium_model = 0
        do ir=1, n_gr
           rmaj(ir)   = 0.0
           q(ir)      = 0.0
           b_unit(ir) = 0.0
           do ip=1,n_grid_exp
              rmaj(ir)   = rmaj(ir)   + rmaj_p(ip)
              q(ir)      = q(ir)      + q_exp(ir)
              b_unit(ir) = b_unit(ir) + b_unit_p(ip)
           enddo
           rmaj(ir)   = rmaj(ir)   / (1.0 * n_grid_exp)
           q(ir)      = q(ir)      / (1.0 * n_grid_exp)
           b_unit(ir) = b_unit(ir) / (1.0 * n_grid_exp)
        enddo

     case(1)
        ! use the input equilibrium parameters (miller)
        equilibrium_model = 2

     case(2)
        ! use the input equilibrium parameters (general)
        equilibrium_model = 3
 
    end select

     ! Compute rho/a for species 1 using dimensional quantities
     ! mass(i) is thus measured in units of deuterium mass.
     rho(:) = sqrt(temp(1,:) * temp_norm_fac &
             * mass(1) * mass_deuterium) &
             / (charge_norm_fac * b_unit(:)) &
             * 1.0e-4 / a_meters

     ! Compute collision frequency

     ! Numerical coefficient (relative to M_D)
     cc = sqrt(2.0) * pi * charge_norm_fac**4 &
          * 1.0 / (4.0 * pi * 8.8542)**2 &
          * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
          * 1e9

     do is=1,n_species
        do ir=1, n_gr

           loglam = 24.0 - log(sqrt(ne_ade(ir)*1e13)/(te_ade(ir)*1000))

           ! Collision rate (1/sec)
           nu(is,ir) = cc * loglam * dens(is,ir) * z(is)**4 &
                / (sqrt(mass(is)) * temp(is,ir)**1.5)

           ! Express in local dimensionless GKCOLL units:
           nu(is,ir) = nu(is,ir)/vth_norm(ir)

        enddo
     enddo

     do is=1,n_species
        do ir=1, n_gr
           ! Normalize N and T to value at r for species 1.
           dens(is,ir) = dens(is,ir)/dens_norm(ir)
           temp(is,ir) = temp(is,ir)/temp_norm(ir)
           if(is == 1) then
              te_ade(ir) = te_ade(ir) / temp_norm(ir)
              ne_ade(ir) = ne_ade(ir) / dens_norm(ir)
           endif
           ! Normalize vth/a
           vth(is,ir) = vth(is,ir)/vth_norm(ir)
        enddo
     enddo

     call PROFILE_EXP_alloc(0)

  end select

  ! Print the re-mapped equilibrium data
  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//'out.gkcoll.equil',status='replace')
     do ir=1,n_gr
        write (io,'(e16.8,$)') r(ir)
        write (io,'(e16.8,$)') q(ir)
        write (io,'(e16.8,$)') rho(ir)
        write (io,'(e16.8,$)') rmaj(ir)
        write (io,'(e16.8,$)') dens_norm(ir)
        write (io,'(e16.8,$)') temp_norm(ir)
        write (io,'(e16.8,$)') vth_norm(ir)
        do is=1,n_species
           write (io,'(e16.8,$)') dens(is,ir)
           write (io,'(e16.8,$)') temp(is,ir)
           write (io,'(e16.8,$)') dlnndr(is,ir)
           write (io,'(e16.8,$)') dlntdr(is,ir)
           write (io,'(e16.8,$)') nu(is,ir)
        enddo
        write (io,*)
     enddo
     close(io)
  end if

end subroutine gkcoll_make_profiles
