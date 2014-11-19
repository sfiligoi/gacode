subroutine cgyro_check

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is
  logical :: lfe

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_run,file=trim(path)//runfile,&
          status='old',position='append')
  endif

  !-----------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_xi,2) /= 0) then 
     call cgyro_error('n_xi must be even.')
     return
  endif

  if (n_species > 6) then
     call cgyro_error('n_species <= 6.')
     return
  endif

  ! Field consistency checks
  if (n_field > 1) then
     if (abs(betae_unit) < epsilon(0.0)) then
        call cgyro_error('BETAE_UNIT must be nonzero for electromagnetic simulation.')
        return
     endif
     if (ae_flag == 1) then
        call cgyro_error('Electrons must be gyrokinetic for electromagnetic simulation.')
        return
     endif
  endif

  ! Fields
  select case(n_field)

  case (1)
     call cgyro_info('Electrostatic fluctuations (Phi)')

  case (2)
     call cgyro_info('Transverse electromagnetic fluctuations (Phi,A_par)')

  case (3)
     call cgyro_info('Compressional electromagnetic fluctuations (Phi,A_par,B_par)')
     stop
  case default
     call cgyro_error('Invalid value for n_field.')
     return
  end select

  if (collision_model == 1 .and. ae_flag == 1) then
     call cgyro_error('collision_model=1 requires kinetic electrons')
     return
  endif

  ! Collision model

  select case (collision_model)  

  case(0)
     call cgyro_info('Collision model = none')

  case (1) 
     call cgyro_info('Collision model = CONNOR EE+EI LORENTZ only')

  case (2) 
     call cgyro_info('Collision model = Connor')

  case (3) 
     call cgyro_info('Collision model = Reduced Hirshman-Sigmar')

  case (4) 
     call cgyro_info('Collision model = Ad hoc Fokker-Planck')

  case default

     call cgyro_error('Invalid value for collision_model')
     return

  end select

  select case (collision_mom_restore)
     case(0)
        call cgyro_info('Collision momentum restoring = not included')
     case(1)
        call cgyro_info('Collision momentum restoring = included')
     case default
        call cgyro_error('Invalid value for collision_mom_restore')
        return
     end select

  if(collision_model == 4) then
     select case (collision_ene_diffusion)
     case(0)
        call cgyro_info('Collision energy diffusion = not included')
     case(1)
        call cgyro_info('Collision energy diffusion = included')
     case default
        call cgyro_error('Invalid value for collision_ene_diffusion')
        return
     end select
     select case (collision_ene_restore)
     case(0)
        call cgyro_info('Collision energy restoring = not included')
     case(1)
        call cgyro_info('Collision energy restoring = included')
     case default
        call cgyro_error('Invalid value for collision_ene_restore')
        return
     end select
  endif


  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  

  case (0) 
     call cgyro_info('Equlibrium model = s-alpha')

  case (1) 
     call cgyro_error('Invalid value for equilibrium_model.')
     return

  case (2) 
     call cgyro_info('Equlibrium model = Miller')

  case (3) 
     call cgyro_info('Equlibrium model = General (Fourier)')

     if (geo_ny <= 0) then
        call cgyro_error('Fourier geometry coefficients missing.')
        return
     endif

  case default

     call cgyro_error('ERROR: (CGYRO) equilibrium_model invalid')
     return

  end select
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Check profile parameters
  !
  do is=1,n_species
     if (dens(is) <= 0.0) then
        call cgyro_error('Densities must be positive.')
        return
     endif
     if (temp(is) <= 0.0) then
        call cgyro_error('Temperatures must be positive.')
        return
     endif
     if (nu(is) < 0.0) then
        call cgyro_error('Collision frequencies must be non-negative.')
        return
     endif
     if (z(is) == 0) then
        call cgyro_error('Charge must be non-zero.')
        return
     endif
  enddo

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_run,file=trim(path)//runfile,status='old',position='append')
     write(io_run,*)
     write(io_run,'(a)') 'n_radial  n_theta   n_species n_energy  n_xi'
     write(io_run,'(t1,5(i4,6x))') n_radial,n_theta,n_species,n_energy,n_xi

     write(io_run,*) 
     write(io_run,'(a,f6.2,a,f6.2)') '(Lx,Ly)/rho',1/r_length_inv/rho,',',2*pi/(q/rmin)
     write(io_run,*) 
     if (n_toroidal == 1) then
        write(io_run,20) 'ky*rho',ky
     else
        write(io_run,20) 'min(ky*rho)',(q/rmin)*rho
        write(io_run,20) 'max(ky*rho)',(n_toroidal-1)*(q/rmin)*rho
     endif
     write(io_run,*) 
     write(io_run,20) 'min(kx*rho)',2*pi*r_length_inv*rho
     write(io_run,20) 'max(kx*rho)',2*pi*r_length_inv*rho*(n_radial/2-1)
     write(io_run,*) 
     write(io_run,20) 'r/a',rmin
     write(io_run,20) 'R/a',rmaj
     write(io_run,20) 'q',q
     write(io_run,20) 's',s
     write(io_run,20) 'shift',shift
     write(io_run,20) 'kappa',kappa
     write(io_run,20) 's_kappa',s_kappa
     write(io_run,20) 'delta',delta
     write(io_run,20) 's_delta',s_delta
     write(io_run,20) 'zeta',zeta
     write(io_run,20) 's_zeta',s_zeta
     write(io_run,20) 'zmag',zmag
     write(io_run,20) 's_zmag',s_zmag

     write(io_run,*)
     write(io_run,'(a)') 'indx  z    n/n_norm    T/T_norm    m/m_norm     a/Ln        a/Lt        nu'
     do is=1,n_species
        write(io_run,'(t2,i2,2x,i2,2x,6(1pe10.4,2x))') is,z(is),dens(is),temp(is),mass(is),dlnndr(is),dlntdr(is),nu(is)
     enddo

     write(io_run,*)

     close(io_run)
  endif

  ! Check for existence of restart file
  if (restart_mode == 1) then

     inquire(file=trim(path)//runfile_restart,exist=lfe)
     if (lfe .eqv. .false.) then
        call cgyro_error('ERROR: (CGYRO) Missing restart file')
        return
     endif

  endif

20 format(t2,a,':',t13,1pe12.4) 

end subroutine cgyro_check
