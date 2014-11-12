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
     call cgyro_error('ERROR: (CGYRO) n_xi must be even')
     return
  endif

  if (n_species > 6) then
     call cgyro_error('ERROR: (CGYRO) max n_species is 6')
     return
  endif

  ! Field consistency checks
  if (n_field > 1) then
     if (abs(betae_unit) < epsilon(0.0)) then
        call cgyro_error('ERROR: (CGYRO) BETAE_UNIT must be nonzero for electromagnetic simulation.')
        return
     endif
     if (ae_flag == 1) then
        call cgyro_error('ERROR: (CGYRO) Electrons must be gyrokinetic for electromagnetic simulation.')
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
     call cgyro_error('ERROR: (CGYRO) invalid n_field')
     return
  end select

  if (collision_model == 0 .and. ae_flag == 1) then
     call cgyro_error('ERROR: (CGYRO) collision_model=0 requires kinetic electrons')
     return
  endif

  ! Collision model
  !
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

     call cgyro_error('ERROR: (CGYRO) invalid collision_model')
     return

  end select

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  

  case (0) 
     call cgyro_info('Equlibrium model = s-alpha')

  case (1) 
     call cgyro_error('ERROR: (CGYRO) equilibrium_model invalid')
     return

  case (2) 
     call cgyro_info('Equlibrium model = Miller')

  case (3) 
     call cgyro_info('Equlibrium model = General (Fourier)')

     if (geo_ny <= 0) then
        call cgyro_error('ERROR: (CGYRO) geometry coefficients missing')
        return
     endif

  case default

     call cgyro_error('ERROR: (CGYRO) equilibrium_model invalid')
     return

  end select
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Check local profile params

  do is=1,n_species
     if(dens(is) <= 0.0) then
        call cgyro_error('ERROR: (CGYRO) density must be positive')
        return
     end if
     if(temp(is) <= 0.0) then
        call cgyro_error('ERROR: (CGYRO) temperature must be positive')
        return
     end if
     if(nu(is) < 0.0) then
        call cgyro_error('ERROR: (CGYRO) collision frequency must be positive')
        return
     end if
     if(z(is) == 0.0) then
        call cgyro_error('ERROR: (CGYRO) charge must be non-zero')
        return
     end if
  enddo

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_run,file=trim(path)//runfile,status='old',position='append')
     write(io_run,*)
     write(io_run,*) 'GRID DIMENSIONS'
     write(io_run,*) '---------------'
     write(io_run,10) 'n_radial',n_radial
     write(io_run,10) 'n_energy',n_energy
     write(io_run,10) 'n_xi',n_xi
     write(io_run,10) 'n_theta',n_theta

     write(io_run,*) 
     write(io_run,*) 'PHYSICS PARAMETERS'
     write(io_run,*) '------------------'
     write(io_run,20) 'r/R',rmin
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

     do is=1,n_species
        write(io_run,*) 
        write(io_run,'(t2,a,i1)') 'Species ',is
        write(io_run,*) '----------'
        write(io_run,10) 'Z',z(is)
        write(io_run,20) 'dens',dens(is)
        write(io_run,20) 'temp',temp(is)
        write(io_run,20) 'mass',mass(is)
        write(io_run,20) 'a/Ln',dlnndr(is)
        write(io_run,20) 'a/LT',dlntdr(is)
        write(io_run,20) 'nu',nu(is)
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

10 format(t2,a,':',t14,i2)
20 format(t2,a,':',t13,1pe12.4) 

end subroutine cgyro_check
