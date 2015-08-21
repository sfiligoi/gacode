subroutine cgyro_read_input

  use cgyro_globals

  implicit none

  integer :: is
  character (len=1) :: cdummy

  open(unit=1,file=trim(path)//'input.cgyro.gen',status='old')
  read(1,*) n_energy
  read(1,*) n_xi
  read(1,*) n_theta
  read(1,*) n_radial
  read(1,*) n_toroidal
  read(1,*) n_field
  read(1,*) e_max
  read(1,*) delta_t
  read(1,*) max_time
  read(1,*) print_step
  read(1,*) restart_step
  read(1,*) freq_tol
  read(1,*) restart_write
  read(1,*) restart_mode
  read(1,*) up_radial
  read(1,*) up_theta
  read(1,*) nup
  read(1,*) implicit_flag
  read(1,*) constant_wind_flag
  read(1,*) ky
  read(1,*) box_size
  read(1,*) silent_flag
  read(1,*) profile_model
  read(1,*) equilibrium_model
  read(1,*) collision_model
  read(1,*) collision_mom_restore
  read(1,*) collision_ene_restore
  read(1,*) collision_ene_diffusion
  read(1,*) collision_kperp
  read(1,*) collision_field_model
  read(1,*) collision_trap_model
  read(1,*) zf_test_flag
  read(1,*) nonlinear_flag
  read(1,*) nonlinear_method
  read(1,*) te_ade
  read(1,*) ne_ade
  read(1,*) masse_ade
  read(1,*) lambda_debye
  read(1,*) lambda_debye_scale
  read(1,*) test_flag
  read(1,*) h_print_flag
  read(1,*) amp
  read(1,*) gamma_e
  read(1,*) gamma_p
  read(1,*) mach
  read(1,*) gamma_e_scale
  read(1,*) gamma_p_scale
  read(1,*) mach_scale

  read(1,*) rmin
  read(1,*) rmaj
  read(1,*) q
  read(1,*) s
  read(1,*) shift    
  read(1,*) kappa   
  read(1,*) s_kappa  
  read(1,*) delta       
  read(1,*) s_delta
  read(1,*) zeta      
  read(1,*) s_zeta
  read(1,*) zmag       
  read(1,*) s_zmag
  read(1,*) beta_star
  read(1,*) betae_unit

  read(1,*) subroutine_flag

  read(1,*) n_species

  read(1,*) nu_ee

  do is=1,6
     read(1,*) z(is)
     read(1,*) mass(is)
     read(1,*) dens(is)
     read(1,*) temp(is)
     read(1,*) dlnndr(is)
     read(1,*) dlntdr(is)
  enddo

  close(1)

  ! GEO fourier coefficients
  geo_ny_in = 0
  geo_yin_in(:,:) = 0.0
  if (subroutine_flag == 0 .and. equilibrium_model == 3 & 
       .and. profile_model == 1) then
     open(unit=1,file=trim(path)//'input.geo',status='old')
     ! header skip
     do
        read(1,'(a)') cdummy
        if (cdummy /= '#') exit
     enddo
     backspace 1
     ! n_fourier
     read(1,*) geo_ny_in
     ! fourier coefficients
     read(1,*) geo_yin_in(:,0:geo_ny_in)
     close(1)
  endif

end subroutine cgyro_read_input
