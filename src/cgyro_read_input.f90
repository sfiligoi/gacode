subroutine cgyro_read_input

  use cgyro_globals

  implicit none

  integer :: is

  open(unit=1,file=trim(path)//'input.cgyro.gen',status='old')
  read(1,*) n_energy
  read(1,*) n_xi
  read(1,*) n_theta
  read(1,*) n_radial
  read(1,*) e_max
  read(1,*) delta_t
  read(1,*) max_time
  read(1,*) freq_tol
  read(1,*) restart_write
  read(1,*) restart_mode
  read(1,*) rupwind_eps
  read(1,*) rupwind_n
  read(1,*) tupwind_eps
  read(1,*) toroidal_model
  read(1,*) toroidal_num
  read(1,*) rho
  read(1,*) k_theta_rho
  read(1,*) r_length_rho
  read(1,*) rmin
  read(1,*) rmaj
  read(1,*) silent_flag
  read(1,*) equilibrium_model
  read(1,*) collision_model
  read(1,*) profile_model
  read(1,*) ipccw
  read(1,*) btccw
  read(1,*) te_ade
  read(1,*) ne_ade
  read(1,*) lambda_debye
  read(1,*) profile_lambda_debye_scale

  read(1,*) n_species

  do is=1,6
     read(1,*) z(is)
     read(1,*) mass(is)
     read(1,*) dens(is)
     read(1,*) temp(is)
     read(1,*) dlnndr(is)
     read(1,*) dlntdr(is)
     read(1,*) nu(is)
  enddo

  read(1,*) q
  read(1,*) shat
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

  close(1)

  ! GEO fourier coefficients are not yet available to read-in
  geo_ny_in = 0
  geo_yin_in(:,:) = 0.0

end subroutine cgyro_read_input
