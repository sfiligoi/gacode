subroutine read_profile_vugyro(dir)

  use gyro_globals

  implicit none

  real, dimension(:), allocatable :: dum

  character (len=*) :: dir

  open(unit=1,file=trim(dir)//'profile_vugyro.out')
  read(1,*) n_x
  read(1,*) n_theta_section
  read(1,*) n_pass
  read(1,*) n_trap
  read(1,*) n_energy
  read(1,*) n_theta_plot
  read(1,*) n0
  read(1,*) n_n
  read(1,*) d_n
  read(1,*) n_explicit_damp
  read(1,*) nonlinear_flag
  read(1,*) electron_method
  read(1,*) n_field
  read(1,*) n_ion
  read(1,*) n_kinetic
  read(1,*) n_spec
  read(1,*) field_r0_flag
  read(1,*) field_r0_grid
  read(1,*) n_grid_exp
  read(1,*) boundary_method

  !-----------------------------------------
  ! Set some dimensions for unused variables 
  n_n_1  = 1
  m_gyro = 1
  m_dx   = 1
  i_gyro = 0
  source_flag = 0
  !----------------------------------------

  call gyro_alloc_profile_sim(1)

  allocate(dum(n_x))

  read(1,*) r(:)
  read(1,*) q(:)
  read(1,*) r_s(:)
  read(1,*) q_s(:)
  read(1,*) dlntdr_s(:,:)
  read(1,*) dlnndr_s(:,:)
  read(1,*) tem_s(:,:)
  read(1,*) den_s(:,:)
  read(1,*) dum(:)       ! Need to add new variable
  read(1,*) delta_s(:)
  read(1,*) zeta_s(:)
  read(1,*) kappa_s(:)
  read(1,*) drmaj_s(:)
  read(1,*) shat_s(:)
  read(1,*) s_delta_s(:)
  read(1,*) s_zeta_s(:)
  read(1,*) s_kappa_s(:)
  read(1,*) zmag_s(:)
  read(1,*) dzmag_s(:)
  read(1,*) beta_unit_s(:)
  read(1,*) gamma_e_s(:) 
  read(1,*) gamma_p_s(:) 
  read(1,*) mach_s(:) 
  read(1,*) b_unit_s(:)
  read(1,*) dr_eodr(:)
  read(1,*) z_eff_s(:)
  read(1,*) nu_s(n_spec,:)
  read(1,*) w0_s(:)

  read(1,*) box_multiplier

  close(1)

end subroutine read_profile_vugyro
