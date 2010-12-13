!---------------------------------------------------------------
! gyro_select_methods.f90 [caller BigScience]
!
! PURPOSE:
!
!  Run through all method selection parameters.  As much 
!  consistency checking as possible should be done here.
!---------------------------------------------------------------

subroutine gyro_select_methods

  use gyro_globals

  !--------------------------
  implicit none
  character (len=40) :: mstr
  !--------------------------

  !-------------------------------------------------------------
  ! Various consistency and error checks
  ! 
  ! Special diffusivity check:
  !
  ! Grid offset checks:
  !
  if (n_x_offset /= 0 .and. boundary_method == 1) then
     call catch_error('ERROR: RADIAL_GRID_OFFSET must be zero.')
  endif
  if (abs(n_x_offset) > n_x/2-1) then
     call catch_error('ERROR: RADIAL_GRID_OFFSET too large.')
  endif
  !
  ! Study mode check:
  !
  if (n_study < 0 .or. n_study > n_n) then
     call catch_error('ERROR: N_STUDY out of bounds.')
  endif
  !
  ! Toroidal grid check:
  !
  if (nonlinear_flag == 1 .and. n_n < 2) then
     call catch_error('ERROR: Nonlinear runs require TOROIDAL_GRID > 1.')
  endif
  !
  ! Distribution function print check: 
  !
  if (dist_print == 1) then
     if  (n_n > 1) then
        call catch_error('ERROR: DIST_PRINT_FLAG=1 requires TOROIDAL_GRID=1.')
     endif
     if  (n_proc > 1) then
        call catch_error('ERROR: DIST_PRINT_FLAG=1 on one processor only.')
     endif
  endif
  !
  ! Radial grid check:
  !
  if (n_x < 4) then
     call catch_error('ERROR: You must have RADIAL_GRID >= 4.')
  endif
  !
  ! Radial gyroaverage check:
  !
  if (m_gyro > n_x/2) then
     call catch_error('ERROR: RADIAL_GYRO_BAND cannot exceed half RADIAL_GRID. ')
  endif
  !
  ! Buffer check:
  !
  if (boundary_method == 2 .and. n_x-2*n_explicit_damp <= 0) then
     call catch_error('ERROR: Buffers are too large (see EXPLICIT_DAMP_GRID).')
  endif
  !
  ! Source check:
  !
  if (n_source <= 0) then
     call catch_error('ERROR: require N_SOURCE > 0.')
  endif
  !
  ! Check for validity of transfer diagnostic
  !
  if (boundary_method == 2 .and. nonlinear_transfer_flag == 1) then
     call catch_error('ERROR: Need BOUNDARY_METHOD=1 for NONLINEAR_TRANSFER_FLAG=1.')
  endif
  !
  ! Check for solver consistency
  !
  if (linsolve_method > 1 .and. nonlinear_flag == 1) then
     call catch_error('ERROR: Can only have NONLINEAR_FLAG=1 with LINSOLVE_METHOD=1.')
  endif
  !
  if (radial_upwind < 0.0) then
     call catch_error('ERROR: RADIAL_UPWIND must be > 0.')
  endif
  !---------------------------------------------------

  !---------------------------------------------------
  ! Total number of pitch-angles
  !
  n_lambda = n_pass+n_trap
  !---------------------------------------------------

  !---------------------------------------------------
  ! n_moment: number of traditional moments (n,E)
  ! p_moment: number of new moments (n,E,parallel,k,heating)

  n_moment = 2
  p_moment = 4
  !---------------------------------------------------

  !---------------------------------------------------
  ! Figure out species numbers:
  !
  n_spec = 1
  do while (n_vec(n_spec) > 0.0) 
     n_spec = n_spec+1
  enddo
  n_ion = n_spec-1
  if (electron_method == 4) n_ion = n_ion+1
  !---------------------------------------------------

  !---------------------------------------------------
  ! Sort out species indices.  There may seem to be an
  ! awful lot of dimensions here, but all are needed 
  ! to support the 4 different ways of treating the 
  ! species (electron_method):
  !
  ! n_spec    : total number of ions and electrons
  ! n_ion     : total number of ions and gk electrons
  ! n_gk      : number of gyrokinetic species
  ! n_kinetic : number of kinetic species
  ! indx_e    : the electron index
  !
  select case (electron_method) 

  case (1)

     ! Adiabatic electrons

     n_kinetic = n_ion
     indx_e    = n_spec
     n_gk      = n_ion

  case (2)

     ! Drift-kinetic electrons

     n_kinetic = n_spec
     indx_e    = n_spec
     n_gk      = n_ion

  case (3)

     ! Adiabatic ions

     n_kinetic = 1
     indx_e    = 1
     n_gk      = 1

  case (4)

     ! Gyrokinetic electrons

     n_kinetic = n_ion
     indx_e    = n_ion
     n_gk      = n_ion

  end select
  !---------------------------------------------------

  !---------------------------------------------------
  ! Manage orbit upwind in special cases
  !
  select case (electron_method)
  case (2)
     ! No orbit upwind for implicit electrons
     orbit_upwind_vec(0) = 0.0  
  case (3)
     ! Permute i-e for adiabatic ions
     orbit_upwind_vec(0:1) = orbit_upwind_vec(1:0:-1) 
  case (4)
     ! Map electrons to final ion
     orbit_upwind_vec(n_ion) = orbit_upwind_vec(0)
  end select

  if (linsolve_method == 2) then
     orbit_upwind_vec(:) = 0.0
  endif
  !---------------------------------------------------

  !---------------------------------------------------
  ! Check on number of fields:
  !
  if (electron_method == 1) n_field = 1
  if (ampere_scale == 0.0) n_field = 1
  if (radial_profile_method /= 3 .and. betae_unit == 0.0) n_field = 1
  !
  if (n_field == 3) then
     call send_message(&
          'INFO: Rotation effects with delta B_parallel are experimental.')
  endif
  !---------------------------------------------------

  !---------------------------------------------------
  ! Set eparallel_plot_flag
  ! 
  eparallel_plot_flag = 1
  !---------------------------------------------------

  !---------------------------------------------------
  ! Determine if collisions will be included
  !
  if (electron_method == 1) then
     ! No electron equation, no electron collisions
     nu_ei_scale = 0.0
  endif
  if (electron_method == 3) then
     ! No ion equation, no ion collisions
     nu_ii_scale = 0.0
  endif
  !
  if (radial_profile_method == 3) then
     if (nu_ei_scale > 0.0 .or. nu_ii_scale > 0.0) then
        collision_flag = 1
     else
        collision_flag = 0
     endif
  else
     if (nu_ei*nu_ei_scale > 0.0 .or. nu_ei*nu_ii_scale > 0.0) then
        collision_flag = 1
     else
        collision_flag = 0
     endif
  endif
  !
  ! Catch some collision-related errors
  !
  if (collision_flag == 1) then

     if (n_field == 3 .and. collision_method /= 2) then
        call catch_error(&
             'ERROR: This collision model incompatible with delta B_parallel')
     endif

     if (n_pass < 4) then
        call catch_error(&
             'ERROR: N_PASS >= 4 required for collisions.')
     endif

     if (variable_egrid_flag == 1) then
        call catch_error(&
             'ERROR: Collisions not supported with variable energy grid')
     endif

  endif
  !---------------------------------------------------

  call send_line(separator)

  !----------------------------------------------------
  ! BASIC OPERATIONAL MODE
  !
  select case (linsolve_method)

  case (1)

     call send_line('operational mode     : INITIAL VALUE')
     fakefield_flag = 1

  case (2)

     call send_line('operational mode     : FULL GK EIGENVALUE')
     fakefield_flag = 1

  case (3)

     call send_line('operational mode     : FIELD EIGENVALUE')
     fakefield_flag = 0

  end select
  !----------------------------------------------------

  !----------------------------------------------------
  ! BOUNDARY CONDITIONS:
  !
  select case (boundary_method)

  case (1)

     call send_line('boundary conditions  : PERIODIC')

     ! Set damping region size to zero
     n_explicit_damp = 0

     ! We must have flat profile in this case:
     flat_profile_flag = 1

     ! Don't want a source
     source_flag = 0

     if (m_gyro == n_x/2) then
        ! Fully pseudospectral gyroaverages
        i_gyro = 1
     else
        ! Truncated gyroaverages
        i_gyro = 0
     endif

     if (m_dx == n_x/2) then
        ! Full complex form
        i_dx = 1
     else
        i_dx = 0
     endif

     if (radial_profile_method == 2) then

        ! profile_method == 2 only consistent with nonperiodic 
        ! mode of operation.  Avoid resets.

        call catch_error('INVALID: boundary_method/profile_method')

     endif

     if (nonuniform_grid_flag == 1) then
        call catch_error('INVALID: cannot have nonuniform grid.')
     endif

  case (2)

     call send_line('boundary conditions  : NONPERIODIC')

     ! No harmonic operators for nonperiodic case

     i_dx   = 0
     i_gyro = 0

  case default

     call catch_error('INVALID: boundary_method')

  end select
  !----------------------------------------------------

  !----------------------------------------------------
  ! PROFILE TYPE:
  !
  select case (radial_profile_method) 

  case (1)

     call send_line('profiles             : FLAT')

     geometry_method = -1
     flat_profile_flag = 1

  case (3)

     geometry_method = num_equil_flag

     if (flat_profile_flag == 0) then
        call send_line(&
             'profiles             : EXPERIMENTAL')
     else
        call send_line(&
             'profiles             : EXPERIMENTAL (flattened)')
     endif

  case (5)

     geometry_method = num_equil_flag

     call send_line('profiles             : FLAT')

     flat_profile_flag = 1

  case default

     call catch_error('INVALID: profile_method')

  end select
  !----------------------------------------------------

  !----------------------------------------------------
  ! GEOMETRY:
  !
  select case (geometry_method)

  case (-1) 
     call send_line('geometry             : S-ALPHA')

  case (0)
     call send_line('geometry             : MODEL SHAPE')

  case (1)
     call send_line('geometry             : GENERAL SHAPE (from EFIT)')

  case default
     call catch_error('INVALID: num_equil_flag')

  end select
  !----------------------------------------------------

  !----------------------------------------------------
  ! ELECTRON TREATMENT:
  !
  select case (electron_method) 

  case (1)

     call send_line('electron_method      : ADIABATIC ELECTRONS')

  case (2)

     call send_line('electron_method      : IMPLICIT DRIFT-KINETIC ELECTRONS')

  case (3)

     call send_line('electron_method      : ADIABATIC IONS')

  case (4)

     call send_line('electron_method      : ALL SPECIES GYROKINETIC')

  case default

     call catch_error('INVALID: electron_method')

  end select
  !----------------------------------------------------

  !----------------------------------------------------
  ! ELECTROSTATIC or ELECTROMAGNETIC:
  !
  select case (n_field) 

  case (1)

     call send_line('fluctuations         : ELECTROSTATIC')

  case (2)

     call send_line('fluctuations         : EM: (Phi,A_par)')

  case (3)

     call send_line('fluctuations         : EM: (Phi,A_par,B_par)')

  case default

     call catch_error('INVALID: n_field')

  end select
  !----------------------------------------------------

  !----------------------------------------------------
  ! GYRO-AVERAGE METHOD:
  !
  select case (gyro_method) 

  case (1)

     call send_line('gyro_method          : ORIGINAL (INTEGRATED) POLARIZATION')

  case (2)

     call send_line('gyro_method          : GRID-BASED J0^2 POLARIZATION')

  case default

     call catch_error('INVALID: boundary_method')

  end select

  if (gyro_method == 1) then
     mg_dx = m_dx
     ig_dx = i_dx
  else 
     mg_dx = m_gyro
     ig_dx = i_gyro
  endif
  !----------------------------------------------------

  !----------------------------------------------------
  ! FLUCTUATION OUTPUT
  !
  mstr = ''
  if (plot_u_flag == 1) mstr=trim(mstr)//'PHI(+)'
  if (plot_u_flag == 0) mstr=trim(mstr)//'PHI( )'
  if (plot_n_flag == 1) mstr=trim(mstr)//' N(+)'
  if (plot_n_flag == 0) mstr=trim(mstr)//' N( )'
  if (plot_e_flag == 1) mstr=trim(mstr)//' E(+)'
  if (plot_e_flag == 0) mstr=trim(mstr)//' E( )'
  if (plot_v_flag == 1) mstr=trim(mstr)//' V(+)'
  if (plot_v_flag == 0) mstr=trim(mstr)//' V( )'

  call send_line('fluctuation output   : '//mstr)
  !----------------------------------------------------

  !----------------------------------------------------
  ! GRID SPACE METHOD:
  !
  select case (nonuniform_grid_flag) 

  case (0)

     call send_line('radial grid          : UNIFORM')

  case (1)

     if (flat_profile_flag == 1) then
        nonuniform_grid_flag = 0 
        call send_line('radial grid          : UNIFORM [reset]')
     else
        call send_line('radial grid          : NONUNIFORM (not recommended)')
     endif

  case default

     call catch_error('INVALID: nonuniform_grid_flag')

  end select
  !----------------------------------------------------

  !----------------------------------------------------
  ! INTEGRATOR:
  !
  if (electron_method == 2) then

     select case (integrator_method) 

     case (1)

        call send_line('integrator_method    : IMEX2')

     case (2)

        call send_line('integrator_method    : IMEX2+RK4')

     case (3)

        call send_line('integrator_method    : IMEX2+2(RK4)')

     case (4)

        call send_line('integrator_method    : IMEX2+3(RK4)')

     case default

        call catch_error('INVALID: integrator_method')

     end select

  else

     select case (integrator_method) 

     case (1)

        call send_line('integrator_method    : RK2')

     case (2)

        call send_line('integrator_method    : RK4')

     case default

        call catch_error('INVALID: integrator_method')

     end select

  endif
  !----------------------------------------------------

  !-------------------------------------------------------
  ! SPARSE MATRIX PACKAGE:
  !
  select case (sparse_method)

  case (1)

     call send_line('field solve          : UMFPACK (SEQUENTIAL)')

  case (2)

     call send_line('field solve          : MUMPS (DISTRIBUTED)')

  case default

     call catch_error('INVALID: sparse_method')

  end select
  !-------------------------------------------------------

  !----------------------------------------------------
  ! SOURCE METHOD:
  !
  select case (source_flag) 

  case (0)

     call send_line('source_flag          : SOURCE OFF')

  case (1)

     call send_line('source_flag          : SOURCE ON')

     ! There are 4 more partial lumps than full lumps.
     n_lump = 4+n_source

     if (n_lump < 5) then
        call catch_error('INVALID: N_SOURCE TOO SMALL')
     endif

  case default

     call catch_error('INVALID: source_flag')

  end select
  !----------------------------------------------------

  !-------------------------------------------------------
  ! COLLISION FLAG:
  !
  if (collision_flag == 0) then

     call send_line('pitch-angle coll.    : OFF')

  else

     select case (ord_rbf)

     case (3)

        call send_line('pitch-angle coll.    : RBF r^3')

     case (5)

        call send_line('pitch-angle coll.    : RBF r^5')

     case (7)

        call send_line('pitch-angle coll.    : RBF r^7')

     case (9)

        call send_line('pitch-angle coll.    : RBF r^9')

     case (11)

        call send_line('pitch-angle coll.    : RBF r^11')

     case default

        call catch_error('INVALID: ord_rbf')

     end select

  endif

  !-------------------------------------------------------

  !-------------------------------------------------------
  ! KROOK COLLISIONS FLAG:
  !
  if (nu_i_krook > 0.0) then

     krook_flag = 1
     call send_line('ion Krook coll.      : ON for all n')

  else

     krook_flag = 0
     !     call send_line('Ion Krook coll.      : OFF')

  endif
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! ROTATION METHOD:
  !
  select case (rotation_theory_method)

  case (1)
     call send_line('rotation parameters  : CANDY METHOD')
  case (2)
     call send_line('rotation effects     : WALTZ METHOD')
  case default
     call catch_error('INVALID: rotation_theory_method')
  end select
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! NONLINEAR FLAG:
  !
  select case (nonlinear_flag)

  case (0)

     call send_line('nonlinear dynamics   : OFF')

     !-------------------------------------------------------
     ! LINDIFF_METHOD:
     !
     select case (lindiff_method)

     case (1)

        call send_line('linear diffusivity   : OFF')

     case (2)

        call send_line('linear diffusivity   : D/CHI ')

     case (3)

        call send_line('linear diffusivity   : RESPONSE FUNCTION')

     case (4)

        call send_line('linear diffusivity   : phi_squared QL weight form')

     case (5)

        call send_line('linear diffusivity   : local phi_squared QL weight form')

     case (6)

        call send_line('linear diffusivity   : local g_squared QL(1) |g|**2 weight form')

     case (7)

        call send_line('linear diffusivity   : local g_squared QL(2) |phi|**2 weight form')

     case (8)

        call send_line('linear diffusivity   : local q_squared QL(3) |g|*|phi| weight form')

     case default

        call catch_error('INVALID: lindiff_method')

     end select
     !-------------------------------------------------------

  case (1)

     call send_line('nonlinear dynamics   : ON')

     !-------------------------------------------------------
     ! NONLINEAR METHOD:
     !
     select case (nl_method)

     case (1)

        call send_line('nonlinear bracket    : DIRECT')

     case (2)

        call send_line('nonlinear bracket    : REAL FFT')

     case default

        call catch_error('INVALID: nl_method')

     end select
     !-------------------------------------------------------

  case default

     call catch_error('INVALID: nonlinear_flag')

  end select
  !-------------------------------------------------------

  call send_line(separator)


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[select_methods done]'
  endif

end subroutine gyro_select_methods
