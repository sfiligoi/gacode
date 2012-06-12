subroutine tgyro_write_input

  use mpi
  use tgyro_globals

  implicit none

  integer :: i_ion

  !----------------------------------------------------------------
  ! Trap miscellaneous errors
  !
  if (i_bc < 0) then
     call tgyro_catch_error('ERROR: (TGYRO) Problem with matching radius.')
  endif
  !
  ! Advanced iteration methods cannot do zero iterations:
  !
  if (tgyro_iteration_method >= 4 .and. tgyro_relax_iterations == 0) then
     call tgyro_catch_error('ERROR: (TGYRO) TGYRO_ITERATION_METHOD > 4 requires TGYRO_RELAX_ITERATIONS > 0.')
  endif
  !----------------------------------------------------------------

  if (i_proc_global == 0) then

     open(unit=1,file=trim(runfile),position='append')

     write(1,*) 'Job control'
     write(1,*) 

     if (tgyro_mode == 2) then

        write(1,30) 'TGYRO_STAB_NSEARCH',tgyro_stab_nsearch
        write(1,30) 'TGYRO_STAB_NKY',tgyro_stab_nky
        write(1,20) 'TGYRO_STAB_KYMIN',tgyro_stab_kymin

        ! Ooh the dreaded goto ...

        goto 100

     endif

     !--------------------------------------------------------
     select case (loc_restart_flag)

     case (0)

        write(1,10) 'LOC_RESTART_FLAG','New simulation'

     case (1)

        write(1,10) 'LOC_RESTART_FLAG','Restarting old simulation'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_RESTART_FLAG'

     end select
     !--------------------------------------------------------

     write(1,*)
     write(1,*) 'Iteration control'
     write(1,*) 

     select case (tgyro_iteration_method)

     case (1)

        write(1,10) 'TGYRO_ITERATION_METHOD','Local Residual Minimization (original scheme)'

     case (2)

        write(1,10) 'TGYRO_ITERATION_METHOD','Levenberg-Marquardt Minimization Algorithm'

        if (loc_residual_method == 2) then
           error_flag = 1
           error_msg = 'Error: TGYRO_ITERATION_METHOD=2 not compatible with LOC_RESIDUAL_METHOD=2'
        endif

        write(1,20) '*LM_BOOST (damp inc)',lm_boost
        write(1,20) '*LM_DROP (damp decr)',lm_drop

     case (3)

        write(1,10) 'TGYRO_ITERATION_METHOD','Newton with Linesearch'

     case (4)

        write(1,10) 'TGYRO_ITERATION_METHOD','Block method (serial)'

     case (5)

        write(1,10) 'TGYRO_ITERATION_METHOD','Block method (parallel)'

     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_ITERATION_METHOD'

     end select

     !--------------------------------------------------------
     if (tgyro_iteration_method == 2 .or. tgyro_iteration_method == 3) then

        select case (tgyro_backtrack_method)

        case (0)

           write(1,10) '*TGYRO_BACKTRACK_METHOD','Backtracking Off'

        case (1)

           write(1,10) '*TGYRO_BACKTRACK_METHOD','Golden Ratio Search'

        case (2)

           write(1,10) '*TGYRO_BACKTRACK_METHOD','Optimal Parabolic Backtrack'

        case (3)

           write(1,10) '*TGYRO_BACKTRACK_METHOD','3'

        case default

           error_flag = 1
           error_msg = 'Error: TGYRO_BACKTRACK_METHOD'

        end select

     endif
     !--------------------------------------------------------

     if (tgyro_global_newton_flag == 0) then
        write(1,10) 'TGYRO_GLOBAL_NEWTON_FLAG','Block-diagonal flux Jacobian'
     else
        write(1,10) 'TGYRO_GLOBAL_NEWTON_FLAG','Exact flux Jacobian'
     endif

     !--------------------------------------------------------
     select case (loc_residual_method)

     case (1)

        write(1,10) 'LOC_RESIDUAL_METHOD','(f-g)^2/(f^2+g^2)'

     case (2)

        write(1,10) 'LOC_RESIDUAL_METHOD','abs(f-g)'

     case (3)

        write(1,10) 'LOC_RESIDUAL_METHOD','(f-g)^2'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_RESIDUAL_METHOD'

     end select
     !--------------------------------------------------------

     write(1,20) 'LOC_DX (Jacobian dx)',loc_dx
     write(1,20) 'LOC_DX_GYRO (GYRO Jacobian dx)',loc_dx_gyro
     write(1,20) 'LOC_DX_MAX (maximum dx)',loc_dx_max
     write(1,20) 'LOC_RELAX (conv. relaxation)',loc_relax

     write(1,*)
     write(1,*) 'Scenario control'
     write(1,*) 

     !--------------------------------------------------------
     select case (loc_ti_feedback_flag)

     case (0)

        write(1,10) 'LOC_TI_FEEDBACK_FLAG','Ti evolution OFF'

     case (1)

        write(1,10) 'LOC_TI_FEEDBACK_FLAG','Ti evolution ON'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_TI_FEEDBACK_FLAG'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_te_feedback_flag)

     case (0)

        write(1,10) 'LOC_TE_FEEDBACK_FLAG','Te evolution OFF'

     case (1)

        write(1,10) 'LOC_TE_FEEDBACK_FLAG','Te evolution ON'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_TE_FEEDBACK_FLAG'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_ne_feedback_flag)

     case (0)

        write(1,10) 'LOC_NE_FEEDBACK_FLAG','ne evolution OFF'

     case (1)

        write(1,10) 'LOC_NE_FEEDBACK_FLAG','ne evolution ON'

        if (loc_residual_method == 1) then
           error_flag = 1
           error_msg = 'Error: LOC_RESIDUAL_METHOD=1 not compatible with LOC_NE_FEEDBACK_FLAG=1.'
        endif

        select case (loc_pflux_method)

        case (1)

           write(1,10) 'LOC_PFLUX_METHOD','gamma=0'

        case (2)

           write(1,10) 'LOC_PFLUX_METHOD','gamma=beam'

        case default

           error_flag = 1
           error_msg = 'Error: LOC_PFLUX_METHOD'

        end select

     case default

        error_flag = 1
        error_msg = 'Error: LOC_NE_FEEDBACK_FLAG'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_er_feedback_flag)

     case (0)

        write(1,10) 'LOC_ER_FEEDBACK_FLAG','Er evolution OFF'

     case (1)

        write(1,10) 'LOC_ER_FEEDBACK_FLAG','Er evolution ON'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_ER_FEEDBACK_FLAG'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_scenario)

     case (1)

        write(1,10) 'LOC_SCENARIO','Experimental, static exchange'

     case (2)

        write(1,10) 'LOC_SCENARIO','Experimental, dynamic exchange'

     case (3)

        write(1,10) 'LOC_SCENARIO','Reactor with input aux. power'

     case (4)

        write(1,10) 'LOC_SCENARIO','Reactor with model aux. power'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_SCENARIO'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_sawtooth_model)

     case (1)

        write(1,10) 'LOC_SAWTOOTH_MODEL','No corrections for sawteeth'

     case (2)

        write(1,10) 'LOC_SAWTOOTH_MODEL','Enhanced qe for q < 1'

     case (3)

        write(1,10) 'LOC_SAWTOOTH_MODEL','Enhanced qi and qe for q < 1.'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_SAWTOOTH_MODEL'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (loc_neo_method)

     case (1)

        if (loc_chang_hinton == 1) then
           write(1,10) 'LOC_NEO_METHOD','Hinton-Hazeltine theory (Chang-Hinton chii)'
        else
           write(1,10) 'LOC_NEO_METHOD','Hinton-Hazeltine theory'
        endif

     case (2)

        write(1,10) 'LOC_NEO_METHOD','NEO code'
        if (loc_n_ion >= 4) then
           write(1,10) '        * INFO','Using reduced energy resolution to cope with so many ions.'
        endif

     case default

        error_flag = 1
        error_msg = 'Error: LOC_NEO_METHOD'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (neo_gv_flag)

     case (0)

        write(1,10) 'NEO_GV_FLAG','OFF (no gyroviscosity)'

     case (1)

        write(1,10) 'NEO_GV_FLAG','ON (include gyroviscosity)'

     case default

        error_flag = 1
        error_msg = 'Error: NEO_GV_FLAG'

     end select
     !--------------------------------------------------------

     !--------------------------------------------------------
     select case (tgyro_tglf_revision)

     case (1)

        write(1,10) 'TGYRO_TGLF_REVISION','APS07'

     case (2)

        write(1,10) 'TGYRO_TGLF_REVISION','TGLF-09'

     case (3)

        write(1,10) 'TGYRO_TGLF_REVISION','in development'

     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_TGLF_REVISION'

     end select
     !--------------------------------------------------------

     write(1,20) 'Pivot radius',r(i_bc)/r_min

100  continue 

     write(1,*)
     write(1,*) 'Physics parameters'
     write(1,*) 

     write(1,20) 'LOC_RMIN',r(2)/r_min
     write(1,20) 'LOC_RMAX',r(n_r)/r_min
     write(1,20) 'LOC_NU_SCALE',loc_nu_scale
     write(1,20) 'LOC_BETAE_SCALE',loc_betae_scale

     !--------------------------------------------------------
     select case (loc_zeff_flag)

     case (0)

        write(1,10) 'LOC_ZEFF_FLAG','Z_eff=1'

     case (1)

        write(1,10) 'LOC_ZEFF_FLAG','Z_eff from data'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_ZEFF_FLAG'

     end select
     !--------------------------------------------------------
     !--------------------------------------------------------
     select case (loc_circ_flag)

     case (0)

        write(1,10) 'LOC_CIRC_FLAG','Use standard Miller expressions'

     case (1)

        write(1,10) 'LOC_CIRC_FLAG','Use Miller circle'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_CIRC_FLAG'

     end select
     !--------------------------------------------------------
     !--------------------------------------------------------
     select case (tgyro_rotation_theory_method)

     case (1)

        write(1,10) 'TGYRO_ROTATION_THEORY_METHOD','Candy method'

     case (2)

        write(1,10) 'TGYRO_ROTATION_THEORY_METHOD','Waltz method'

     case default

        error_flag = 1
        error_msg = 'Error: TGYRO_ROTATION_THEORY_METHOD'

     end select
     !--------------------------------------------------------

     if (loc_betae_scale > 0.0) then
        write(1,10) 'Fluctuations','Electromagnetic'
     else
        write(1,10) 'Fluctuations','Electrostatic'
     endif


     !--------------------------------------------------------
     select case (loc_quasineutral_flag)

     case (0)

        write(1,10) 'LOC_QUASINEUTRAL_FLAG','Densities unaltered'

     case (1)

        write(1,10) 'LOC_QUASINEUTRAL_FLAG','Forced quasineutrality'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_QUASINEUTRAL_FLAG'

     end select
     !--------------------------------------------------------
     !--------------------------------------------------------
     select case (loc_lock_profile_flag)

     case (0)

        write(1,10) 'LOC_LOCK_PROFILE_FLAG','Initial profiles from gradients'

     case (1)

        write(1,10) 'LOC_LOCK_PROFILE_FLAG','Initial profiles from data'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_LOCK_PROFILE_FLAG'

     end select
     !--------------------------------------------------------
      !--------------------------------------------------------
     select case (loc_evolve_grad_only_flag)

     case (0)

        write(1,10) 'LOC_EVOLVE_GRAD_ONLY_FLAG','Self-consistently evolve profiles with gradients'

     case (1)

        write(1,10) 'LOC_EVOLVE_GRAD_ONLY_FLAG','Evolve gradients but *NOT* corresponding profiles'

     case default

        error_flag = 1
        error_msg = 'Error: LOC_LOCK_PROFILE_FLAG'

     end select
     !--------------------------------------------------------
    do i_ion=1,loc_n_ion
        write(1,20) 'ion '//trim(ion_tag(i_ion))//' mass',mi_vec(i_ion)
        write(1,20) 'ion '//trim(ion_tag(i_ion))//' charge ',zi_vec(i_ion)
     enddo

     write(1,*)'-----------------------------------------------------------------'

     close(1)

  endif

  call MPI_BCAST(error_flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(error_msg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (error_flag == 1) then
     call tgyro_catch_error(trim(error_msg))
  endif

10 format(t2,a,t33,':',t35,a)
20 format(t2,a,t33,':',t35,f8.4)
30 format(t2,a,t33,':',t35,i4)

end subroutine tgyro_write_input

