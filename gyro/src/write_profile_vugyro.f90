!-----------------------------------------------------------
! write_profile_vugyro.f90 [caller make_profiles]
!
! PURPOSE:
!  Write data relevant for post-simulation analysis.
!
! NOTES:
!  See documentation for complete definitions of
!  variables. 
!-----------------------------------------------------------

subroutine write_profile_vugyro(datafile,io)

  use gyro_globals
  use gyro_profile_exp
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !---------------------------------------------------

  select case (output_flag)

  case (1)

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='replace')

        ! Essential scalars
        write(io,20) n_x
        write(io,20) n_theta_section
        write(io,20) n_pass
        write(io,20) n_trap
        write(io,20) n_energy
        write(io,20) n_theta_plot
        write(io,20) n0
        write(io,20) n_n
        write(io,20) d_n
        write(io,20) n_explicit_damp
        write(io,20) nonlinear_flag
        write(io,20) electron_method
        write(io,20) n_field
        write(io,20) n_ion
        write(io,20) n_kinetic
        write(io,20) n_spec
        write(io,20) field_r0_flag
        write(io,20) field_r0_grid
        write(io,20) n_grid_exp
        write(io,20) boundary_method

        ! Basic profile data (** CORRECT IN VUGYRO **)
        write(io,10) r(:)
        write(io,10) q(:)
        write(io,10) r_s(:)
        write(io,10) q_s(:)
        write(io,10) dlntdr_s(:,:)
        write(io,10) dlnndr_s(:,:)
        write(io,10) tem_s(:,:)
        write(io,10) den_s(:,:)
        write(io,10) rmaj_s(:)/r_s(:)
        write(io,10) delta_s(:)
        write(io,10) zeta_s(:)
        write(io,10) kappa_s(:)
        write(io,10) drmaj_s(:)
        write(io,10) shat_s(:)
        write(io,10) s_delta_s(:)
        write(io,10) s_zeta_s(:)
        write(io,10) s_kappa_s(:)
        write(io,10) zmag_s(:)
        write(io,10) dzmag_s(:)
        write(io,10) beta_unit_s(:)
        write(io,10) gamma_e_s(:)
        write(io,10) gamma_p_s(:) 
        write(io,10) mach_s(:) 
        write(io,10) b_unit_s(:)
        write(io,10) dr_eodr(:)
        write(io,10) z_eff_s(:)
        write(io,10) nu_s(n_spec,:)
        write(io,10) w0_s(:)

        write(io,10) box_multiplier

        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)
        write(io,10) 0.0*r(:)

        write(io,10) lambda(ir_norm,:)
        write(io,10) energy(:,1)
        write(io,10) lambda_tp(ir_norm)
        write(io,10) krho_collect(:)
        write(io,10) rhos_norm
        write(io,10) z(:)

        ! Recent additions
        write(io,20) n_theta_plot*n_theta_mult

        ! Number of flux moments
        write(io,20) n_moment

     endif

  end select

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[write_profile_vugyro done]'
  endif

10 format(1pe15.8)
20 format(i4)

end subroutine write_profile_vugyro
