!-----------------------------------------------------------
! write_profile_sim.f90 [caller make_profiles]
!
! PURPOSE:
!  Write all potentially relevant/useful profile
!  quantities on the SIMULATION GRID (_s).
!
! NOTES:
!  See documentation for complete definitions of
!  variables. 
!-----------------------------------------------------------

subroutine write_profile_sim(datafile,io)

  use gyro_globals
  use math_constants
  use GEO_interface

  !---------------------------------------------------
  implicit none
  !
  real :: gamma_t_s(n_x)
  real :: gamma_n_s(n_x)
  real :: d_f(n_x)
  real :: hbfac
  real :: vol(n_x),volp(n_x)
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  character (len=12), dimension(7) :: tag
  character (len=84) :: sep
  !---------------------------------------------------

  select case (output_flag)

  case (1)

     if (i_proc == 0) then

        sep = '-----------------------------------------------------------------------------------'

        open(unit=io,file=datafile,status='replace')

        tag(1) = 'r'
        tag(2) = 'dlntdr_s'
        tag(3) = 'dlnndr_s'
        tag(4) = 'tem_s'
        tag(5) = 'den_s'
        tag(6) = 'alpha_s'
        tag(7) = '(null)'

        do is=1,n_spec
           write(io,'(a)') sep
           if (is < n_spec) then
              write(io,'(t2,a,i2)') '=> ion #',is
           else
              write(io,'(a,i2)') '=> electrons'
           endif
           write(io,10) tag(:)
           write(io,'(a)') sep
           do i=1,n_x
              write(io,20) r(i),&
                   dlntdr_s(is,i),&
                   dlnndr_s(is,i),&
                   tem_s(is,i),&
                   den_s(is,i),&
                   alpha_s(is,i),&
                   0.0
           enddo ! i
        enddo ! is

        tag(1) = 'r'
        tag(2) = 'q_s'
        tag(3) = 'shat_s'
        tag(4) = 'kappa_s'
        tag(5) = 's_kappa_s'
        tag(6) = 'delta_s'
        tag(7) = 's_delta_s'

        write(io,'(a)') sep
        write(io,'(t2,a)') '=> Geometry'
        write(io,10) tag(:)
        write(io,'(a)') sep
        do i=1,n_x
           write(io,20) r(i),&
                q_s(i),&
                shat_s(i),&
                kappa_s(i),&
                s_kappa_s(i),&
                delta_s(i),&
                s_delta_s(i)
        enddo ! i

        tag(1) = 'r'
        tag(2) = 'R/r'
        tag(3) = 'drmaj_s'
        tag(4) = 'b_unit_s'
        tag(5) = 'beta_unit_s'
        tag(6) = 'rhosda_s'
        tag(7) = '(null)'

        write(io,'(a)') sep
        write(io,'(t2,a)') '=> Geometry'
        write(io,10) tag(:)
        write(io,'(a)') sep
        do i=1,n_x
           write(io,20) r(i),&
                rmaj_s(i)/r_s(i),&
                drmaj_s(i),&
                b_unit_s(i),&
                beta_unit_s(i),&
                rhosda_s(i),&
                0.0
        enddo ! i

        tag(1) = 'r'
        tag(2) = 'r_s'
        tag(3) = 'r_e'
        tag(4) = 'dr_eodr'
        tag(5) = 'q'
        tag(6) = '(null)'
        tag(7) = '(null)'

        write(io,'(a)') sep
        write(io,'(t2,a)') '=> Geometry and Flows'
        write(io,10) tag(:)
        write(io,'(a)') sep
        do i=1,n_x
           write(io,20) r(i),&
                r_s(i),&
                r_e(i),&
                dr_eodr(i),&
                q(i),&
                0.0,&
                0.0
        enddo ! i

        tag(1) = 'r'
        tag(2) = 'gamma_e_s'
        tag(3) = 'gamma_p_s'
        tag(4) = 'mach_s'
        tag(5) = 'z_eff_s'
        tag(6) = 'nu_s'
        tag(7) = '(null)'

        write(io,'(a)') sep
        write(io,'(t2,a)') '=> Geometry and Flows'
        write(io,10) tag(:)
        write(io,'(a)') sep
        do i=1,n_x
           write(io,20) r(i),&
                gamma_e_s(i),&
                gamma_p_s(i),&
                mach_s(i),&
                z_eff_s(i),&
                nu_s(n_spec,i),&
                0.0
        enddo ! i

        ! Various shearing rate diagnostics:

        tag(1) = 'r'
        tag(2) = 'gamma_t_s'
        tag(3) = 'gamma_n_s'
        tag(4) = 'gamma_eb_s'
        tag(5) = 'gamma_HB'
        tag(6) = 'w0_s'
        tag(7) = 'w0p_s'

        ! SPECIES: Diagnostics below based on ion species 1. 

        !------------------------------------------------
        call bound_deriv(d_f,q(:)/r(:)*dlntdr_s(1,:),r,n_x)
        !
        ! d_f = (q/r * dlnTdr)'
        !
        gamma_t_s(:) = rhos_norm*r(:)/q(:)*d_f(:)
        !------------------------------------------------

        !--------------------------------------------------------
        call bound_deriv(d_f,q(:)/r(:)*tem_s(1,:)*dlnndr_s(1,:),r,n_x)
        !
        ! d_f = (q/r * T * dlnndr)'
        !
        ! ExB shear due to diamagnetic doppler
        !
        gamma_n_s(:) = rhos_norm*r(:)/q(:)*d_f(:)
        !-------------------------------------------------------

        write(io,'(a)') sep
        write(io,'(t2,a)') '=> Geometry and Flows'
        write(io,10) tag(:)
        write(io,'(a)') sep
        do i=1,n_x

           call gyro_to_geo(i)
           call GEO_interp(0.0)
           hbfac = GEO_grad_r/GEO_gq
           vol(i) = GEO_volume
           volp(i) = GEO_volume_prime

           write(io,20) r(i),&
                gamma_t_s(i),&
                gamma_n_s(i),&
                gamma_e_s(i),&
                gamma_e_s(i)*hbfac,&
                w0_s(i),&
                w0p_s(i)

        enddo ! i

        ! Various shearing rate diagnostics:

        tag(1) = 'r'
        tag(2) = '(null)'
        tag(3) = '(null)'
        tag(4) = '(null)'
        tag(5) = 'zmag_s'
        tag(6) = 'V (volume)'
        tag(7) = 'dV/dr'

        write(io,'(a)') sep
        write(io,'(t2,a)') '=> Geometry and Flows'
        write(io,10) tag(:)
        write(io,'(a)') sep
        do i=1,n_x
           write(io,20) r(i),&
                0.0,&
                0.0,&
                0.0,&
                zmag_s(i),&
                vol(i),&
                volp(i)
        enddo ! i

        close(io)

     endif

  end select

10 format(t2,7(a))
20 format(7(1pe11.4,1x))

end subroutine write_profile_sim
