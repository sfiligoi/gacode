!-----------------------------------------------------------
! write_geometry_arrays.f90 [caller make_geometry_arrays]
!
! PURPOSE:
!  Write 2-D data for Miller geometry functions.
!
! NOTES:
!  Theta resolution is determined by parameters 
!  N_THETA_PLOT and N_THETA_MULT.
!-----------------------------------------------------------

subroutine write_geometry_arrays(datafile,io)

  use gyro_globals
  use math_constants
  use GEO_interface

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer :: n_fine
  integer :: io_mode
  real :: theta
  real :: dr
  !---------------------------------------------------

  if (i_proc == 0) then

     select case (output_flag)

     case (0)

        io_mode = 0

     case (1)

        if (io /= 6) then
           io_mode = 1
        else
           io_mode = 0
        endif

     end select

     !------------------------------------------------------------------
     ! Output to unit io:
     !
     if (io_mode == 1) then
        open(unit=io,file=datafile,status='replace')
     endif

     n_fine = n_theta_plot*n_theta_mult

     do i=1,n_x

        if (flat_profile_flag == 0) then

           ! All profiles are global and radial variation 
           ! is consistent

           GEO_rmin_in      = r_s(i)
           GEO_rmaj_in      = rmaj_s(i)
           GEO_drmaj_in     = drmaj_s(i)
           GEO_zmag_in      = zmag_s(i)
           GEO_dzmag_in     = dzmag_s(i)
           GEO_q_in         = q_s(i)
           GEO_s_in         = shat_s(i)
           GEO_kappa_in     = kappa_s(i)
           GEO_s_kappa_in   = s_kappa_s(i)
           GEO_delta_in     = delta_s(i)
           GEO_s_delta_in   = s_delta_s(i)
           GEO_zeta_in      = zeta_s(i)
           GEO_s_zeta_in    = s_zeta_s(i)
           GEO_beta_star_in = beta_star_s(i)

        else

           ! Profiles are flat and so some parameters need
           ! to be linearly extrapolated.

           dr = r(i)-r(ir_norm)

           GEO_rmin_in  = r(i)
           GEO_rmaj_in  = rmaj_s(ir_norm)+drmaj_s(ir_norm)*dr
           GEO_drmaj_in = drmaj_s(ir_norm)
           GEO_zmag_in  = zmag_s(ir_norm)+dzmag_s(ir_norm)*dr
           GEO_dzmag_in = dzmag_s(ir_norm)
           GEO_q_in     = q(i)
           GEO_s_in     = shat_s(ir_norm)
           GEO_kappa_in = kappa_s(ir_norm)+&
                      kappa_s(ir_norm)*s_kappa_s(ir_norm)/r(ir_norm)*dr
           GEO_s_kappa_in = s_kappa_s(ir_norm)
           GEO_delta_in   = delta_s(ir_norm)+s_delta_s(ir_norm)/r(ir_norm)*dr
           GEO_s_delta_in = s_delta_s(ir_norm)
           GEO_zeta_in    = zeta_s(ir_norm)+s_zeta_s(ir_norm)/r(ir_norm)*dr
           GEO_s_zeta_in  = s_zeta_s(ir_norm)
           GEO_beta_star_in = beta_star_s(ir_norm)

        endif

        GEO_fourier_in(:,:) = a_fourier_geo_s(:,0:n_fourier_geo,i)
        call GEO_do()

        do j=1,n_fine

           theta = -pi+(j-1)*pi_2/n_fine

           ! Test for special case
           if (n_fine == 1) theta = 0.0

           call GEO_interp(theta)

           if (io_mode == 1) then
              write(io,10) GEO_nu
              write(io,10) GEO_gsin
              write(io,10) GEO_gcos1
              write(io,10) GEO_gcos2
              write(io,10) GEO_usin
              write(io,10) GEO_ucos
              write(io,10) GEO_b
              write(io,10) GEO_g_theta
              write(io,10) GEO_grad_r
              write(io,10) GEO_gq
              write(io,10) GEO_captheta
           endif

        enddo ! j

     enddo ! i

     if (io_mode == 1) close(io)
     !------------------------------------------------------------------

  endif ! i_proc=0 test

10 format(12(es11.4,1x))

end subroutine write_geometry_arrays
