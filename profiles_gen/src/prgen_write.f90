!---------------------------------------------------------
! prgen_write.f90
!
! PURPOSE:
!  Generate input.profiles and input.profiles.extra
!----------------------------------------------------------

subroutine prgen_write

  use prgen_globals
  use expro

  !---------------------------------------------------------------
  implicit none
  !
  integer :: i
  real :: qd
  !---------------------------------------------------------------

  expro_rvbv = 0.0
  expro_ipa  = 0.0

  ! Ensure correct sign of toroidal flux (Bt)
  !
  select case (format_type)

  case (0)
     ! (nothing but gfile)
     expro_torfluxa = -btccw*abs(expro_torfluxa)
  case (1)
     ! onetwo statefile
     expro_torfluxa = -0.5*btccw*abs(onetwo_Btor)*onetwo_rho_grid(nx)**2
     expro_rvbv = onetwo_R0*onetwo_Btor
     expro_ipa = -ipccw*abs(ip_tot)
  case (2)
     ! plasmastate 
     expro_torfluxa = -btccw*abs(plst_phit(nx))/(2*pi)
  case (3)
     ! pfile
     expro_torfluxa = -btccw*abs(peqdsk_torfluxa)
  case (5)
     ! corsica
     expro_torfluxa = -btccw*abs(corsica_torfluxa)
  case (6)
     ! ufile
     expro_torfluxa = -btccw*abs(ufile_torfluxa)
  case (7)
     ! ufile
     expro_torfluxa = -btccw*abs(expro_torfluxa)
  end select

  !-------------------------------------------------------------------------------------
  ! Now we need to write the profile data
  !
  expro_head_original  = '#  *original : '//trim(date)
  expro_head_statefile = '# *statefile : '//trim(file_state)
  expro_head_gfile     = '#     *gfile : '//trim(file_g)
  expro_head_cerfile   = '#   *cerfile : '//trim(file_cer)

  call expro_write('input.gacode')
  print '(a)','INFO: (prgen_write) Wrote input.gacode.'
  !-------------------------------------------------------------------------------------

  !call expro_read('input.gacode')
  !qd = 0.0
  !do i=1,nx
  !   print *,i,expro_vol(i),(pi*expro_rmin(i)**2*expro_kappa(i))*2*pi*expro_rmaj(i)
  !enddo
  !do i=1,nx-1
  !   qd = qd+0.5*(expro_vol(i+1)-expro_vol(i))*(onetwo_qdelt(i+1)+onetwo_qdelt(i))
  !   print *,i,expro_pow_ei(i+1),qd*(-1e-6)
  !enddo
 
end subroutine prgen_write

