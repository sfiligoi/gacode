!---------------------------------------------------------
! prgen_write.f90
!
! PURPOSE:
!  Generate input.profiles and input.profiles.extra
!----------------------------------------------------------

subroutine prgen_write

  use prgen_globals
  use expro

  implicit none
  
  expro_rcentr = rcentr
  expro_bcentr = -btccw*abs(bcentr)
  expro_current = -ipccw*abs(current)

  ! Ensure correct sign of toroidal flux (Bt)
  !
  select case (format_type)

  case (0,1,2,3,4,5,6)
     expro_torfluxa = -btccw*abs(torfluxa)
  case (7,8)
     ! GACODE/LEGACY
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
  
end subroutine prgen_write

