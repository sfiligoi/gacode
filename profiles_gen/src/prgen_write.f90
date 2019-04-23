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
  !---------------------------------------------------------------

  expro_header(1) = '#  original : '//trim(date)
  expro_header(2) = '# statefile : '//trim(file_state)
  expro_header(3) = '#     gfile : '//trim(file_g)
  expro_header(4) = '#   cerfile : '//trim(file_cer)

  !write(1,'(a)')     '#                      IONS :  Name       Z    Mass'
  !do i=1,expro_n_ion
  !   write(1,'(a,t32,a,t41,i3,t47,f5.1,t53,a)') '#',&
  !        expro_name(i),nint(expro_z(i)),expro_mass(i),expro_type(i)
  !enddo
  !write(1,'(a,i2)')  '#                     IPCCW : ',ipccw
  !write(1,'(a,i2)')  '#                     BTCCW : ',btccw

  expro_rvbv  = 0.0
  expro_ip_exp = 0.0

  ! Ensure correct sign of toroidal flux (Bt)
  !
  select case (format_type)

  case (0)
     ! (nothing but gfile)
     expro_torfluxa = -btccw*abs(expro_torfluxa)
  case (1)
     ! onetwo statefile
     expro_torfluxa = -0.5*btccw*abs(onetwo_Btor)*onetwo_rho_grid(nx)**2
     expro_rvbv  = onetwo_R0*onetwo_Btor
     expro_ip_exp = -ipccw*abs(ip_tot)
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
  call expro_write_all('input.gacode')
  print '(a)','INFO: (prgen_write) Wrote input.gacode.'
  !-------------------------------------------------------------------------------------

end subroutine prgen_write

