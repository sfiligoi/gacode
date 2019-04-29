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
  real,dimension(nx) :: sr,sv
  real :: rat
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

  call expro_read('input.gacode')
!  sv(1) = (expro_pow_ei(2)-expro_pow_ei(1))/(expro_vol(2)-expro_vol(1))
!  sv(2) = sv(1)
!  sv(2) = 0.25*(expro_pow_ei(2)-expro_pow_ei(1))/(expro_vol(2)-expro_vol(1))+&
!       0.25*(expro_pow_ei(3)-expro_pow_ei(2))/(expro_vol(3)-expro_vol(2))
!  sv(1) = -sv(2)+2*(expro_pow_ei(2)-expro_pow_ei(1))/(expro_vol(2)-expro_vol(1))

  sv(1) = -1.88e-1
  sv(2) = -1.87e-1
  do i=2,nx-1
     sv(i+1) = -sv(i)+2*(expro_pow_ei(i+1)-expro_pow_ei(i))/(expro_vol(i+1)-expro_vol(i))
  enddo

  sr(1) = (expro_pow_ei(2)-expro_pow_ei(1))/(expro_vol(2)-expro_vol(1))
  sr(2) = sr(1)
  do i=2,nx-1
     rat = (expro_pow_ei(i+1)-expro_pow_ei(i))/(expro_vol(i+1)-expro_vol(i))
     sr(i+1) = rat
  enddo

  open(unit=1,file='out',status='replace')
  do i=1,nx
     print '(i3,2x,3(1pe12.5,1x))',i,sv(i),-1e-6*onetwo_qdelt(i)
     write(1,'(4(1pe12.5,1x))') rho(i),-1e6*sv(i),onetwo_qdelt(i),-1e6*sr(i)
  enddo
  close(1)
 
end subroutine prgen_write

