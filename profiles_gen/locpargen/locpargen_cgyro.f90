subroutine locpargen_cgyro

  use locpargen_globals
  use EXPRO_interface
  use EXPRO_locsim_interface

  implicit none

  open(unit=1,file='input.cgyro.locpargen',status='replace')
  write(1,10) 'RMIN=',r0
  write(1,10) 'RMAJ=',rmaj_loc
  write(1,*)
  write(1,10) 'SHIFT=',shift_loc
  write(1,10) 'ZMAG=',zmag_loc
  write(1,10) 'DZMAG=',dzmag_loc
  write(1,10) 'Q=',q_loc
  write(1,10) 'S=',s_loc
  write(1,10) 'KAPPA=',kappa_loc
  write(1,10) 'S_KAPPA=',s_kappa_loc
  write(1,10) 'DELTA=',delta_loc
  write(1,10) 'S_DELTA=',s_delta_loc
  write(1,10) 'ZETA=',zeta_loc
  write(1,10) 'S_ZETA=',s_zeta_loc
  write(1,*)
  write(1,10) 'GAMMA_E=',gamma_e_loc*a/cs_loc
  write(1,10) 'GAMMA_P=',gamma_p_loc*a/cs_loc
  write(1,10) 'MACH=',mach_loc/cs_loc
  write(1,*)
  write(1,10) 'NU_EE=',nu_ee*a/cs_loc
  write(1,10) 'BETAE_UNIT=',betae_unit

  !---------------------------------------------------------
  ! Species data
  
  write(1,*)
  write(1,11) 'N_SPECIES=',ise

  do is=1,ise
     write(1,*)
     write(1,11) 'Z_'//tag(is)//'=',int(z_loc(is))
     ! Deuteron mass normalization
     write(1,10) 'MASS_'//tag(is)//'=',mass_loc(is)/2.0
     write(1,10) 'DENS_'//tag(is)//'=',dens_loc(is)/dens_loc(ise)
     write(1,10) 'TEMP_'//tag(is)//'=',temp_loc(is)/temp_loc(ise)
     write(1,10) 'DLNNDR_'//tag(is)//'=',dlnndr_loc(is)
     write(1,10) 'DLNTDR_'//tag(is)//'=',dlntdr_loc(is)
     write(1,10) 'SDLNNDR_'//tag(is)//'=',sdlnndr_loc(is)
     write(1,10) 'SDLNTDR_'//tag(is)//'=',sdlntdr_loc(is)
  enddo
  close(1)

10 format(a,sp,1pe12.5)
11 format(a,i0)

end subroutine locpargen_cgyro
