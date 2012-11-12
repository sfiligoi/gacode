!------------------------------------------------------------
! exprodump.f90
!
! PURPOSE:
!  Dump a given EXPRO interface variable as a function of 
!  rmin and rho.
!------------------------------------------------------------

program locpargen

  use EXPRO_interface

  implicit none

  character(len=20) :: var
  integer :: i
  real, dimension(3) :: z


  open(unit=1,file='input.exprodump',status='old')
  read(1,*) var
  read(1,*) z(1)
  read(1,*) z(2)
  read(1,*) z(3)
  close(1)

  EXPRO_ctrl_density_method = 1
  EXPRO_ctrl_z(1:3) = z(1:3)
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_signq = (-1)*(-1)
  EXPRO_ctrl_signb = -(-1)
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_alloc('./',1) 
  call EXPRO_read

  select case(trim(var))

  case ('list')

     print *,'rho'
     print *,'rmin'
     print *,'rmaj'
     print *,'q'
     print *,'kappa'

     print *,'delta'
     print *,'te'
     print *,'ne'
     print *,'z_eff'
     print *,'w0'

     print *,'flow_mom'
     print *,'pow_e'
     print *,'pow_i'
     print *,'pow_ei'
     print *,'zeta'

     print *,'flow_beam'
     print *,'flow_wall'
     print *,'zmag'
     print *,'ptot'
     print *,'poloidalfluxpver2pi'

     print *,'bunit'
     print *,'s'
     print *,'drmaj'
     print *,'dzmag'
     print *,'sdelta'

  case ('rho')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_rho(i)
     enddo

  case ('rmin')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_rmin(i)
     enddo

  case ('rmaj')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_rmaj(i)
     enddo

  case ('q')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_q(i)
     enddo

  case ('kappa')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_kappa(i)
     enddo

  case ('delta')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_delta(i)
     enddo

  case ('te')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_te(i)
     enddo

  case ('ne')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_ne(i)
     enddo

  case ('z_eff')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_z_eff(i)
     enddo

  case ('w0')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_w0(i)
     enddo

  case ('flow_mom')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_flow_mom(i)
     enddo

  case ('pow_e')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_pow_e(i)
     enddo

  case ('pow_i')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_pow_i(i)
     enddo

  case ('pow_ei')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_pow_ei(i)
     enddo

  case ('zeta')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_zeta(i)
     enddo

  case ('flow_beam')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_flow_beam(i)
     enddo

  case ('flow_wall')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_flow_wall(i)
     enddo

  case ('zmag')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_zmag(i)
     enddo

  case ('ptot')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_ptot(i)
     enddo

  case ('poloidalfluxover2pi')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_poloidalfluxover2pi(i)
     enddo

  case ('bunit')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_bunit(i)
     enddo

  case ('s')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_s(i)
     enddo

  case ('drmaj')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_drmaj(i)
     enddo

  case ('dzmag')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_dzmag(i)
     enddo

  case ('sdelta')

     do i=1,EXPRO_n_exp
        print 10,EXPRO_rho(i),EXPRO_sdelta(i)
     enddo

  end select

  call EXPRO_alloc('./',0) 

10 format(t2,2(1pe12.5,1x))

end program locpargen
