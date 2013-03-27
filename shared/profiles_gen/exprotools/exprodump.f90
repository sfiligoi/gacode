!------------------------------------------------------------
! exprodump.f90
!
! PURPOSE:
!  Dump a given EXPRO interface variable as a function of 
!  rmin and rho.
!------------------------------------------------------------

program exprodump

  use EXPRO_interface

  implicit none

  character(len=20) :: var,x
  integer :: i
  real, dimension(3) :: z


  open(unit=1,file='input.exprodump',status='old')
  read(1,*) var
  read(1,*) x
  read(1,*) z(1)
  read(1,*) z(2)
  read(1,*) z(3)
  close(1)

  EXPRO_ctrl_density_method = 1
  EXPRO_ctrl_z(1:3) = z(1:3)
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_alloc('./',1) 
  call EXPRO_read

  if (trim(var) == 'list') then

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

     print *,'vol'
     print *,'volp'
     print *,'cs'
     print *,'rhos'
     print *,'drdrho'

     print *,'grad_r0'
     print *,'ave_grad_r'

     print *,'THIS LIST IS INCOMPLETE.  PLEASE SEND EMAIL TO candy@fusion.gat.com'
     print *,'if you want me to finish the coding.'

  else   

     do i=1,EXPRO_n_exp

        ! Column 1
        select case (trim(x)) 
        case ('none')
           write(*,'(a)',advance='no') ' ' 
        case ('r')
           write(*,10,advance='no') EXPRO_rmin(i)/EXPRO_rmin(EXPRO_n_exp)
        case ('rho')
           write(*,10,advance='no') EXPRO_rho(i)
        case ('psi')
           write(*,10,advance='no') EXPRO_poloidalfluxover2pi(i)/EXPRO_poloidalfluxover2pi(EXPRO_n_exp)
        end select

        ! Column 2
        select case (trim(var))

        case ('rho')
           print 20,EXPRO_rho(i)
        case ('rmin')
           print 20,EXPRO_rmin(i)
        case ('rmaj')
           print 20,EXPRO_rmaj(i)
        case ('q')
           print 20,EXPRO_q(i)
        case ('kappa')
           print 20,EXPRO_kappa(i)
        case ('delta')
           print 20,EXPRO_delta(i)
        case ('te')
           print 20,EXPRO_te(i)
        case ('ne')
           print 20,EXPRO_ne(i)
        case ('z_eff')
           print 20,EXPRO_z_eff(i)
        case ('w0')
           print 20,EXPRO_w0(i)
        case ('flow_mom')
           print 20,EXPRO_flow_mom(i)
        case ('pow_e')
           print 20,EXPRO_pow_e(i)
        case ('pow_i')
           print 20,EXPRO_pow_i(i)
        case ('pow_ei')
           print 20,EXPRO_pow_ei(i)
        case ('zeta')
           print 20,EXPRO_zeta(i)
        case ('flow_beam')
           print 20,EXPRO_flow_beam(i)
        case ('flow_wall')
           print 20,EXPRO_flow_wall(i)
        case ('zmag')
           print 20,EXPRO_zmag(i)
        case ('ptot')
           print 20,EXPRO_ptot(i)
        case ('poloidalfluxover2pi')
           print 20,EXPRO_poloidalfluxover2pi(i)
        case ('bunit')
           print 20,EXPRO_bunit(i)
        case ('s')
           print 20,EXPRO_s(i)
        case ('drmaj')
           print 20,EXPRO_drmaj(i)
        case ('dzmag')
           print 20,EXPRO_dzmag(i)
        case ('sdelta')
           print 20,EXPRO_sdelta(i)
        case ('vol')
           print 20,EXPRO_vol(i)
        case ('volp')
           print 20,EXPRO_volp(i)
        case ('cs')
           print 20,EXPRO_cs(i)
        case ('rhos')
           print 20,EXPRO_rhos(i)
        case ('drdrho')
           print 20,EXPRO_drdrho(i)
        case ('grad_r0')
           print 20,EXPRO_grad_r0(i)
        case ('ave_grad_r')
           print 20,EXPRO_ave_grad_r(i)
        case default
           print 20,0.0
        end select

     enddo

  endif

  call EXPRO_alloc('./',0) 

10 format(t2,1pe12.5)
20 format(1x,1pe12.5)

end program exprodump
