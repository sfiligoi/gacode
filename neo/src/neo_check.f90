subroutine neo_check

  use neo_globals

  implicit none

  integer :: ir, is, flag

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
  endif

  !-----------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_theta,2) == 0) then 
     call neo_error('ERROR: (NEO) n_theta must be odd')
     return
  endif
  !
  if (n_species > 11) then
     call neo_error('ERROR: (NEO) max n_species is 6')
     return
  endif

  if (rho_in < 0) then
     call neo_error('ERROR: (NEO) rho_unit must be positive')
     return
  endif

  !-----------------------------------------------------------
  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_neoout,*) 'CONTROL PARAMETERS'
     write(io_neoout,*) '------------------'
  end if
  !------------------------------------------------------------

  ! Simulation model
  select case (sim_model)
  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','THEORY'
     endif
  case(1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','NUMERICAL (with theory and nclass)'
     endif
  case(2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','NUMERICAL (with theory)'
     endif
  case(3)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','THEORY (with nclass)'
     endif
  case(4)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','NEURAL NETWORK'
     endif
  case(5)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sim_model','THEORY (modified Sauter)'
     endif
  case default   
     call neo_error('ERROR: (NEO) invalid sim_model')
     return
  end select

  ! Electron model
  select case (ae_flag)
  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'ae_flag','KINETIC ELECTRONS'
     endif
  case(1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'ae_flag','ADIABATIC ELECTRONS'
     endif
  case default
     call neo_error('ERROR: (NEO) invalid ae_flag')
     return
  end select
     
  ! Collision model
  !
  select case (collision_model)  

  case (1) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','CONNOR'
     endif

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','REDUCED HIRSHMAN-SIGMAR'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','FULL HIRSHMAN-SIGMAR'
     endif

  case (4) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','FULL LINEARIZED FOKKER-PLANCK'
     end if

  case (5) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','FULL LINEARIZED FOKKER-PLANCK WITH AD-HOC FIELD PARTICLE TERMS' 
     end if

  case default

     call neo_error('ERROR: (NEO) invalid collision_model')
     return

  end select

  select case (coll_uncoupledei_model) 

  case (0) 
     ! full e-i coupling

  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','Reduced e-i collisional coupling'
     endif

  case (2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','Reduced e-i collisional coupling'
     endif

  case default
     call neo_error('ERROR: (NEO) invalid coll_uncoupledei_model')
     return

  end select

  select case (coll_uncoupledaniso_model) 

  case (0) 
     ! full aniso collision coupling

  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'collision_model','Uncoupled aniso collisions'
     endif

  case default
     call neo_error('ERROR: (NEO) invalid coll_uncoupledaniso_model')
     return

  end select

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Equilibrium model
  !
  if(profile_model == 1) then
     select case (equilibrium_model)  


     case (0) 

        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'equilibrium_model','S-ALPHA'
        end if

     case (1) 

        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'equilibrium_model','LARGE-ASPECT-RATIO'
        end if

     case (2) 

        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'equilibrium_model','MILLER'
        end if

     case default

        call neo_error('ERROR: (NEO) equilibrium_model invalid')
        return

     end select
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Profile model
  !
  select case (profile_model)  
  case (1) 
     if(n_radial > 1) then
        call neo_error('ERROR: (NEO) profile_model=1 must be run with n_radial = 1')
        return
     endif
     ir=1
     do is=1,n_species
        if (dens(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) density must be positive')
           return
        endif
        if (temp(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) temperature must be positive')
           return
        endif
        if (nu(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) collision frequency must be positive')
           return
        endif
        if (abs(z(is))< epsilon(0.0)) then
           call neo_error('ERROR: (NEO) charge must be non-zero')
           return
        endif
     enddo

     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'profile_model','LOCAL'
     endif

  case (2) 

     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'profile_model','GLOBAL PROFILE'
     endif

     select case (profile_erad0_model)
     case(0)
        if (silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'profile_erad0_model','ERAD0 NOT INCLUDED'
        endif
     case (1)
        if (silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'profile_erad0_model','ERAD0 INCLUDED'
        endif
     case default
        call neo_error('ERROR: (NEO) invalid profile_erad0_model')
        return
     end select

  case default

     call neo_error('ERROR: (NEO) invalid profile_model')
     return

  end select

  !-----------------------------------------------------------
  ! Sign of B checks
  if (sign_q > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_q','POSITIVE'
     endif
  else
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_q','NEGATIVE'
     endif
  endif

  if (sign_bunit > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_bunit','POSITIVE (BT CW)'
     endif
  else
     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'sign_bunit','NEGATIVE (BT CCW)'
     endif
  endif


  !-----------------------------------------------------------
  ! Rotation parameter checks
  !
  select case (rotation_model)     
  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'rotation model','ROTATION EFFECTS NOT INCLUDED'
     endif
  case (2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'rotation model','ROTATION EFFECTS INCLUDED'
     endif
  case default
     call neo_error('ERROR: (NEO) invalid rotation_model')
     return
  end select

  select case (spitzer_model)
  case(0)
     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'spitzer_model:','OFF (NEOCLASSICAL TRANSPORT)'
     endif
  case (1)
     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'spitzer_model:','ON (SOLVE SPITZER PROBLEM)'
     endif
  case default
     call neo_error('ERROR: (NEO) invalid spitzer_model')
     return
  end select

  select case(threed_model)
  case(0)
     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'threed_model:','AXISYMMETRIC EQUILIBRIUM'
     endif
  case (1)
     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'threed_model:','NON-AXISYMMETRIC EQUILIBRIUM (LE3)'
     endif
     if(profile_model == 2) then
        call neo_error('ERROR: (NEO) threed_model not available with global profiles')
     endif

     select case(threed_exb_model)
     case(0)
        if (silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'threed_exb_model:','NO HIGHER-ORDER V_EXB DRIFT'
        endif
     case (1)
        if (silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'threed_exb_model:','HIGHER-ORDER V_EXB DRIFT'
        endif
     case default
        call neo_error('ERROR: (NEO) invalid threed_exb_model')
        return
     end select

     select case(threed_drift_model)
     case(0)
        if (silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'threed_drift_model:','NO HIGHER-ORDER V_DRIFT'
        endif
     case (1)
        if (silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,30) 'threed_drift_model:','HIGHER-ORDER V_DRIFT'
        endif
     case default
        call neo_error('ERROR: (NEO) invalid threed_drift_model')
        return
     end select

  case default
     call neo_error('ERROR: (NEO) invalid threed_model')
     return
  end select

  !-----------------------------------------------------------
  ! Anisotropic species checks
  !
  flag = 0
  do is=1,n_species
     select case (aniso_model(is))     
     case (1)
     case (2)
        flag=1
     case default
        call neo_error('ERROR: (NEO) invalid aniso_model')
        return
     end select
  enddo

  if (flag == 1) then
     if (profile_model == 2) then
        call neo_error('ERROR: (NEO) aniso_model not available with global profiles')
     endif
     if (rotation_model == 1) then
        call neo_error('ERROR: (NEO) aniso_model requires rotation_model=2')
     endif

     if (silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,30) 'aniso model:','ANISOTROPIC SPECIES INCLUDED (with poloidal asymmetry)'
     endif
  endif

  !------------------------------------------------------------

  if (silent_flag == 0 .and. i_proc == 0) then

     write(io_neoout,*)
     write(io_neoout,*) 'GRID DIMENSIONS'
     write(io_neoout,*) '---------------'
     write(io_neoout,10) 'n_radial',n_radial
     write(io_neoout,10) 'n_energy',n_energy
     write(io_neoout,10) 'n_xi',n_xi
     if (threed_model == 0) then
        write(io_neoout,10) 'n_theta',n_theta
     else
        write(io_neoout,10) 'n_theta (modes)',n_tptheta
        write(io_neoout,10) 'n_phi   (modes)',n_tpvarphi
        write(io_neoout,10) 'n_tp',tpmatsize
        write(io_neoout,"(t2,a,t19,': ',i5)") 'n_row',n_row
     endif

     write(io_neoout,*)
     write(io_neoout,*) 'PHYSICS PARAMETERS'
     write(io_neoout,*) '------------------'
     do ir=1,n_radial
        write(io_neoout,20) '      r/a:',r(ir)
        write(io_neoout,20) '      R/a:',rmaj(ir), '    shift:',shift(ir)
        write(io_neoout,20) '        q:',q(ir),    '        s:',shear(ir)
        write(io_neoout,20) '     zmag:',zmag(ir), '    dzmag:',s_zmag(ir)
        write(io_neoout,20) '    kappa:',kappa(ir),'  s_kappa:',s_kappa(ir)
        write(io_neoout,20) '    delta:',delta(ir),'  s_delta:',s_delta(ir)
        write(io_neoout,20) '     zeta:',zeta(ir), '   s_zeta:',s_zeta(ir)
        write(io_neoout,20) '     sin3:',shape_sin(3,ir), '   s_sin3:',shape_s_sin(3,ir)
        write(io_neoout,20) '     sin4:',shape_sin(4,ir), '   s_sin4:',shape_s_sin(4,ir)
        write(io_neoout,20) '     sin5:',shape_sin(5,ir), '   s_sin5:',shape_s_sin(5,ir)
        write(io_neoout,20) '     sin6:',shape_sin(6,ir), '   s_sin6:',shape_s_sin(6,ir)
        write(io_neoout,20) '     cos0:',shape_cos(0,ir), '   s_cos0:',shape_s_cos(0,ir)
        write(io_neoout,20) '     cos1:',shape_cos(1,ir), '   s_cos1:',shape_s_cos(1,ir)
        write(io_neoout,20) '     cos2:',shape_cos(2,ir), '   s_cos2:',shape_s_cos(2,ir)
        write(io_neoout,20) '     cos3:',shape_cos(3,ir), '   s_cos3:',shape_s_cos(3,ir)
        write(io_neoout,20) '     cos4:',shape_cos(4,ir), '   s_cos4:',shape_s_cos(4,ir)
        write(io_neoout,20) '     cos5:',shape_cos(5,ir), '   s_cos5:',shape_s_cos(5,ir)
        write(io_neoout,20) '     cos6:',shape_cos(6,ir), '   s_cos6:',shape_s_cos(6,ir)
        
        write(io_neoout,*)
        write(io_neoout,20) ' dphi0/dr:',dphi0dr(ir)
        write(io_neoout,20) '    omega:',omega_rot(ir), &
             'domega/dr:',omega_rot_deriv(ir)
        write(io_neoout,20) '    epar0:',epar0(ir)

        write(io_neoout,*)
        write(io_neoout,'(a)') &
             'indx   z       n/n_norm     T/T_norm     m/m_norm     a/Ln         a/Lt         nu'
        do is=1,n_species
           write(io_neoout,'(t2,i2,2x,f7.3,2x,6(1pe11.4,2x))') &
                is,Z(is),dens(is,ir),temp(is,ir),mass(is),dlnndr(is,ir),dlntdr(is,ir),nu(is,ir)
        enddo
        if (ae_flag == 1) then
           write(io_neoout,'(t2,a3,1x,f7.3,2x,2(1pe11.4,2x),1x,a2,10x,2(1pe11.4,2x),1x,a2)') &
                'ade',-1.000,dens_ae(ir),temp_ae(ir),'--',dlnndr_ae(ir),dlntdr_ae(ir),'--'
        endif

        flag = 0
        do is=1,n_species
           if (aniso_model(is) >= 2) then
              flag = 1
           endif
        enddo
        if (flag == 1) then
           write(io_neoout,*)
           write(io_neoout,'(a)') 'Anisotropic Species'
           write(io_neoout,'(a)') &
                'indx   z       Tpar/T_norm  Tperp/T_norm a/Lt_par     a/Lt_perp'
           do is=1,n_species
              if (aniso_model(is) >= 2) then
                 write(io_neoout,'(t2,i2,2x,f7.3,2x,4(1pe11.4,2x))') &
                      is,Z(is),temp_para(is,ir),temp_perp(is,ir), &
                      dlntdr_para(is,ir),dlntdr_perp(is,ir)
              endif
           enddo
        endif

        write(io_neoout,*) '------------------'

     enddo
     write(io_neoout,*)

  endif

  if (silent_flag == 0 .and. i_proc == 0) then
     close(io_neoout)
  endif

10 format(t2,a,t19,': ',i3)
20 format(t2,2(a,1x,1pe11.4,4x))
30 format(t2,a,t32,': ',a)

end subroutine neo_check
