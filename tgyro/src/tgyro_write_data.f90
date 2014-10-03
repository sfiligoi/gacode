!-----------------------------------------------------------
! tgyro_write_data.f90
!
! PURPOSE:
!  Manage I/O for everything but iteration diagnostics.
!
! NOTES:
!  Files are opened and closed at each iteration.
!----------------------------------------------------------

subroutine tgyro_write_data(i_print)

  use mpi
  use tgyro_globals

  implicit none

  integer :: i
  integer :: ip
  integer :: i_print
  integer :: i_ion
  integer :: p
  integer, parameter :: trinity_flag=0
  real, dimension(2:n_r,4) :: res2,relax2

  !--------------------------------------------------------------------------------
  ! First, generate and write TGLF linear growth rates
  !
  call tgyro_stab_driver
  !--------------------------------------------------------------------------------

  if (i_proc_global > 0) return

  if (i_print == 0) then

     open(unit=1,file='out.tgyro.geometry.1',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.geometry.2',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.nu_rho',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.power_e',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.power_i',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.alpha',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.flux_target',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.mflux_target',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.gyrobohm',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.gradient',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.residual',status='replace')
     close(1)

     do i_ion=1,loc_n_ion
        open(unit=1,file='out.tgyro.flux_i'//trim(ion_tag(i_ion)),status='replace')
        close(1)
        open(unit=1,file='out.tgyro.mflux_i'//trim(ion_tag(i_ion)),status='replace')
        close(1)
        open(unit=1,file='out.tgyro.profile'//trim(ion_tag(i_ion)),status='replace')
        close(1)
     enddo

     open(unit=1,file='out.tgyro.flux_e',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.mflux_e',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.profile',status='replace')
     close(1)

     if (trinity_flag == 1) then
        open(unit=1,file='out.tgyro.trinity.eflux.out',status='replace')
        close(1)
     endif

     open(unit=1,file='out.tgyro.control',status='replace')
     write(1,*) n_r
     write(1,*) n_evolve
     write(1,*) 0
     close(1)

     ! Special case: iteration status file (see loc_write_intermediate)
     open(unit=1,file='out.tgyro.iterate',status='replace')
     close(1)

     open(unit=1,file='out.tgyro.prec',status='replace')
     close(1)

     return

  endif

  ! Geometry 1 [** constant **]

  open(unit=1,file='out.tgyro.geometry.1',status='old',position='append')

  write(1,20) 'r/a','rho','q','s','kappa','s_kappa','delta','s_delta','shift','rmaj/a','b_unit'
  write(1,20) '','','','','','','','','','','(T)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          rho(i),&
          q(i),&
          s(i),&
          kappa(i),&
          s_kappa(i),&
          delta(i),&
          s_delta(i),&
          shift(i),&
          r_maj(i)/r_min,&
          b_unit(i)/1e4
  enddo

  close(1)

  ! Geometry 2 [** constant **]

  open(unit=1,file='out.tgyro.geometry.2',status='old',position='append')

  write(1,20) 'r/a','zmag/a','dzmag','zeta','s_zeta','volume','d(vol)/dr','<|grad_r|>','rmin'
  write(1,20) '','','','','','(m^3)','(m^2)','','(cm)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          zmag(i)/r_min,&
          dzmag(i),&
          zeta(i),&
          s_zeta(i),&
          vol(i)/1e6,&
          volp(i)/1e4,&
          ave_grad_r(i),&
          r(i)

  enddo

  close(1)

  ! Collisions and gyroradii (out.tgyro.nu_rho)

  open(unit=1,file='out.tgyro.nu_rho',status='old',position='append')

  write(1,20) 'r/a','(a/cs)/t_ii','(a/cs)/t_ee','1/nue_star','(a/cs)nu_exch','rho_i/a','rho_s/a','frac_ae'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          nui(1,i)*r_min/c_s(i),&
          nue(i)*r_min/c_s(i),&
          nue_star(i),&
          nu_exch(i)*r_min/c_s(i),&
          rho_i(i)/r_min,&
          rho_s(i)/r_min,&
          frac_ae(i)
  enddo

  close(1)

  ! Electron powers (erg/s = 1e-7 W = 1e-7*1e-6 MW)

  open(unit=1,file='out.tgyro.power_e',status='old',position='append')

  write(1,20) 'r/a','p_e_fus','p_e_aux','p_brem','p_sync','p_line','p_exch[-]','p_expwd[-]','p_e_tot'
  write(1,20) '','(MW)','(MW)','(MW)','(MW)','(MW)','(MW)','(MW)','(MW)'
  do i=1,n_r
     write(1,10) &
          r(i)/r_min,&
          +p_e_fus(i)*1e-7*1e-6,&
          +p_e_aux_in(i)*1e-7*1e-6,&
          -p_brem(i)*1e-7*1e-6,&
          -p_sync(i)*1e-7*1e-6,&
          -p_line_in(i)*1e-7*1e-6, &
          -p_exch(i)*1e-7*1e-6,&
          -p_expwd(i)*1e-7*1e-6,&
          +p_e(i)*1e-7*1e-6
  enddo

  close(1) 

  ! Ion powers (erg/s = 1e-7 W = 1e-7*1e-6 MW)

  open(unit=1,file='out.tgyro.power_i',status='old',position='append')

  write(1,20) 'r/a','p_i_fus','p_i_aux','p_exch','p_expwd','p_i_tot'
  write(1,20) '','(MW)','(MW)','(MW)','(MW)','(MW)'
  do i=1,n_r
     write(1,10) &
          r(i)/r_min,&
          p_i_fus(i)*1e-7*1e-6,&
          p_i_aux_in(i)*1e-7*1e-6,&
          p_exch(i)*1e-7*1e-6,&
          p_expwd(i)*1e-7*1e-6,&
          p_i(i)*1e-7*1e-6 
  enddo

  close(1) 

  ! Alpha parameters (out.tgyro.alpha) (erg/s = 1e-7 W = 1e-7*1e-6 MW)

  open(unit=1,file='out.tgyro.alpha',status='old',position='append')

  write(1,20) 'r/a','s_alpha','s_alpha_i','s_alpha_e','frac_ai','frac_ae','E_alpha/E_c'
  write(1,20) '','(MW/cm^3)','(MW/cm^3)','(MW/cm^3)','','',''
  do i=1,n_r
     write(1,10) &
          r(i)/r_min,&
          (s_alpha_i(i)+s_alpha_e(i))*1e-7*1e-6,&
          s_alpha_i(i)*1e-7*1e-6,&
          s_alpha_e(i)*1e-7*1e-6,&
          frac_ai(i),&
          frac_ae(i),&
          e_alpha/e_cross(i)
  enddo

  close(1) 

  ! Transport+target fluxes (flux_target.out)

  open(unit=1,file='out.tgyro.flux_target',status='old',position='append')

  write(1,20) 'r/a','eflux_i_tot','eflux_i_target','eflux_e_tot','eflux_e_target','pflux_e_tot','pflux_e_target'
  write(1,20) '','(GB)','(GB)','(GB)','(GB)','(GB)','(GB)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          eflux_i_tot(i),&
          eflux_i_target(i),&
          eflux_e_tot(i),&
          eflux_e_target(i),&
          pflux_e_tot(i),&
          pflux_e_target(i)
  enddo

  close(1)

  ! Transport+target fluxes for momentum (mflux_target.out)

  open(unit=1,file='out.tgyro.mflux_target',status='old',position='append')

  write(1,20) 'r/a','mflux_tot','mflux_target'
  write(1,20) '','(GB)','(GB)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          mflux_tot(i),&
          mflux_target(i)
  enddo

  close(1)

  ! Ion particle and energy fluxes (flux_i.out)

  open(unit=1,file='out.tgyro.flux_i',status='old',position='append')

  write(1,20) 'r/a','pflux_i_neo','pflux_i_tur','eflux_i_neo','eflux_i_tur'
  write(1,20) '','(GB)','(GB)','(GB)','(GB)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          pflux_i_neo(1,i),&
          pflux_i_tur(1,i),&
          eflux_i_neo(1,i),&
          eflux_i_tur(1,i)
  enddo

  close(1)

  ! Ion momentum fluxes and exchange powers (mflux_i.out)

  open(unit=1,file='out.tgyro.mflux_i',status='old',position='append')

  write(1,20) 'r/a','mflux_i_neo','mflux_i_tur','expwd_i_tur'
  write(1,20) '','(GB)','(GB)','(GB)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          mflux_i_neo(1,i),&
          mflux_i_tur(1,i),&
          expwd_i_tur(1,i)
  enddo

  close(1)

  ! Electron particle and energy fluxes (flux_e.out)

  open(unit=1,file='out.tgyro.flux_e',status='old',position='append')

  write(1,20) 'r/a','pflux_e_neo','pflux_e_tur','eflux_e_neo','eflux_e_tur'
  write(1,20) '','(GB)','(GB)','(GB)','(GB)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          pflux_e_neo(i),&
          pflux_e_tur(i),&
          eflux_e_neo(i),&
          eflux_e_tur(i)
  enddo

  close(1)

  ! Electron momentum fluxes and exchange powers (mflux_e.out)

  open(unit=1,file='out.tgyro.mflux_e',status='old',position='append')

  write(1,20) 'r/a','mflux_e_neo','mflux_e_tur','expwd_e_tur'
  write(1,20) '','(GB)','(GB)','(GB)'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          mflux_e_neo(i),&
          mflux_e_tur(i),&
          expwd_e_tur(i)
  enddo

  close(1)

  ! gyroBohm factors in physical units

  open(unit=1,file='out.tgyro.gyrobohm',status='old',position='append')

  write(1,20) 'r/a','Chi_GB','Q_GB','Gamma_GB','Pi_GB','c_s'
  write(1,20) '','m^2/s','MW/m^2','10^19/m^2/s','J/m^2','m/s'
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          chi_gb(i)*1e-4,&
          q_gb(i)*1e-7*1e-6/1e-4,&
          gamma_gb(i)*1e-19/1e-4,&
          pi_gb(i)*1e-7/1e-4,&
          c_s(i)/100.0
  enddo

  close(1)

  ! Temperature and density profiles

  open(unit=1,file='out.tgyro.profile',status='old',position='append')

  write(1,20) 'r/a','ni','ne','ti','te','ti/te','betae_unit','M=wR/cs'
  write(1,20) '','(1/cm^3)','(1/cm^3)','(keV)','(keV)','','',''
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          ni(1,i),&
          ne(i),&
          ti(1,i)/1e3,&
          te(i)/1e3,&
          ti(1,i)/te(i),&
          betae_unit(i),&
          u00(i)/c_s(i)
  enddo

  close(1)

  ! Temperature and density gradient profiles

  open(unit=1,file='out.tgyro.gradient',status='old',position='append')

  write(1,20) 'r/a','a/Lni','a/Lne','a/LTi','a/LTe','a/Lp','a*gamma_e/cs','a*gamma_p/cs'
  write(1,20) '','','','','',''
  do i=1,n_r
     write(1,10) r(i)/r_min,&
          r_min*dlnnidr(1,i),&
          r_min*dlnnedr(i),&
          r_min*dlntidr(1,i),&
          r_min*dlntedr(i),&
          r_min*dlnpdr(i),&
          r_min/c_s(i)*gamma_eb(i),&
          r_min/c_s(i)*gamma_p(i)
  enddo

  close(1)

  do i_ion=2,loc_n_ion

     ! Ion 2 particle and energy fluxes (flux_i*.out)

     open(unit=1,&
          file='out.tgyro.flux_i'//trim(ion_tag(i_ion)),&
          status='old',position='append')

     write(1,20) 'r/a','pflux_i_neo','pflux_i_tur','eflux_i_neo','eflux_i_tur'
     write(1,20) '','(GB)','(GB)','(GB)','(GB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             pflux_i_neo(i_ion,i),&
             pflux_i_tur(i_ion,i),&
             eflux_i_neo(i_ion,i),&
             eflux_i_tur(i_ion,i)
     enddo

     close(1)

     ! Ion 2 momentum fluxes and exchange powers (mflux_i*.out)

     open(unit=1,&
          file='out.tgyro.mflux_i'//trim(ion_tag(i_ion)),&
          status='old',position='append')

     write(1,20) 'r/a','mflux_i_neo','mflux_i_tur','expwd_i_tur'
     write(1,20) '','(GB)','(GB)','(GB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             mflux_i_neo(i_ion,i),&
             mflux_i_tur(i_ion,i),&
             expwd_i_tur(i_ion,i)
     enddo

     close(1)

     ! Impurity profiles

     open(unit=1,&
          file='out.tgyro.profile'//trim(ion_tag(i_ion)),&
          status='old',position='append')

     write(1,20) 'r/a','ni','a/Lni','Ti','a/LTi'
     write(1,20) '','(1/cm^3)','','(keV)',''
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             ni(i_ion,i),&
             r_min*dlnnidr(i_ion,i),&
             ti(i_ion,i)/1e3,&
             r_min*dlntidr(i_ion,i)
     enddo

     close(1)

  enddo ! i_ion

  ! Residuals

  open(unit=1,file='out.tgyro.residual',status='old',position='append')

  if (tgyro_relax_iterations == 0) then
     write(1,30) 'ITERATION*: ',i_tran,sum(res)/size(res),flux_counter*n_worker*n_inst
  else 
     write(1,30) 'ITERATION : ',i_tran,sum(res)/size(res),flux_counter*n_worker*n_inst
  endif

  res2(:,:)   = 0.0
  relax2(:,:) = 0.0 

  p = 0
  do i=2,n_r
     if (loc_ti_feedback_flag == 1) then
        p  = p+1
        res2(i,1) = res(p)
        relax2(i,1) = relax(p)
     endif
     if (loc_te_feedback_flag == 1) then
        p  = p+1
        res2(i,2) = res(p)
        relax2(i,2) = relax(p)
     endif
     if (loc_ne_feedback_flag == 1) then
        p  = p+1
        res2(i,3) = res(p)
        relax2(i,3) = relax(p)
     endif
     if (loc_er_feedback_flag == 1) then
        p  = p+1
        res2(i,4) = res(p)
        relax2(i,4) = relax(p)
     endif
     write(1,40) &
          r(i)/r_min,(res2(i,ip),relax2(i,ip),ip=1,4)
  enddo

  close(1)

  ! Control (control.out)

  open(unit=1,file='out.tgyro.control',status='old',position='append')
  backspace(1)
  write(1,*) i_tran
  close(1)

  open(unit=1,file='out.tgyro.prec',status='old',position='append')
  write(1,*) sum(abs(eflux_i_tot(:))+abs(eflux_e_tot(:)))
  close(1)

  !--------------------------------------------------------------------------------
  ! Trinity-type fluxes
  ! Electron particle and energy fluxes (flux_e.out)

  if (trinity_flag == 1) then
     open(unit=1,file='out.tgyro.trinity.eflux.out',status='old',position='append')

     write(1,20) 'r/a','eflux_i_neo','eflux_e_neo','eflux_i_tur','eflux_e_tur'
     write(1,20) '','(TGB)','(TGB)','(TGB)','(TGB)'
     do i=1,n_r
        write(1,10) r(i)/r_min,&
             eflux_i_neo(1,i)*q_gb(i)/q_tgb(i),&
             eflux_e_neo(i)*q_gb(i)/q_tgb(i),&
             eflux_i_tur(1,i)*q_gb(i)/q_tgb(i),&
             eflux_e_tur(i)*q_gb(i)/q_tgb(i)
     enddo
     close(1)
  endif
  !--------------------------------------------------------------------------------

  ! Write progress to screen
  if (i_tran < 10) then
     print '(a,i1,a,1pe10.3,a)', 'INFO: (TGYRO) Finished iteration ',i_tran,' [',sum(res)/size(res),']'
  else if (i_tran < 100) then
     print '(a,i2,a,1pe10.3,a)', 'INFO: (TGYRO) Finished iteration ',i_tran,' [',sum(res)/size(res),']'
  else
     print '(a,i3,a,1pe10.3,a)', 'INFO: (TGYRO) Finished iteration ',i_tran,' [',sum(res)/size(res),']'
  endif

  ! Data
10 format(t1,11(1pe13.6,2x))
  ! Text headers
20 format(t2,a,t17,a,t32,a,t47,a,t62,a,t77,a,t92,a,t107,a,t122,a,t137,a,t152,a)
  ! Residual header
30 format(t2,a,i3,1pe12.5,2x,'[',i6,']')
  ! Residuals
40 format(t2,f8.6,4(2x,2(1pe10.3,1x)))

end subroutine tgyro_write_data
