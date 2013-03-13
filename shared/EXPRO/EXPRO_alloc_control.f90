subroutine EXPRO_alloc_control(i_proc,flag)

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer, intent(in) :: i_proc
  integer, intent(in) :: flag

  integer :: ierr
  integer, parameter :: io=1


  if (flag == 1) then

     if (i_proc == 0 .and. EXPRO_n_exp == 0) then

        !--------------------------------------------------------------
        ! Read input.profiles.gen if it exists and EXPRO_n_exp not set
        ! 
        open(unit=io,&
             file=trim(path)//'input.profiles.gen',&
             status='old',&
             iostat=ierr)

        if (ierr == 0) then
           read(io,*) EXPRO_ncol
           read(io,*) EXPRO_nblock
           read(io,*) EXPRO_n_exp
           close(io)
           open(unit=io,file=trim(path)//trim(runfile),status='replace')
           write(io,'(a)') 'INFO: (EXPRO) input.profiles.gen found'
           close(io)
        else
           close(io)
           open(unit=io,file=trim(path)//trim(runfile),status='replace')
           write(io,'(a)') 'ERROR: (EXPRO) input.profiles.gen does not exist'
           close(io)
           stop
        endif
        !--------------------------------------------------------------

     endif

     if (i_proc == 0) then

        !--------------------------------------------------------------
        ! Read geometry Fourier harmonics
        ! 
        open(unit=io,&
             file=trim(path)//'input.profiles.geo',&
             status='old',&
             iostat=ierr)

        if (ierr /= 0) then
           EXPRO_nfourier = -1
           close(io)
        else
           call EXPRO_skip_header(io)
           read(io,*) EXPRO_nfourier
           close(io)
           open(unit=io,file=trim(path)//trim(runfile),status='replace')
           write(io,'(a)') 'INFO: (EXPRO) input.profiles.geo found'
           close(io)
        endif

        !--------------------------------------------------------------

     endif

     allocate(EXPRO_rho(EXPRO_n_exp))
     allocate(EXPRO_rmin(EXPRO_n_exp))
     allocate(EXPRO_rmaj(EXPRO_n_exp))
     allocate(EXPRO_q(EXPRO_n_exp))
     allocate(EXPRO_kappa(EXPRO_n_exp))

     allocate(EXPRO_delta(EXPRO_n_exp))
     allocate(EXPRO_te(EXPRO_n_exp))
     allocate(EXPRO_ne(EXPRO_n_exp))
     allocate(EXPRO_z_eff(EXPRO_n_exp))
     allocate(EXPRO_w0(EXPRO_n_exp))

     allocate(EXPRO_flow_mom(EXPRO_n_exp))
     allocate(EXPRO_pow_e(EXPRO_n_exp))
     allocate(EXPRO_pow_i(EXPRO_n_exp))
     allocate(EXPRO_pow_ei(EXPRO_n_exp))
     allocate(EXPRO_zeta(EXPRO_n_exp))

     allocate(EXPRO_flow_beam(EXPRO_n_exp))
     allocate(EXPRO_flow_wall(EXPRO_n_exp))
     allocate(EXPRO_zmag(EXPRO_n_exp))
     allocate(EXPRO_ptot(EXPRO_n_exp))
     allocate(EXPRO_dlnptotdr(EXPRO_n_exp))

     allocate(EXPRO_ni(nion_max,EXPRO_n_exp))
     allocate(EXPRO_ti(nion_max,EXPRO_n_exp))
     allocate(EXPRO_vtor(nion_max,EXPRO_n_exp))
     allocate(EXPRO_vpol(nion_max,EXPRO_n_exp))

     allocate(EXPRO_bunit(EXPRO_n_exp))
     allocate(EXPRO_s(EXPRO_n_exp))
     allocate(EXPRO_drmaj(EXPRO_n_exp))
     allocate(EXPRO_dzmag(EXPRO_n_exp))
     allocate(EXPRO_sdelta(EXPRO_n_exp))
     allocate(EXPRO_skappa(EXPRO_n_exp))
     allocate(EXPRO_szeta(EXPRO_n_exp))
     allocate(EXPRO_dlnnedr(EXPRO_n_exp))
     allocate(EXPRO_dlntedr(EXPRO_n_exp))

     allocate(EXPRO_dlnnidr(nion_max,EXPRO_n_exp))
     allocate(EXPRO_dlntidr(nion_max,EXPRO_n_exp))

     allocate(EXPRO_vol(EXPRO_n_exp))
     allocate(EXPRO_volp(EXPRO_n_exp))

     allocate(EXPRO_cs(EXPRO_n_exp))
     allocate(EXPRO_rhos(EXPRO_n_exp))

     allocate(EXPRO_w0p(EXPRO_n_exp))
     allocate(EXPRO_gamma_e(EXPRO_n_exp))
     allocate(EXPRO_gamma_p(EXPRO_n_exp))
     allocate(EXPRO_mach(EXPRO_n_exp))

     if (EXPRO_nfourier > 0) then  
        allocate(EXPRO_geo(4,0:EXPRO_nfourier,EXPRO_n_exp))
        allocate(EXPRO_dgeo(4,0:EXPRO_nfourier,EXPRO_n_exp))
     endif

     allocate(EXPRO_ni_new(EXPRO_n_exp))
     allocate(EXPRO_dlnnidr_new(EXPRO_n_exp))
     allocate(EXPRO_grad_r0(EXPRO_n_exp))
     allocate(EXPRO_ave_grad_r(EXPRO_n_exp))
     allocate(EXPRO_drdrho(EXPRO_n_exp))
     allocate(EXPRO_bp0(EXPRO_n_exp))
     allocate(EXPRO_bt0(EXPRO_n_exp))
     allocate(EXPRO_poloidalfluxover2pi(EXPRO_n_exp))

  else

     deallocate(EXPRO_rho)
     deallocate(EXPRO_rmin)
     deallocate(EXPRO_rmaj)
     deallocate(EXPRO_q)
     deallocate(EXPRO_kappa)

     deallocate(EXPRO_delta)
     deallocate(EXPRO_te)
     deallocate(EXPRO_ne)
     deallocate(EXPRO_z_eff)
     deallocate(EXPRO_w0)

     deallocate(EXPRO_flow_mom)
     deallocate(EXPRO_pow_e)
     deallocate(EXPRO_pow_i)
     deallocate(EXPRO_pow_ei)
     deallocate(EXPRO_zeta)

     deallocate(EXPRO_flow_beam)
     deallocate(EXPRO_flow_wall)
     deallocate(EXPRO_zmag)
     deallocate(EXPRO_ptot)
     deallocate(EXPRO_poloidalfluxover2pi)

     deallocate(EXPRO_ni)
     deallocate(EXPRO_ti)
     deallocate(EXPRO_vtor)
     deallocate(EXPRO_vpol)

     deallocate(EXPRO_bunit)
     deallocate(EXPRO_s)
     deallocate(EXPRO_drmaj)
     deallocate(EXPRO_dzmag)
     deallocate(EXPRO_sdelta)
     deallocate(EXPRO_skappa)
     deallocate(EXPRO_szeta)
     deallocate(EXPRO_dlnnedr)
     deallocate(EXPRO_dlntedr)

     deallocate(EXPRO_dlnnidr)
     deallocate(EXPRO_dlntidr)

     deallocate(EXPRO_vol)
     deallocate(EXPRO_volp)

     deallocate(EXPRO_cs)
     deallocate(EXPRO_rhos)

     deallocate(EXPRO_w0p)
     deallocate(EXPRO_gamma_e)
     deallocate(EXPRO_gamma_p)
     deallocate(EXPRO_mach)

     if (allocated(EXPRO_geo)) deallocate(EXPRO_geo)
     if (allocated(EXPRO_dgeo)) deallocate(EXPRO_dgeo)

     deallocate(EXPRO_ni_new)
     deallocate(EXPRO_dlnnidr_new)
     deallocate(EXPRO_grad_r0)
     deallocate(EXPRO_ave_grad_r)
     deallocate(EXPRO_drdrho)
     deallocate(EXPRO_bp0)
     deallocate(EXPRO_bt0)
     deallocate(EXPRO_dlnptotdr)

  endif

end subroutine EXPRO_alloc_control
