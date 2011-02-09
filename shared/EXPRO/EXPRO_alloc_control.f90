subroutine EXPRO_alloc_control(i_proc,path,flag)

  use EXPRO_interface

  implicit none

  character(len=*) :: path

  integer, intent(in) :: i_proc
  integer, intent(in) :: flag
  integer :: ierr
  integer, parameter :: io=1


  if (flag == 1) then

     if (i_proc == 0) then

        !--------------------------------------------------------------
        ! Read number of experimental gridpoints:
        ! 
        open(unit=io,&
             file=trim(path)//'input.profiles.gen',&
             status='old',&
             iostat=ierr)

        if (ierr /= 0) then
           print '(a)','ERROR: input.profiles.gen does not exist'
           stop
        endif

        read(io,*) EXPRO_ncol
        read(io,*) EXPRO_nblock
        read(io,*) EXPRO_n_exp

        close(io)
        !--------------------------------------------------------------

        !--------------------------------------------------------------
        ! Read number of experimental gridpoints:
        ! 
        open(unit=io,&
             file=trim(path)//'input.profiles.geo',&
             status='old',&
             iostat=ierr)

        if (ierr /= 0) then
           EXPRO_nfourier = -1
        else
           print '(a)','INFO: input.profiles.geo found'
           read(io,*) EXPRO_nfourier
        endif

        close(io)
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
     deallocate(EXPRO_bp0)
     deallocate(EXPRO_bt0)
     deallocate(EXPRO_poloidalfluxover2pi)

  endif

end subroutine EXPRO_alloc_control
