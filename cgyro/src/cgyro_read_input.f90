subroutine cgyro_read_input

  use mpi
  use cgyro_globals

  implicit none

  integer :: is
  integer :: xint
  real :: x
  character (len=1) :: cdummy

  if (i_proc == 0) open(unit=1,file=trim(path)//'input.cgyro.gen',status='old')

  call readbc_int(n_energy)
  call readbc_int(n_xi)
  call readbc_int(n_theta)
  call readbc_int(n_radial)
  call readbc_int(n_toroidal)
  call readbc_int(n_field)
  call readbc_real(e_max)
  call readbc_int(e_method)
  call readbc_real(delta_t)
  call readbc_real(max_time)
  call readbc_int(print_step)
  call readbc_int(restart_step)
  call readbc_real(freq_tol)
  call readbc_int(restart_write)
  call readbc_int(restart_mode)
  call readbc_real(up_radial)
  call readbc_real(up_theta)
  call readbc_int(nup_theta)
  call readbc_int(nup_radial)
  call readbc_int(implicit_flag)
  call readbc_int(constant_wind_flag)
  call readbc_real(ky)
  call readbc_int(box_size)
  call readbc_real(ipccw)
  call readbc_real(btccw)
  call readbc_int(silent_flag)
  call readbc_int(profile_model)
  call readbc_int(equilibrium_model)
  call readbc_int(collision_model)
  call readbc_int(collision_mom_restore)
  call readbc_int(collision_ene_restore)
  call readbc_int(collision_ene_diffusion)
  call readbc_int(collision_kperp)
  call readbc_int(collision_field_model)
  call readbc_int(zf_test_flag)
  call readbc_int(nonlinear_flag)
  call readbc_int(nonlinear_method)
  call readbc_real(te_ade)
  call readbc_real(ne_ade)
  call readbc_real(masse_ade)
  call readbc_real(lambda_debye)
  call readbc_int(test_flag)
  call readbc_int(h_print_flag)
  call readbc_real(amp0)
  call readbc_real(amp)
  call readbc_real(gamma_e)
  call readbc_real(gamma_p)
  call readbc_real(mach)
  call readbc_real(error_tol)
  call readbc_int(kxfilter_flag)
  call readbc_real(gamma_e_decay)
  call readbc_int(hiprec_flag)
  call readbc_int(udsymmetry_flag)

  call readbc_real(rmin)
  call readbc_real(rmaj)
  call readbc_real(q)
  call readbc_real(s)
  call readbc_real(shift)    
  call readbc_real(kappa)   
  call readbc_real(s_kappa)  
  call readbc_real(delta)       
  call readbc_real(s_delta)
  call readbc_real(zeta)      
  call readbc_real(s_zeta)
  call readbc_real(zmag)       
  call readbc_real(dzmag)
  call readbc_real(beta_star)
  call readbc_real(betae_unit)

  call readbc_int(subroutine_flag)

  call readbc_int(n_species)

  call readbc_real(nu_ee)

  do is=1,6
     call readbc_int(xint) ; z(is) = xint
     call readbc_real(x)   ; mass(is) = x
     call readbc_real(x)   ; dens(is) = x
     call readbc_real(x)   ; temp(is) = x
     call readbc_real(x)   ; dlnndr(is) = x
     call readbc_real(x)   ; dlntdr(is) = x
  enddo

  call readbc_real(lambda_debye_scale)
  call readbc_real(gamma_e_scale)
  call readbc_real(gamma_p_scale)
  call readbc_real(mach_scale)
  call readbc_real(q_scale)
  call readbc_real(s_scale)
  do is=1,6
     call readbc_real(x)   ; dlnndr_scale(is) = x
     call readbc_real(x)   ; dlntdr_scale(is) = x
  enddo

  if (i_proc == 0) close(1)

  ! GEO fourier coefficients
  geo_ny_in = 0
  geo_yin_in(:,:) = 0.0
  if (subroutine_flag == 0 .and. equilibrium_model == 3 & 
       .and. profile_model == 1) then
     if (i_proc == 0) then
        open(unit=1,file=trim(path)//'input.geo',status='old')
        ! header skip
        do
           read(1,'(a)') cdummy
           if (cdummy /= '#') exit
        enddo
        backspace 1
        ! n_fourier
        read(1,*) geo_ny_in
        ! fourier coefficients
        read(1,*) geo_yin_in(:,0:geo_ny_in)
        close(1)
     endif
     call MPI_BCAST(geo_ny_in,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)
     call MPI_BCAST(geo_yin_in,size(geo_yin_in),MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)
  endif

end subroutine cgyro_read_input

!------------------------------------------------------------
! Service routines: 
!
! (1) read and broadcast an integer:
!
subroutine readbc_int(p)

  use mpi
  use cgyro_globals

  implicit none
  integer, intent(inout) :: p

  if (i_proc == 0) read(1,*) p

  call MPI_BCAST(p,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

end subroutine readbc_int
!
! (2) read and broadcast a real:
!
subroutine readbc_real(x)
  
  use mpi
  use cgyro_globals

  implicit none
  real, intent(inout) :: x

  if (i_proc == 0) read(1,*) x

  call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)

end subroutine readbc_real
!------------------------------------------------------------
