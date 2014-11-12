subroutine cgyro_make_profiles

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is,ir,ix 
  integer :: j
  integer :: num_ele
  integer, parameter :: io=20

  !---------------------------------------------------
  ! Manage electrons
  !
  num_ele = 0
  do is=1,n_species
     if (z(is) == -1) then
        num_ele = num_ele + 1
        is_ele = is
     endif
  enddo

  if (num_ele == 0) then

     ! Adiabatic electrons

     ae_flag = 1
     call cgyro_info('Using adiabatic electron model.')

     dens_ele = ne_ade
     temp_ele = te_ade

  else if (num_ele == 1) then

     ! GK electrons

     ae_flag = 0
     call cgyro_info('Using gyrokinetic electrons.')

     dens_ele = dens(is_ele)
     temp_ele = temp(is_ele)

  else

     call cgyro_error('Only one electron species allowed.')
     return

  endif
  !---------------------------------------------------

  ! Standard local simulation (one point)

  q = abs(q) 

  if (zf_test_flag == 1) then

     ! Zonal flow (n=0) test

     k_theta      = 0.0
     rho          = ky/(q/rmin)
     r_length_inv = s*ky/box_size

     call cgyro_info('Triggered zonal flow test.')

     if (n_radial /= 1) then
        call cgyro_error('For zonal flow test, set n_radial=1.')
        return
     endif

  else if (n_toroidal == 1) then

     ! Single linear mode (assume n=1, compute rho)

     k_theta      = q/rmin
     rho          = ky/k_theta
     r_length_inv = s*ky/box_size

     print *,rho

     call cgyro_info('Single-mode linear analysis.')

  else

     call cgyro_error('Nonlinear not implemented')
     return

  endif

  ! general geometry -- accessible only from interface 
  ! via parameters geo_ny_in and geo_yin_in
  geo_numeq_flag = 0
  geo_ny = 0
  allocate(geo_yin(8,0:geo_ny))
  geo_yin(:,:) = 0.0
  if (equilibrium_model == 3) then
     geo_numeq_flag = 1
     geo_ny = geo_ny_in  
     deallocate(geo_yin)
     allocate(geo_yin(8,0:geo_ny))
     do j=0,geo_ny
        geo_yin(:,j) = geo_yin_in(:,j)
     enddo
  endif

  ! Species-dependent quantities
  do is=1,n_species

     ! thermal velocity
     vth(is) = sqrt(temp(is)/mass(is))

     ! collision frequency
     nu(is) = nu_1_in *(1.0*z(is))**4/(1.0*z(1))**4 &
          * dens(is) / dens(1) &
          * sqrt(mass(1)/mass(is)) * (temp(1)/temp(is))**1.5

  enddo

  ! Fourier index mapping
  allocate(indx_xi(n_xi))
  do ix=1,n_xi
     indx_xi(ix) = ix-1
  enddo
  allocate(indx_r(n_radial))
  do ir=1,n_radial
     indx_r(ir) = -n_radial/2 + (ir-1)
  enddo
  if (zf_test_flag == 1) indx_r(1) = 1

  ! Print the re-mapped equilibrium data
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//'out.cgyro.equil',status='replace')
     write (io,'(e16.8)',advance='no') rmin
     write (io,'(e16.8)',advance='no') rmaj
     write (io,'(e16.8)',advance='no') q
     write (io,'(e16.8)',advance='no') s
     write (io,'(e16.8)',advance='no') rho
     write (io,'(e16.8)',advance='no') ky
     do is=1,n_species
        write (io,'(e16.8)',advance='no') dens(is)
        write (io,'(e16.8)',advance='no') temp(is)
        write (io,'(e16.8)',advance='no') dlnndr(is)
        write (io,'(e16.8)',advance='no') dlntdr(is)
        write (io,'(e16.8)',advance='no') nu(is)
     enddo
     write (io,*)
     close(io)
  endif

end subroutine cgyro_make_profiles
