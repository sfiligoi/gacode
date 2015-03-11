!------------------------------------------------------------
! tgyro_stab_driver.f90
!
! PURPOSE:
!  Manage linear stability calculation for GYRO/TGLF.
!------------------------------------------------------------

subroutine tgyro_stab_driver

  use mpi
  use tgyro_globals
  use gyro_interface
  use tglf_interface
  use glf23_interface

  implicit none

  integer :: iky
  integer :: i_err
  integer :: i
  integer :: imax(1)
  real :: wr,wi,err

  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wi_ion
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wi_elec
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wr_ion
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wr_elec
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wi_ion_gather
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wi_elec_gather
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wr_ion_gather
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wr_elec_gather

  real :: wi_ion_loc 
  real :: wr_ion_loc 
  real :: wi_elec_loc 
  real :: wr_elec_loc

  real :: wi_ion_vec(n_worker) 
  real :: wr_ion_vec(n_worker) 
  real :: wi_elec_vec(n_worker) 
  real :: wr_elec_vec(n_worker) 

  real, dimension(0:tgyro_stab_nky-1) :: ky

  integer, parameter :: i_print = 0

  wi_elec(:,:) = 0.0
  wr_elec(:,:) = 0.0
  wi_ion(:,:)  = 0.0
  wr_ion(:,:)  = 0.0
  wi_elec_gather(:,:) = 0.0
  wr_elec_gather(:,:) = 0.0
  wi_ion_gather(:,:)  = 0.0
  wr_ion_gather(:,:)  = 0.0

  call tgyro_profile_functions

  select case (flux_method)

  case (2)

     ! Map TGYRO parameters to TGLF
     call tgyro_tglf_map

     ! Set use to linear stability, not transport
     tglf_use_transport_model_in = .false.
     tglf_nmodes_in = 2

     do iky=0,tgyro_stab_nky-1

        wi_ion_loc  = 0.0
        wi_elec_loc = 0.0
        wr_ion_loc  = 0.0
        wr_elec_loc = 0.0

        ky(iky) = tgyro_stab_kymin + iky*tgyro_stab_deltaky

        tglf_ky_in = ky(iky)

        call tglf_run_mpi()

        ! wr and wi are now the COMPUTED eigenvalues

        do i=1,tglf_nmodes_in

           wr = real(tglf_eigenvalue_out(i))
           wi = aimag(tglf_eigenvalue_out(i))

           ! electron mode
           if (wr > 0.0 .and. wi > wi_elec_loc) then
              wi_elec_loc = wi
              wr_elec_loc = wr
           endif

           ! ion mode
           if (wr < 0.0 .and. wi > wi_ion_loc) then
              wi_ion_loc = wi
              wr_ion_loc = wr
           endif

        enddo

        wi_elec(iky,i_r-1) = wi_elec_loc
        wr_elec(iky,i_r-1) = wr_elec_loc

        wi_ion(iky,i_r-1) = wi_ion_loc     
        wr_ion(iky,i_r-1) = wr_ion_loc

        if (worker == 0 .and. i_print == 1) then
           print '(2(f4.2,1x),2(1pe13.5,1x))', &
                r(i_r)/r_min,ky(iky),wr_elec(iky,i_r-1),wi_elec(iky,i_r-1)
           print '(2(f4.2,1x),2(1pe13.5,1x))', &
                r(i_r)/r_min,ky(iky),wr_ion(iky,i_r-1),wi_ion(iky,i_r-1)
        endif

     enddo

   case (3)
     ! Map TGYRO parameters to GLF23
     call tgyro_glf23_map
     ! Set use to linear stability, not transport
     glf23_use_transport_model_in = .false.

     glf23_nmodes_in = 2

     do iky=0,tgyro_stab_nky-1

        wi_ion_loc  = 0.0
        wi_elec_loc = 0.0
        wr_ion_loc  = 0.0
        wr_elec_loc = 0.0

        ky(iky) = tgyro_stab_kymin + iky*tgyro_stab_deltaky

        glf23_ky_in = ky(iky)
 
        call glf23_run()

        ! wr and wi are now the COMPUTED eigenvalues

        do i=1,glf23_nmodes_in

           wr = real(tglf_eigenvalue_out(i))
           wi = aimag(tglf_eigenvalue_out(i))

           ! electron mode
           if (wr > 0.0 .and. wi > wi_elec_loc) then
              wi_elec_loc = wi
              wr_elec_loc = wr
           endif

           ! ion mode
           if (wr < 0.0 .and. wi > wi_ion_loc) then
              wi_ion_loc = wi
              wr_ion_loc = wr
           endif

        enddo

        wi_elec(iky,i_r-1) = wi_elec_loc
        wr_elec(iky,i_r-1) = wr_elec_loc

        wi_ion(iky,i_r-1) = wi_ion_loc     
        wr_ion(iky,i_r-1) = wr_ion_loc

        if (worker == 0 .and. i_print == 1) then
           print '(2(f4.2,1x),2(1pe13.5,1x))', &
                r(i_r)/r_min,ky(iky),wr_elec(iky,i_r-1),wi_elec(iky,i_r-1)
           print '(2(f4.2,1x),2(1pe13.5,1x))', &
                r(i_r)/r_min,ky(iky),wr_ion(iky,i_r-1),wi_ion(iky,i_r-1)
        endif

     enddo

  case (4)

     ! Map TGYRO parameters to GYRO
     call tgyro_gyro_map

     ! Use FIELDEIGEN
     gyro_linsolve_method_in = 3

     ! Use GK electron (only possibility with FIELDEIGEN)
     gyro_electron_method_in = 4

     ! Set this so was can specify L_Y=ky*rhos exactly.
     gyro_box_multiplier_in = -1.0

     ! Silent running
     gyro_silent_flag_in = 1

     do iky=0,tgyro_stab_nky-1

        wi_ion_loc  = 0.0
        wi_elec_loc = 0.0
        wr_ion_loc  = 0.0
        wr_elec_loc = 0.0

        ! Set wavenumber
        ky(iky) = tgyro_stab_kymin + iky*tgyro_stab_deltaky
        ! Pass to GYRO
        gyro_l_y_in = ky(iky)

        ! Set wr for initial guess
        if (n_worker == 1) then
           wr = 0.0
        else
           wr = ky(iky)*(-1.0+worker*(2.0/(n_worker-1)))
        endif

        ! Pass initial guess to GYRO
        ! NOTE: gyro_fieldeigen_wi_in taken from input.gyro (should be about 0.1)
        gyro_fieldeigen_wr_in = wr

        call gyro_run(gyrotest_flag,gyro_restart_method,transport_method)

        ! wr and wi are now the COMPUTED eigenvalues

        wr  = real(gyro_fieldeigen_omega_out)
        wi  = aimag(gyro_fieldeigen_omega_out)
        err = gyro_fieldeigen_error_out

        ! electron mode
        if (wr > 0.0 .and. wi > 0.0 .and. err < 0.9) then
           wi_elec_loc = wi
           wr_elec_loc = wr
        endif

        ! ion mode
        if (wr < 0.0 .and. wi > 0.0 .and. err < 0.9) then
           wi_ion_loc = wi
           wr_ion_loc = wr
        endif

        call MPI_ALLGATHER(wi_elec_loc,1,MPI_DOUBLE_PRECISION,&
             wi_elec_vec,1,MPI_DOUBLE_PRECISION,gyro_rad,i_err)
        call MPI_ALLGATHER(wr_elec_loc,1,MPI_DOUBLE_PRECISION,&
             wr_elec_vec,1,MPI_DOUBLE_PRECISION,gyro_rad,i_err)
        call MPI_ALLGATHER(wi_ion_loc,1,MPI_DOUBLE_PRECISION,&
             wi_ion_vec,1,MPI_DOUBLE_PRECISION,gyro_rad,i_err)
        call MPI_ALLGATHER(wr_ion_loc,1,MPI_DOUBLE_PRECISION,&
             wr_ion_vec,1,MPI_DOUBLE_PRECISION,gyro_rad,i_err)

        imax = maxloc(wi_elec_vec)
        wi_elec(iky,i_r-1) = wi_elec_vec(imax(1))     
        wr_elec(iky,i_r-1) = wr_elec_vec(imax(1))     

        imax = maxloc(wi_ion_vec)
        wi_ion(iky,i_r-1) = wi_ion_vec(imax(1))     
        wr_ion(iky,i_r-1) = wr_ion_vec(imax(1))     

        if (worker == 0 .and. i_print == 1) then
           print '(2(f4.2,1x),2(1pe13.5,1x))', &
                r(i_r)/r_min,ky(iky),wr_elec(iky,i_r-1),wi_elec(iky,i_r-1)
           print '(2(f4.2,1x),2(1pe13.5,1x))', &
                r(i_r)/r_min,ky(iky),wr_ion(iky,i_r-1),wi_ion(iky,i_r-1)
        endif

     enddo

  case default

     call tgyro_catch_error('ERROR: Bad value in tgyro_stab_driver.')

  end select

  call MPI_ALLGATHER(wi_elec(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wi_elec_gather,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  call MPI_ALLGATHER(wi_ion(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wi_ion_gather,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  call MPI_ALLGATHER(wr_elec(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wr_elec_gather,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  call MPI_ALLGATHER(wr_ion(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wr_ion_gather,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  if (i_proc_global == 0) then

     if (i_tran == 0) then
        open(unit=1,file='out.tgyro.wi_elec',status='replace')
     else
        open(unit=1,file='out.tgyro.wi_elec',position='append')
     endif

     write(1,'(t11,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wi_elec_gather(:,i)
     enddo

     close(1)

     if (i_tran == 0) then
        open(unit=1,file='out.tgyro.wi_ion',status='replace')
     else
        open(unit=1,file='out.tgyro.wi_ion',position='append')
     endif

     write(1,'(t11,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wi_ion_gather(:,i)
     enddo

     close(1)

     if (i_tran == 0) then
        open(unit=1,file='out.tgyro.wr_elec',status='replace')
     else
        open(unit=1,file='out.tgyro.wr_elec',position='append')
     endif

     write(1,'(t11,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wr_elec_gather(:,i)
     enddo

     close(1)

     if (i_tran == 0) then
        open(unit=1,file='out.tgyro.wr_ion',status='replace')
     else
        open(unit=1,file='out.tgyro.wr_ion',position='append')
     endif

     write(1,'(t11,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wr_ion_gather(:,i)
     enddo

     close(1)

  endif

10 format(f5.3,3x,10(1pe12.5,1x))
20 format(t2,a,5x,10(f7.3,6x))

end subroutine tgyro_stab_driver
