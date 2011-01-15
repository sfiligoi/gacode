!------------------------------------------------------------
! tgyro_stab_driver.f90
!
! PURPOSE:
!  Manage linear stability calculation
!------------------------------------------------------------

subroutine tgyro_stab_driver

  use tgyro_globals
  use gyro_interface

  implicit none

  integer :: i_search
  integer :: iky
  integer :: i_err
  integer :: i
  integer :: imax(1)
  real :: wr,wi,err

  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wi_ion
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wi_elec
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wr_ion
  real, dimension(0:tgyro_stab_nky-1,n_inst) :: wr_elec

  real :: wi_ion_loc 
  real :: wr_ion_loc 
  real :: wi_elec_loc 
  real :: wr_elec_loc

  real :: wi_ion_vec(n_worker) 
  real :: wr_ion_vec(n_worker) 
  real :: wi_elec_vec(n_worker) 
  real :: wr_elec_vec(n_worker) 

  character (len=7) :: ktag
  real, dimension(0:tgyro_stab_nky-1) :: ky

  include 'mpif.h'

  ! Map TGYRO parameters to GYRO
  call tgyro_gyro_map

  ! Use FIELDEIGEN
  gyro_linsolve_method_in = 3

  ! Use GK electron (only possibility with FIELDEIGEN)
  gyro_electron_method_in = 4

  ! Set this so was can specify L_Y=ky*rhos exactly.
  gyro_box_multiplier_in = -1.0

  wi_elec(:,:) = 0.0
  wr_elec(:,:) = 0.0
  wi_ion(:,:)  = 0.0
  wr_ion(:,:)  = 0.0

  do iky=0,tgyro_stab_nky-1

     wi_ion_loc  = 0.0
     wi_elec_loc = 0.0
     wr_ion_loc  = 0.0
     wr_elec_loc = 0.0

     ky(iky) = tgyro_stab_kymin + iky*0.1

     if (n_worker == 1) then
        wr = 0.0
     else
        wr = ky(iky)*(-1.0+worker*(2.0/(n_worker-1)))
     endif

     gyro_l_y_in = ky(iky)
     gyro_fieldeigen_wr_in = wr
     
     ! NOTE: gyro_fieldeigen_wi_in taken from input.gyro (should be about 0.1)

     call gyro_run(gyrotest_flag, gyro_restart_method, &
          transport_method, gyro_exit_status(i_r), gyro_exit_message(i_r))

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

     if (worker == 0) then
        print '(2(f4.2,1x),2(1pe13.5,1x))', &
             r(i_r)/r_min,ky(iky),wr_elec(iky,i_r-1),wi_elec(iky,i_r-1)
        print '(2(f4.2,1x),2(1pe13.5,1x))', &
             r(i_r)/r_min,ky(iky),wr_ion(iky,i_r-1),wi_ion(iky,i_r-1)
     endif

  enddo

  call MPI_ALLGATHER(wi_elec(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wi_elec,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  call MPI_ALLGATHER(wi_ion(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wi_ion,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  call MPI_ALLGATHER(wr_elec(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wr_elec,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  call MPI_ALLGATHER(wr_ion(:,i_r-1),&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       wr_ion,&
       tgyro_stab_nky,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  print *

  ktag = 'ky='

  if (i_proc_global == 0) then

     open(unit=1,file='wi_elec.out',status='replace')

     write(1,'(t10,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wi_elec(:,i)
     enddo

     close(1)

     open(unit=1,file='wi_ion.out',status='replace')

     write(1,'(t10,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wi_ion(:,i)
     enddo

     close(1)

     open(unit=1,file='wr_elec.out',status='replace')

     write(1,'(t10,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wr_elec(:,i)
     enddo

     close(1)

     open(unit=1,file='wr_ion.out',status='replace')

     write(1,'(t10,a)') 'ky rho'
     write(1,20) 'r/a',ky(:)
     do i=1,n_inst
        write(1,10) r(i+1)/r_min,wr_ion(:,i)
     enddo

     close(1)

  endif

10 format(f5.3,3x,10(1pe12.5,1x))
20 format(t2,a,5x,10(f5.3,8x))

end subroutine tgyro_stab_driver
