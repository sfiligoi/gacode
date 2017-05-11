!------------------------------------------------------------------------
! globalmap.f90
!
! PURPOSE: 
!------------------------------------------------------------------------

program globalmap

  use mpi
  use globalmap_globals
  use cgyro_globals
  use cgyro_io
  use GEO_interface

  implicit none

  integer :: supported
  integer :: ir,it
  real :: rho_star,dr,ri,lr

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator
  call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,i_err)
  !-----------------------------------------------------------------

  ! Path is cwd:
  path= './'

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  CGYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,i_err)
  !-----------------------------------------------------------------

  call cgyro_read_input
  open(unit=1,file='input.globalmap',status='old')
  read(1,*) rho_star
  close(1)
  
  test_flag   = 1
  silent_flag = 1
  call cgyro_kernel

  lr = rho_star/rho*length
  print *,'Lr',lr

  d_theta = (2*pi/n_theta)
  do it=1,n_theta
     theta(it) = -pi+(it-1)*d_theta
  enddo

  do ir=1,n_radial

       ri = rmin + (ir-(n_radial+1.0)/2)*lr/(n_radial-1)
       dr = ri-rmin
       
       print *,ri,dr

      GEO_rmin_in  = ri
      GEO_rmaj_in  = rmaj+shift*dr
      GEO_drmaj_in = shift
      GEO_zmag_in  = zmag+dzmag*dr
      GEO_dzmag_in = dzmag
      GEO_q_in     = q*(1+s*dr/rmin)
      GEO_s_in     = s
      GEO_kappa_in = kappa*(1+s_kappa*dr/rmin)
      GEO_s_kappa_in = s_kappa
      !  GEO_delta_in   = delta_s(ir_norm)+s_delta_s(ir_norm)/r(ir_norm)*dr
      !  GEO_s_delta_in = s_delta_s(ir_norm)
      !  GEO_zeta_in    = zeta_s(ir_norm)+s_zeta_s(ir_norm)/r(ir_norm)*dr
      !  GEO_s_zeta_in  = s_zeta_s(ir_norm)
      !  GEO_beta_star_in = beta_star_s(ir_norm)

     !rmin = ra + (ir-1)*(rb-ra)/(n_radial

     GEO_rmin_in      = rmin
     GEO_rmaj_in      = rmaj
     GEO_drmaj_in     = shift
     GEO_zmag_in      = zmag
     GEO_dzmag_in     = dzmag
     GEO_q_in         = q
     GEO_s_in         = s
     GEO_kappa_in     = kappa
     GEO_s_kappa_in   = s_kappa
     GEO_delta_in     = delta
     GEO_s_delta_in   = s_delta
     GEO_zeta_in      = zeta
     GEO_s_zeta_in    = s_zeta
     GEO_beta_star_in = beta_star(0)

     call GEO_do()  

     do it=1,n_theta

        call GEO_interp(theta(it))     

        !print *,GEO_bigr

     enddo
  enddo

  call MPI_FINALIZE(i_err)

end program globalmap
