program pneo

  use mpi
  use pneo_globals
  use neo_interface
  use EXPRO_interface

  implicit none

  integer :: p,k,j
  integer :: i1,i2,i3,i4,i5,i6,i7,i8
  integer,parameter :: n1=14, n2=10, n3=10, n4=4, n5=6, n6=6, n7=10, n8=10
  real, dimension(n1) :: rmin_over_rmaj = (/ 0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.65,0.70,0.75 /)
  real, dimension(n2) :: q = (/ 1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5 /)
  real, dimension(n3) :: nu_ee = (/ 1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1e0,5e0 /)
  real, dimension(n4) :: ni_over_ne = (/ 0.7,0.8,0.9,0.99 /)
  real, dimension(n5) :: ti_over_te  = (/ 1.0,1.2,1.4,1.6,1.8,2.0 /)
  real, dimension(n6) :: delta   = (/ 0.0,0.1,0.2,0.3,0.4,0.5 /)
  real, dimension(n7) :: s_delta = (/ 0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8 /)
  real, dimension(n8) :: s_kappa = (/ 0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8 /)
  integer, parameter ::  n_sub = 5
  real, dimension(n_sub) :: r_sub
  integer, dimension(n_sub) :: i_sub
  integer :: values(1:8)
  integer, dimension(:), allocatable :: seed

  !---------------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------------

  ! Path is cwd:
  path= './'

  call neo_init_serial(path)

  ! pointers

  ntot = n1*n2*n3*n4*n5*n6*n7*n8
  allocate(ic1(ntot))
  allocate(ic2(ntot))
  allocate(ic3(ntot))
  allocate(ic4(ntot))
  allocate(ic5(ntot))
  allocate(ic6(ntot))
  allocate(ic7(ntot))
  allocate(ic8(ntot))

  p = 0
  do i1=1,n1
     do i2=1,n2
        do i3=1,n3
           do i4=1,n4
              do i5=1,n5
                 do i6=1,n6
                    do i7=1,n7
                       do i8=1,n8

                          p = p+1
                          ic1(p) = i1
                          ic2(p) = i2
                          ic3(p) = i3
                          ic4(p) = i4
                          ic5(p) = i5
                          ic6(p) = i6
                          ic7(p) = i7
                          ic8(p) = i8
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo

  ! Subset the data to be computed by computing random p values
  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)
  call random_number(r_sub)
  i_sub = floor(r_sub*ntot)+1

  allocate(indata_loc(8,n_sub))
  allocate(indata(8,n_sub))
  allocate(ingeodata_loc(12,n_sub))
  allocate(ingeodata(12,n_sub))
  allocate(outdata_loc(6,n_sub))
  allocate(outdata(6,n_sub))

  ! Fixed NEO subroutine inputs
  neo_silent_flag_in = 1
  neo_n_energy_in = 6
  neo_n_xi_in     = 17
  neo_n_theta_in  = 17
  neo_collision_model_in = 4   ! Full FP collisions
  neo_rho_star_in = 0.001      ! arbitrary

  ! geometry
  neo_equilibrium_model_in = 2 ! Miller equilibrium
  neo_rmaj_over_a_in = 1.0     ! anorm = rmaj
  neo_kappa_in = 1.0
  neo_zeta_in  = 0.0
  neo_s_zeta_in = 0.0
  neo_zmag_over_a_in = 0.0
  neo_s_zmag_in = 0.0

  ! 3 species: ele + D + C
  neo_n_species_in = 3
  ! electrons
  neo_z_in(1)     = -1
  neo_mass_in(1)  = 0.0002724486
  neo_temp_in(1)  = 1.0        ! Tnorm = Te
  neo_dens_in(1)  = 1.0        ! nnorm = ne
  ! D
  neo_z_in(2)     = 1
  neo_mass_in(2)  = 1.0
  ! C
  neo_z_in(3)     = 6
  neo_mass_in(3)  = 6.0

!!!!!! FOR TESTING, USE THEORY sim_model=0; else use sim_model=1
!  neo_sim_model_in = 0
!!!!!!

  do j=1+i_proc,n_sub,n_proc

     p = i_sub(j)
     print *,p

     i1 = ic1(p) 
     i2 = ic2(p)
     i3 = ic3(p) 
     i4 = ic4(p)
     i5 = ic5(p) 
     i6 = ic6(p)
     i7 = ic7(p) 
     i8 = ic8(p)

     neo_rmin_over_a_in = rmin_over_rmaj(i1)
     neo_q_in       = q(i2)
     neo_nu_1_in    = nu_ee(i3)
     neo_dens_in(2) = ni_over_ne(i4)
     neo_dens_in(3) = (1.0-neo_z_in(2)*neo_dens_in(2))/(1.0*neo_z_in(3))
     neo_temp_in(2) = ti_over_te(i5)
     neo_temp_in(3) = ti_over_te(i5)
     neo_delta_in   = delta(i6)
     neo_s_delta_in = s_delta(i7)
     neo_s_kappa_in = s_kappa(i8)

     indata_loc(1,j) = neo_rmin_over_a_in
     indata_loc(2,j) = neo_q_in
     indata_loc(3,j) = neo_nu_1_in
     indata_loc(4,j) = neo_dens_in(2)
     indata_loc(5,j) = neo_temp_in(2)
     indata_loc(6,j) = neo_delta_in
     indata_loc(7,j) = neo_s_delta_in
     indata_loc(8,j) = neo_s_kappa_in

     ! Cne
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(:) = 0.0; neo_dlnndr_in(1) = 1.0
     !call neo_run()
     outdata_loc(1,j) = neo_jpar_dke_out/(abs(neo_z_in(1))*neo_dens_in(1))
     ! CTe
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlntdr_in(1) = 1.0 
     !call neo_run()
     outdata_loc(2,j) = neo_jpar_dke_out/(abs(neo_z_in(1))*neo_dens_in(1))
     ! Cni1
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlnndr_in(2) = 1.0 
     !call neo_run()
     outdata_loc(3,j) = neo_jpar_dke_out/(abs(neo_z_in(2))*neo_dens_in(2))
     ! CTi1
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlntdr_in(2) = 1.0  
     !call neo_run()
     outdata_loc(4,j) = neo_jpar_dke_out/(abs(neo_z_in(2))*neo_dens_in(2))
     ! Cni2
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlnndr_in(3) = 1.0 
     !call neo_run()
     outdata_loc(5,j) = neo_jpar_dke_out/(abs(neo_z_in(3))*neo_dens_in(3))
     ! CTi2 
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlntdr_in(3) = 1.0    
     !call neo_run()
     outdata_loc(6,j) = neo_jpar_dke_out/(abs(neo_z_in(3))*neo_dens_in(3))

     ! store geometry parameters
     do k=1,12
        ingeodata_loc(k,j) = neo_geoparams_out(k)
     enddo

     ! renormalize the coefficients by 1/I_div_psiprime/rho_star
     outdata_loc(:,j) = outdata_loc(:,j) &
          / (neo_geoparams_out(1)*neo_rho_star_in)

  enddo

  ! Collect all data 
  call MPI_ALLREDUCE(indata_loc,indata,size(indata), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(ingeodata_loc,ingeodata,size(ingeodata), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(outdata_loc,outdata,size(outdata), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)

  if (i_proc == 0) then

     print *,'NTOT=',ntot
     print *,'NSUB=',n_sub

     open(unit=1,file='indata.dat',status='replace')
     write(1,10) indata(:,:)
     close(1)

     open(unit=1,file='ingeodata.dat',status='replace')
     write(1,10) ingeodata(:,:)
     close(1)

     open(unit=1,file='outdata.dat',status='replace')
     write(1,10) outdata(:,:)
     close(1)

  endif

  call MPI_finalize(i_err)

10 format(1pe12.5)

end program pneo
