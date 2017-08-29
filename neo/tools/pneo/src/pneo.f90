program pneo

  use mpi
  use pneo_globals
  use neo_interface
  use EXPRO_interface

  implicit none

  integer :: p, i1,i2,i3,i4,i5,i6,i7,i8
  integer :: ni1=4, ni2=4, ni3=4, ni4=4, ni5=4, ni6=4, ni7=4, ni8=4
  real, dimension(4) :: rmin_over_rmaj = (/ 0.1,0.2,0.3,0.4 /)
  real, dimension(4) :: q = (/ 1.0,2.0,3.0,4.0 /)
  real, dimension(4) :: nu_ee = (/ 1.0e-4,1.0e-3,1.0e-2,1.0e-1 /)
  real, dimension(4) :: ni_over_ne = (/ 0.25,0.5,0.75,1.0 /)
  real, dimension(4) :: ti_over_te  = (/ 0.75,1.0,1.5,2.0 /)
  real, dimension(4) :: delta   = (/ 0.0,0.1,0.2,0.3 /)
  real, dimension(4) :: s_delta = (/ 0.0,0.1,0.2,0.3 /)
  real, dimension(4) :: s_kappa = (/ 0.0,0.1,0.2,0.3 /)
  
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
  
  ntot = ni1*ni2*ni3*ni4*ni5*ni6*ni7*ni8
  allocate(ic1(ntot))
  allocate(ic2(ntot))
  allocate(ic3(ntot))
  allocate(ic4(ntot))
  allocate(ic5(ntot))
  allocate(ic6(ntot))
  allocate(ic7(ntot))
  allocate(ic8(ntot))
  allocate(indata_vec(8,ntot))
  allocate(indata_tot(8,ntot))
  allocate(outdata_vec(6,ntot))
  allocate(outdata_tot(6,ntot))

  p = 0
  do i1=1,ni1
     do i2=1,ni2
        do i3=1,ni3
           do i4=1,ni4
              do i5=1,ni5
                 do i6=1,ni6
                    do i7=1,ni7
                       do i8=1,ni8
                      
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
  
  do p=1+i_proc,ntot,n_proc

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

     indata_vec(1,p) = neo_rmin_over_a_in
     indata_vec(2,p) = neo_q_in
     indata_vec(3,p) = neo_nu_1_in
     indata_vec(4,p) = neo_dens_in(2)
     indata_vec(5,p) = neo_temp_in(2)
     indata_vec(6,p) = neo_delta_in
     indata_vec(7,p) = neo_s_delta_in
     indata_vec(8,p) = neo_s_kappa_in
     
     ! Cne
     neo_dlnndr_in(:) = 0.0; neo_dlnndr_in(1) = 1.0
     !call neo_run()
     outdata_vec(1,p) = neo_jpar_dke_out
     ! CTe
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(1) = 1.0 
     !call neo_run()
     outdata_vec(2,p) = neo_jpar_dke_out
     ! Cni1
     neo_dlnndr_in(:) = 0.0; neo_dlnndr_in(2) = 1.0 
     !call neo_run()
     outdata_vec(3,p) = neo_jpar_dke_out
     ! CTi1
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(2) = 1.0      
     !call neo_run()
     outdata_vec(4,p) = neo_jpar_dke_out
     ! Cni2
     neo_dlnndr_in(:) = 0.0; neo_dlnndr_in(3) = 1.0 
     !call neo_run()
     outdata_vec(5,p) = neo_jpar_dke_out
     ! CTi2 
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(3) = 1.0    
     !call neo_run()
     outdata_vec(6,p) = neo_jpar_dke_out
     
  enddo

  ! Collect all data 
  call MPI_ALLREDUCE(indata_vec,indata_tot,size(indata_tot), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(outdata_vec,outdata_tot,size(outdata_tot), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)

  if (i_proc == 0) then
     print *,indata_tot(8,:)
  endif

  call MPI_finalize(i_err)

end program pneo
