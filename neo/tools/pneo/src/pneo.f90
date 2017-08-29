program pneo

  use mpi
  use pneo_globals
  use neo_interface
  use EXPRO_interface

  implicit none

  integer :: p, k, i1,i2,i3,i4,i5,i6,i7,i8
  integer,parameter :: n1=4, n2=4, n3=1, n4=1, n5=1, n6=1, n7=1, n8=1
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
  
  ntot = n1*n2*n3*n4*n5*n6*n7*n8
  allocate(ic1(ntot))
  allocate(ic2(ntot))
  allocate(ic3(ntot))
  allocate(ic4(ntot))
  allocate(ic5(ntot))
  allocate(ic6(ntot))
  allocate(ic7(ntot))
  allocate(ic8(ntot))
  allocate(indata_loc(8,ntot))
  allocate(indata(8,ntot))
  allocate(ingeodata_loc(12,ntot))
  allocate(ingeodata(12,ntot))
  allocate(outdata_loc(6,ntot))
  allocate(outdata(6,ntot))

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
  neo_sim_model_in = 0
  !!!!!!
  
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

     indata_loc(1,p) = neo_rmin_over_a_in
     indata_loc(2,p) = neo_q_in
     indata_loc(3,p) = neo_nu_1_in
     indata_loc(4,p) = neo_dens_in(2)
     indata_loc(5,p) = neo_temp_in(2)
     indata_loc(6,p) = neo_delta_in
     indata_loc(7,p) = neo_s_delta_in
     indata_loc(8,p) = neo_s_kappa_in
     
     ! Cne
     neo_dlnndr_in(:) = 0.0; neo_dlnndr_in(1) = 1.0
     call neo_run()
     outdata_loc(1,p) = neo_jpar_dke_out/(abs(neo_z_in(1))*neo_dens_in(1))
     ! CTe
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(1) = 1.0 
     call neo_run()
     outdata_loc(2,p) = neo_jpar_dke_out/(abs(neo_z_in(1))*neo_dens_in(1))
     ! Cni1
     neo_dlnndr_in(:) = 0.0; neo_dlnndr_in(2) = 1.0 
     call neo_run()
     outdata_loc(3,p) = neo_jpar_dke_out/(abs(neo_z_in(2))*neo_dens_in(2))
     ! CTi1
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(2) = 1.0      
     call neo_run()
     outdata_loc(4,p) = neo_jpar_dke_out/(abs(neo_z_in(2))*neo_dens_in(2))
     ! Cni2
     neo_dlnndr_in(:) = 0.0; neo_dlnndr_in(3) = 1.0 
     call neo_run()
     outdata_loc(5,p) = neo_jpar_dke_out/(abs(neo_z_in(3))*neo_dens_in(3))
     ! CTi2 
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(3) = 1.0    
     call neo_run()
     outdata_loc(6,p) = neo_jpar_dke_out/(abs(neo_z_in(3))*neo_dens_in(3))

     ! store geometry parameters
     do k=1,12
        ingeodata_loc(k,p) = neo_geoparams_out(k)
     enddo

     ! renormalize the coefficients by 1/I_div_psiprime/rho_star
     outdata_loc(:,p) = outdata_loc(:,p) &
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
     do p=1,ntot
        print *, p, indata(1,p), ingeodata(2,p)
     enddo
  endif

  call MPI_finalize(i_err)

end program pneo
