!--------------------------------------------------------------------
! vgen.f90
!
! PURPOSE: 
!  Driver for the vgen (velocity-generation) capability of NEO.  This 
!  will write a new input.profiles with NEO-computer electric field 
!  and/or velocities, and also generate a new input.profiles.extra.
!--------------------------------------------------------------------

program vgen

  use mpi
  use vgen_globals
  use neo_interface
  use EXPRO_interface

  !---------------------------------------------------
  implicit none

  character(len=80) :: vtag

  integer :: i
  integer :: j
  integer :: ix
  integer :: ia,ib
  integer :: adiabatic_ele_model
  integer :: num_ele
  integer :: indx_ele
  integer :: n_ions
  integer :: rotation_model 
  real :: grad_p
  real :: vtor_er
  real :: ya
  real :: yb
  real :: vtor_diff
  real :: er0
  real :: omega
  real :: omega_deriv

  real, dimension(:), allocatable :: er_exp
  real, dimension(:), allocatable :: jbs_neo
  real, dimension(:), allocatable :: jbs_sauter
  real, dimension(:), allocatable :: jbs_koh
  real, dimension(:), allocatable :: jbs_nclass
  real :: jbs_norm_fac
  real, dimension(:), allocatable :: pflux_sum
  integer :: zfac

  !---------------------------------------------------

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Path is cwd:
  !
  path= './'
  !-----------------------------------------------------------------

  !-------------------------------------------------
  ! Obscure definition of number tags 
  do ix=1,100
     write (tag(ix),fmt) ix-1
  enddo
  !-------------------------------------------------

  !---------------------------------------------------------------------
  ! Read vgen control parameters
  !
  open(unit=1,file='vgen.dat',status='old')
  read(1,*) er_method
  read(1,*) vel_method
  read(1,*) erspecies_indx
  read(1,*) nth_min
  read(1,*) nth_max
  close(1)

  select case(er_method)
  case(1)
     if(i_proc == 0) then
        print '(a)','INFO: (VGEN) Computing omega0 (Er)  from force balance'
     endif
  case(2)
     if(i_proc == 0) then
        print '(a)', 'INFO: (VGEN) Computing omega0 (Er) from NEO (weak rotation limit)'
     endif
  case(4)
     if(i_proc == 0) then
        print '(a)','INFO: (VGEN) Returning given omega0 (Er)'
     endif
  case default
     if(i_proc == 0) then
        print '(a)','ERROR: Invalid er_method'
     endif
     call MPI_finalize(i_err)
     stop
  end select

  select case(vel_method)
  case(1)
     if(i_proc == 0) then
        print '(a)','INFO: (VGEN) Computing velocities from NEO (weak rotation limit)'
     endif
  case(2)
     if(i_proc == 0) then
        print '(a)','INFO: (VGEN) Computing velocities from NEO (strong rotation limit)'
     endif
  case default
     if(i_proc == 0) then
        print '(a)','ERROR: Invalid vel_method'
     endif
     call MPI_finalize(i_err)
     stop
  end select

  if(i_proc == 0) then
     print '(a,i2,a,i2)','INFO: (VGEN) Using NEO Theta Resolution:',&
          nth_min,',',nth_max
  endif

  !---------------------------------------------------------------------

  call neo_init_serial(path)
  call neo_read_input()
  call map_global2interface()

  neo_n_radial_in = 1
  neo_profile_model_in = 1

  ! Species checks

  num_ele = 0
  indx_ele = 0
  if(neo_n_species_in >= 1 .and. neo_z_1_in == -1) then
     num_ele  = num_ele + 1
     indx_ele = 1
  endif
  if(neo_n_species_in >= 2 .and. neo_z_2_in == -1) then
     num_ele  = num_ele + 1
     indx_ele = 2
  endif
  if(neo_n_species_in >= 3 .and. neo_z_3_in == -1) then
     num_ele  = num_ele + 1
     indx_ele = 3
  endif
  if(neo_n_species_in >= 4 .and. neo_z_4_in == -1) then
     num_ele  = num_ele + 1
     indx_ele = 4
  endif
  if(neo_n_species_in >= 5 .and. neo_z_5_in == -1) then
     num_ele  = num_ele + 1
     indx_ele = 5
  endif
  if(neo_n_species_in >= 6 .and. neo_z_6_in == -1) then
     num_ele  = num_ele + 1
     indx_ele = 6
  endif

  if(num_ele == 0) then
     adiabatic_ele_model = 1
     n_ions = neo_n_species_in
  else if(num_ele == 1) then
     adiabatic_ele_model = 0
     n_ions = neo_n_species_in - 1
     if(indx_ele /= neo_n_species_in) then
        if(i_proc == 0) then
           print '(a)','ERROR: (VGEN) Electron species must be n_species'
        endif
        call MPI_finalize(i_err)
        stop
     endif
  else
     if(i_proc == 0) then
        print '(a)', 'ERROR: (VGEN) Only one electron species allowed'
     endif
     call MPI_finalize(i_err)
     stop
  end if

  if(n_ions < 1) then
     if(i_proc == 0) then
        print '(a)', 'ERROR: (VGEN) There must be at least one ion species'
     endif
     call MPI_finalize(i_err)
     stop
  endif

  if(erspecies_indx > n_ions) then
     if(i_proc == 0) then
        print '(a)', 'ERROR: (VGEN) Invalid species index'
     endif
     call MPI_finalize(i_err)
     stop
  endif

  !---------------------------------------------------------------------
  ! Read experimental profiles using EXPRO library.
  !
  call EXPRO_palloc(MPI_COMM_WORLD,path,1)
  !
  ! Reset species 1 density for quasineutrality
  !
  EXPRO_ctrl_rotation_method = 1
  EXPRO_ctrl_density_method = 2 
  EXPRO_ctrl_z(:) = 0.0
  EXPRO_ctrl_z(1) = 1.0*neo_z_1_in
  !
  if (neo_n_species_in >= 2 .and. neo_z_2_in /= -1) then
     EXPRO_ctrl_z(2) = 1.0*neo_z_2_in
  endif
  if (neo_n_species_in >= 3 .and. neo_z_3_in /= -1) then
     EXPRO_ctrl_z(3) = 1.0*neo_z_3_in
  endif
  if (neo_n_species_in >= 4 .and. neo_z_4_in /= -1) then
     EXPRO_ctrl_z(4) = 1.0*neo_z_4_in
  endif
  if (neo_n_species_in >= 5 .and. neo_z_5_in /= -1) then
     EXPRO_ctrl_z(5) = 1.0*neo_z_5_in
  endif
  !
  ! Set equilibrium option for EXPRO
  !
  if (neo_equilibrium_model_in == 3) then
     EXPRO_ctrl_numeq_flag = 1
  else
     EXPRO_ctrl_numeq_flag = 0
  endif

  call EXPRO_pread

  ! Write the derived quantities to input.profiles.extra

  if (i_proc == 0) call EXPRO_write_derived
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Set sign of btccw and ipccw from sign of b and q from EXPRO
  !
  neo_btccw_in = -EXPRO_signb
  neo_ipccw_in = -EXPRO_signb*EXPRO_signq
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Storage for electric field at theta=0 (er0) 
  allocate(er_exp(EXPRO_n_exp))
  allocate(vtor_measured(EXPRO_n_exp))

  ! Storage for bootstrap current calculations
  allocate(jbs_neo(EXPRO_n_exp))
  allocate(jbs_sauter(EXPRO_n_exp))
  allocate(jbs_koh(EXPRO_n_exp))
  allocate(jbs_nclass(EXPRO_n_exp))
  allocate(pflux_sum(EXPRO_n_exp))
  jbs_neo(:)    = 0.0
  jbs_sauter(:) = 0.0
  jbs_koh(:)    = 0.0
  jbs_nclass(:) = 0.0
  pflux_sum(:)  = 0.0

  if (er_method /= 4) then
     er_exp(:) = 0.0
     EXPRO_w0(:) = 0.0
     EXPRO_w0p(:) = 0.0
  endif
  !---------------------------------------------------------------------

  do j=1,5
     if (j == erspecies_indx) then
        vtor_measured(:) = EXPRO_vtor(j,:)
     endif
     EXPRO_vpol(j,:) = 0.0
     EXPRO_vtor(j,:) = 0.0
  enddo
  jbs_neo(:)    = 0.0
  jbs_sauter(:) = 0.0
  jbs_koh(:)    = 0.0
  jbs_nclass(:) = 0.0
  pflux_sum(:)  = 0.0
  !======================================================================
  ! Four alternatives for Er calculation:
  !
  ! 1. Compute Er from force balance using measured vpol and vtor 
  ! 2. Compute Er by matching vtor_measured with vtor_neo at theta=0 assuming
  !    weak rotation
  ! 3. Compute Er by matching vtor_measured with vtor_neo at theta=0 assuming
  !    strong rotation
  ! 4. Return the given Er

  select case (er_method) 

  case (1,4)

     if (er_method == 1) then

        ! Compute Er from force balance using measured vpol and vtor

        ! Er calculation first

        open(unit=1,file='out.vgen.ercomp'//tag(i_proc+1),status='replace')
        close(1)
        do i=2+i_proc,EXPRO_n_exp-1,n_proc
           if (erspecies_indx == 1) then
              grad_p = -(EXPRO_dlnnidr_new(i) + EXPRO_dlntidr(1,i))
           else
              grad_p = -(EXPRO_dlnnidr(erspecies_indx,i) &
                   + EXPRO_dlntidr(erspecies_indx,i))
           endif
           er_exp(i) = (grad_p * EXPRO_grad_r0(i) &
                * EXPRO_ti(erspecies_indx,i)*temp_norm_fac &
                / (EXPRO_ctrl_z(erspecies_indx) * charge_norm_fac) & 
                + EXPRO_vtor(erspecies_indx,i) * EXPRO_bp0(i) &
                - EXPRO_vpol(erspecies_indx,i) * EXPRO_bt0(i)) &
                / 1000
           open(unit=1,file='out.vgen.ercomp'//tag(i_proc+1),status='old',position='append')
           write(1,'(e16.8)',advance='no') EXPRO_rho(i)
           write(1,'(e16.8)',advance='no') grad_p * EXPRO_grad_r0(i) &
                * EXPRO_ti(erspecies_indx,i)*temp_norm_fac &
                / (EXPRO_ctrl_z(erspecies_indx) * charge_norm_fac) / 1000
           write(1,'(e16.8)',advance='no') EXPRO_vtor(erspecies_indx,i) * EXPRO_bp0(i)/1000
           write(1,'(e16.8)',advance='no') -EXPRO_vpol(erspecies_indx,i) * EXPRO_bt0(i)/1000
           write(1,*) 
           close(1)
        enddo

        call vgen_reduce(er_exp(2:EXPRO_n_exp-1),EXPRO_n_exp-2)

        ! Compute omega and omega_deriv from newly-generated Er:

        do i=2,EXPRO_n_exp-1
           EXPRO_w0(i) = 2.9979e10*EXPRO_q(i)*(er_exp(i)/30.0)/ &
                ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
        enddo

        ! Compute w0p
        call bound_deriv(EXPRO_w0p(2:EXPRO_n_exp-1),EXPRO_w0(2:EXPRO_n_exp-1),&
             EXPRO_rmin,EXPRO_n_exp-2)

     endif

     if (er_method == 4) then

        do i=1,EXPRO_n_exp
           er_exp(i) = 1.0/2.9979e10/EXPRO_q(i)*EXPRO_w0(i)*30.0* &
                ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
        enddo

     endif

     ! Flow calculation based on existing EXPRO_w0 and EXPRO_w0p

     do i=2+i_proc,EXPRO_n_exp-1,n_proc

        if (vel_method == 1) then
           rotation_model = 1   ! weak rotation
        else
           rotation_model = 2   ! strong rotation
        endif
        er0 = er_exp(i)
        omega = EXPRO_w0(i) 
        omega_deriv = EXPRO_w0p(i) 

        call vgen_compute_neo(i,vtor_diff,rotation_model,er0,omega,omega_deriv)

        if (neo_error_status_out > 0) then
           print *,neo_error_message_out
           stop
        endif

        do j=1,n_ions
           EXPRO_vpol(j,i) = neo_vpol_dke_out(j) &
                * vth_norm * EXPRO_rmin(EXPRO_n_exp)
           EXPRO_vtor(j,i) = neo_vtor_dke_out(j) &
                * vth_norm * EXPRO_rmin(EXPRO_n_exp)
        enddo
        jbs_norm_fac = charge_norm_fac*dens_norm*vth_norm &
             *EXPRO_rmin(EXPRO_n_exp)/1e6
        jbs_neo(i)    = neo_jpar_dke_out*jbs_norm_fac
        jbs_sauter(i) = neo_jpar_thS_out*jbs_norm_fac
        jbs_koh(i)    = neo_jpar_thK_out*jbs_norm_fac
        jbs_nclass(i) = neo_jpar_thN_out*jbs_norm_fac
        pflux_sum(i)  = 0.0
        do j=1,neo_n_species_in
           if(j==1) zfac=neo_z_1_in
           if(j==2) zfac=neo_z_2_in
           if(j==3) zfac=neo_z_3_in
           if(j==4) zfac=neo_z_4_in
           if(j==5) zfac=neo_z_5_in
           if(j==6) zfac=neo_z_6_in
           pflux_sum(i) = pflux_sum(i) + zfac*neo_pflux_dke_out(j)
        enddo
        pflux_sum(i) = pflux_sum(i) / neo_rho_star_in**2

        print 10,EXPRO_rho(i),&
             er_exp(i),EXPRO_vtor(1,i)/1e3,EXPRO_vpol(1,i)/1e3

     enddo

     ! Reduce vpol,vtor
     do j=1,n_ions
        call vgen_reduce(EXPRO_vpol(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
        call vgen_reduce(EXPRO_vtor(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
     enddo

  case (2)

     ! Compute Er using NEO (weak rotation limit) 
     ! by matching vtor_measured with vtor_neo at theta=0 

     do i=2+i_proc,EXPRO_n_exp-1,n_proc

        rotation_model = 1
        er0 = 0.0
        omega = 0.0
        omega_deriv = 0.0

        call vgen_compute_neo(i,vtor_diff, rotation_model, er0, omega, &
             omega_deriv)

        if (neo_error_status_out > 0) then
           print *,neo_error_message_out
           stop
        endif

        ! omega = (vtor_measured - vtor_neo_ater0) / R

        er_exp(i) = (vtor_diff/(vth_norm * EXPRO_rmin(EXPRO_n_exp)))  &
             / (neo_rmaj_over_a_in + neo_rmin_over_a_in) &
             * neo_rmin_over_a_in / (abs(neo_q_in) * EXPRO_signq) &
             / (abs(neo_rho_star_in) * EXPRO_signb) &
             * EXPRO_grad_r0(i) &
             * (temp_norm*temp_norm_fac / charge_norm_fac &
             / EXPRO_rmin(EXPRO_n_exp)) / 1000

        ! Part of vtor due to Er
        vtor_er = vtor_diff

        ! Store the new flows
        do j=1,n_ions
           EXPRO_vpol(j,i) = neo_vpol_dke_out(j) &
                * vth_norm * EXPRO_rmin(EXPRO_n_exp)
           EXPRO_vtor(j,i) = neo_vtor_dke_out(j) &
                * vth_norm * EXPRO_rmin(EXPRO_n_exp) + vtor_er
        enddo
        jbs_norm_fac = charge_norm_fac*dens_norm*vth_norm &
             *EXPRO_rmin(EXPRO_n_exp)/1e6
        jbs_neo(i)    = neo_jpar_dke_out*jbs_norm_fac
        jbs_sauter(i) = neo_jpar_thS_out*jbs_norm_fac
        jbs_koh(i)    = neo_jpar_thK_out*jbs_norm_fac
        jbs_nclass(i) = neo_jpar_thN_out*jbs_norm_fac
        pflux_sum(i)  = 0.0
        do j=1,neo_n_species_in
           if(j==1) zfac=neo_z_1_in
           if(j==2) zfac=neo_z_2_in
           if(j==3) zfac=neo_z_3_in
           if(j==4) zfac=neo_z_4_in
           if(j==5) zfac=neo_z_5_in
           if(j==6) zfac=neo_z_6_in
           pflux_sum(i) = pflux_sum(i) + zfac*neo_pflux_dke_out(j)
        enddo
        pflux_sum(i) = pflux_sum(i) / neo_rho_star_in**2

        print 10,EXPRO_rho(i),&
             er_exp(i),EXPRO_vtor(1,i)/1e3,EXPRO_vpol(1,i)/1e3

     enddo

     ! Reduce er,vpol,vtor
     call vgen_reduce(er_exp(2:EXPRO_n_exp-1),EXPRO_n_exp-2)

     if (vel_method == 1) then

        do j=1,n_ions
           call vgen_reduce(EXPRO_vpol(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
           call vgen_reduce(EXPRO_vtor(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
        enddo

     else

        EXPRO_vpol(:,:) = 0.0
        EXPRO_vtor(:,:) = 0.0

        if (i_proc == 0) print '(a)', 'INFO: (VGEN) Recomputing flows using NEO in the strong rotation limit.'

        ! Re-compute the flows using strong rotation
        ! omega and omega_deriv 
        do i=2,EXPRO_n_exp-1
           EXPRO_w0(i) = 2.9979e10*EXPRO_q(i)*(er_exp(i)/30.0)/ &
                ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
        enddo
        ! Compute w0p
        call bound_deriv(EXPRO_w0p(2:EXPRO_n_exp-1),EXPRO_w0(2:EXPRO_n_exp-1),&
             EXPRO_rmin,EXPRO_n_exp-2)

        do i=2+i_proc,EXPRO_n_exp-1,n_proc
           rotation_model = 2  
           er0 = er_exp(i)
           omega = EXPRO_w0(i) 
           omega_deriv = EXPRO_w0p(i) 
           call vgen_compute_neo(i,vtor_diff, rotation_model, er0, omega, &
                omega_deriv)
           if (neo_error_status_out > 0) then
              print *,neo_error_message_out
              stop
           endif
           do j=1,n_ions
              EXPRO_vpol(j,i) = neo_vpol_dke_out(j) &
                   * vth_norm * EXPRO_rmin(EXPRO_n_exp)
              EXPRO_vtor(j,i) = neo_vtor_dke_out(j) &
                   * vth_norm * EXPRO_rmin(EXPRO_n_exp)
           enddo
           jbs_norm_fac = charge_norm_fac*dens_norm*vth_norm &
                *EXPRO_rmin(EXPRO_n_exp)/1e6
           jbs_neo(i)    = neo_jpar_dke_out*jbs_norm_fac
           jbs_sauter(i) = neo_jpar_thS_out*jbs_norm_fac
           jbs_koh(i)    = neo_jpar_thK_out*jbs_norm_fac
           jbs_nclass(i) = neo_jpar_thN_out*jbs_norm_fac
           pflux_sum(i)  = 0.0
           do j=1,neo_n_species_in
              if(j==1) zfac=neo_z_1_in
              if(j==2) zfac=neo_z_2_in
              if(j==3) zfac=neo_z_3_in
              if(j==4) zfac=neo_z_4_in
              if(j==5) zfac=neo_z_5_in
              if(j==6) zfac=neo_z_6_in
              pflux_sum(i) = pflux_sum(i) + zfac*neo_pflux_dke_out(j)
           enddo
           pflux_sum(i) = pflux_sum(i) / neo_rho_star_in**2

           print 10,EXPRO_rho(i),&
                er_exp(i),EXPRO_vtor(1,i)/1e3,EXPRO_vpol(1,i)/1e3

        enddo

        ! Reduce vpol,vtor
        do j=1,n_ions
           call vgen_reduce(EXPRO_vpol(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
           call vgen_reduce(EXPRO_vtor(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
        enddo

     endif

  end select
  !======================================================================

  ! Additional reductions
  call vgen_reduce(pflux_sum(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jbs_neo(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jbs_sauter(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jbs_nclass(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jbs_koh(2:EXPRO_n_exp-1),EXPRO_n_exp-2)

  ! extrapolation for r=0 and r=n_exp boundary points

  if(er_method /= 4) then
     call bound_extrap(ya,yb,er_exp,EXPRO_rmin,EXPRO_n_exp)
     er_exp(1) = ya
     er_exp(EXPRO_n_exp) = yb
  endif

  do j=1,n_ions
     call bound_extrap(ya,yb,EXPRO_vpol(j,:),EXPRO_rmin,EXPRO_n_exp)
     EXPRO_vpol(j,1) = ya
     EXPRO_vpol(j,EXPRO_n_exp) = yb
     call bound_extrap(ya,yb,EXPRO_vtor(j,:),EXPRO_rmin,EXPRO_n_exp)
     EXPRO_vtor(j,1) = ya
     EXPRO_vtor(j,EXPRO_n_exp) = yb
  enddo
  call bound_extrap(ya,yb,jbs_neo,EXPRO_rmin,EXPRO_n_exp)
  jbs_neo(1)           = ya
  jbs_neo(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jbs_sauter,EXPRO_rmin,EXPRO_n_exp)
  jbs_sauter(1)           = ya
  jbs_sauter(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jbs_koh,EXPRO_rmin,EXPRO_n_exp)
  jbs_koh(1)           = ya
  jbs_koh(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jbs_nclass,EXPRO_rmin,EXPRO_n_exp)
  jbs_nclass(1)           = ya
  jbs_nclass(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,pflux_sum,EXPRO_rmin,EXPRO_n_exp)
  pflux_sum(1)           = ya
  pflux_sum(EXPRO_n_exp) = yb

  ! output omega_E, vtor_1 
  ! omega
  if(er_method /= 4) then
     do i=2,EXPRO_n_exp-1
        EXPRO_w0(i) = 2.9979e10*EXPRO_q(i)*(er_exp(i)/30.0)/ &
             ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
     enddo
     call bound_extrap(ya,yb,EXPRO_w0,EXPRO_rmin,EXPRO_n_exp)
     EXPRO_w0(1) = ya
     EXPRO_w0(EXPRO_n_exp) = yb
     ! omega_p
     call bound_deriv(EXPRO_w0p(2:EXPRO_n_exp-1),EXPRO_w0(2:EXPRO_n_exp-1),&
          EXPRO_rmin,EXPRO_n_exp-2)
     call bound_extrap(ya,yb,EXPRO_w0p,EXPRO_rmin,EXPRO_n_exp)
     EXPRO_w0p(1) = ya
     EXPRO_w0p(EXPRO_n_exp) = yb
  endif

  ! Write output on processor 0

  if (i_proc == 0) then

     open(unit=1,file='out.vgen.vel',status='replace')
     do i=1,EXPRO_n_exp
        write(1,'(e16.8)',advance='no') EXPRO_rho(i)
        write(1,'(e16.8)',advance='no') er_exp(i)
        write(1,'(e16.8)',advance='no') EXPRO_w0(i)
        write(1,'(e16.8)',advance='no') EXPRO_w0p(i)
        do j=1,n_ions
           write(1,'(e16.8)',advance='no') EXPRO_vpol(j,i)
           write(1,'(e16.8)',advance='no') EXPRO_vtor(j,i)
           write(1,'(e16.8)',advance='no') EXPRO_vpol(j,i) / EXPRO_bp0(i)
           write(1,'(e16.8)',advance='no') (EXPRO_vtor(j,i) &
                - EXPRO_vpol(j,i)*EXPRO_bt0(i)/ EXPRO_bp0(i)) &
                / (EXPRO_rmaj(i)+EXPRO_rmin(i))
        enddo
        write(1,*)
     enddo
     close(1)

     open(unit=1,file='input.profiles.jbs',status='replace')
     write(1,'(a)') '#'
     write(1,'(a)') '# expro_rho'
     write(1,'(a)') '# sum z*pflux_neo/(c_s n_e) /(rho_s/a_norm)**2'
     write(1,'(a)') '# jbs_neo    (MA/m^2)'
     write(1,'(a)') '# jbs_sauter (MA/m^2)'
     write(1,'(a)') '# jbs_nclass (MA/m^2)'
     write(1,'(a)') '# jbs_koh    (MA/m^2)'
     write(1,'(a)') '# where jbs = < j_parallel B > / B_unit'
     write(1,'(a)') '#'
     do i=1,EXPRO_n_exp
        write(1,'(6(1pe14.7,2x))') EXPRO_rho(i), pflux_sum(i), &
             jbs_neo(i), jbs_sauter(i), jbs_nclass(i), jbs_koh(i)
     enddo
     close(1)

     ! Write the new input.profiles
     vtag = '* This file regenerated by VGEN'
     call EXPRO_write_original(tag)

     ! Compute and write the new input.profiles.extra
     call EXPRO_compute_derived
     call EXPRO_write_derived

     call vgen_getgeo()

  endif

  deallocate(er_exp)
  deallocate(vtor_measured)
  deallocate(jbs_neo)
  deallocate(jbs_sauter)
  deallocate(jbs_koh)
  deallocate(jbs_nclass)
  deallocate(pflux_sum)

  call EXPRO_palloc(MPI_COMM_WORLD,path,0)

  call MPI_finalize(i_err)

10 format('rho=',f6.4,3x,'Er_0(kV/m)=',f10.4,3x,'vtor_1(km/s)=',f10.4,3x,'vpol_1(km/s)=',f10.4)

end program vgen
