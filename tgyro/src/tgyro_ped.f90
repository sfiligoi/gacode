!-------------------------------------------------------

module tgyro_ped

  ! EPED1NN inputs

  implicit none

  integer :: nr

  real :: a_in
  real :: betan_in
  real :: bt_in
  real :: delta_in
  real :: ip_in
  real :: kappa_in
  real :: m_in
  real :: neped_in
  real :: r_in
  real :: zeffped_in

  ! Typical values
  !
  !       a = 0.55
  !   betan = 1.28
  !      bt = 1.69
  !   delta = 0.54
  !      ip = 1.30
  !   kappa = 1.86
  !       m = 2.00
  !   neped = 3.62
  !       r = 1.70
  ! zeffped = 2.07

  real, dimension(:), allocatable :: rmin_eped
  real, dimension(:), allocatable :: polflux_eped
  real, dimension(:), allocatable :: polfluxp_eped

contains

  subroutine tgyro_pedestal

    use mpi
    use tgyro_globals

    implicit none

    real :: t_ped,n_ped
    real :: zt_ped
    real :: beta_ped_NN
    real :: p_ped_NN,w_ped_NN
    real :: p_top_NN,w_top_NN
    real :: t_ped_NN,n_ped_NN
    real :: n_edge_NN,t_edge_NN
    real :: psi_ped(1),r_ped(1),psip_ped(1)
    character(len=1000) :: nn_executable
    character(len=1000) :: nn_files
    integer dummy

    integer, parameter :: print_flag=1

    if (tgyro_ped_model == 1) return

    !-------------------------------------------------------------------------
    ! 1. Initializations
    !
    ! ** temporary **
    neped_in   = 3.62
    zeffped_in = 2.07
    betan_in   = 1.28

    ! All *_in variables set in tgyro_init_profiles

    if (print_flag == 1 .and. i_proc_global == 0) then
       print 10,'     a_in [m]',a_in
       print 10,' betan_in [-]',betan_in
       print 10,'    bt_in [T]',bt_in
       print 10,' delta_in [-]',delta_in
       print 10,'   ip_in [MA]',ip_in
       print 10,' kappa_in [-]',kappa_in
       print 10,'     m_in [-]',m_in
       print 10,' neped_in [-]',neped_in
       print 10,'     r_in [m]',r_in
       print 10,'  zeff_in [-]',zeffped_in
    endif
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 2. Call EPED1NN to get p_ped_NN, w_ped_NN
    !
    n_ped_NN  = neped_in
    n_edge_NN = neped_in*0.25
    t_edge_NN  = 80

    if (print_flag == 1 .and. i_proc_global == 0) then
        !
        ! Write input file for the NN
        !
        open (unit=14, file="input.dat", action="write")
        write (14,*) '1'
        write (14,"(10(f6.3,x))") a_in       ,&
                                  betan_in   ,&
                                  bt_in      ,&
                                  delta_in   ,&
                                  ip_in      ,&
                                  kappa_in   ,&
                                  m_in       ,&
                                  neped_in   ,&
                                  r_in       ,&
                                  zeffped_in
         close(14)

        !
        ! Execute the NN
        !
        call get_environment_variable('BRAINFUSE_RUN',nn_executable)
        call get_environment_variable('EPED1NN',nn_files)
        call system(trim(nn_executable)//' '//trim(nn_files)//' input.dat')
    endif

    ! All processes wait for pedestal calculation to finish and read results
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    !
    ! Read outputs
    !
    open (unit=14, file="output.avg", action="read")
    read(14,*) dummy
    read(14,*) beta_ped_NN,p_ped_NN,p_top_NN,w_ped_NN,w_top_NN
    close(14)

    t_ped_NN=p_ped_NN/n_ped_NN/1.6021766208/2.*1E6 ! From pressure in MPa to temperature in eV

    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 3. Map EPED1NN outputs to TGYRO pedestal parameters
    !
    nr = size(rmin_eped)
    psi_ped(1) = 1-2.0*w_ped_NN
    call cub_spline(polflux_eped,rmin_eped,nr,psi_ped,r_ped,1)
    call cub_spline(polflux_eped,polfluxp_eped,nr,psi_ped,psip_ped,1)

    t_ped  = t_ped_NN*1e-3 ! Convert from eV to keV
    n_ped  = n_ped_NN

    ! z = -1/T dT/dr at r=r_ped
    ! Use dT/dr = dT/dx dx/dr where x = Psi_norm 
    zt_ped = (t_ped_NN-t_edge_NN)/(2*tanh(1.0))*(1/w_ped_NN)*(1-tanh(1.0)**2) &
         *psip_ped(1)/t_ped_NN
    !-------------------------------------------------------------------------

    if (print_flag == 1 .and. i_proc_global == 0) then
       print *
       print 10,'   t_ped [keV]',t_ped
       print 10,' zt_ped [1/cm]',zt_ped
       print 10,'dlntedr [1/cm]',dlnnedr(n_r)
    endif

    !-------------------------------------------------------------------------
    ! 4. Integrate to obtain TGYRO pivot

    ! Integration backward from r_ped to r_star
    ti(:,n_r) = t_ped*exp(0.5*(dlntidr(:,n_r)+zt_ped)*(r_ped(1)-r(n_r)))
    te(n_r)   = t_ped*exp(0.5*(dlntedr(n_r)+zt_ped)*(r_ped(1)-r(n_r)))
    !-------------------------------------------------------------------------

    if (print_flag == 1 .and. i_proc_global == 0) then
       print *
       print 10,'psi_ped [-]',psi_ped(1)
       print 10,' r_ped [cm]',r_ped(1)
       print 10,'r(n_r) [cm]',r(n_r)
       print 10,' t_ped[keV]',t_ped
       print 10,'te_piv[keV]',te(n_r)
       stop
    endif

10  format(a,1pe12.5)

  end subroutine tgyro_pedestal

end module tgyro_ped
