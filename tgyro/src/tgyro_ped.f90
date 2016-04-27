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

  integer :: n_exp
  real, dimension(:), allocatable :: rmin_exp
  real, dimension(:), allocatable :: psi_exp
  real, dimension(:), allocatable :: dpsidr_exp

contains

  subroutine tgyro_pedestal

    use mpi
    use tgyro_globals

    implicit none

    ! Parameters interpolated at top of pedestal
    real :: zn_top,zt_top
    real, dimension(1) :: psi_top
    real, dimension(1) :: n_top,t_top,p_top
    real, dimension(1) :: n_p_top,t_p_top
    real, dimension(1) :: r_top,dpsidr_top

    character(len=1000) :: nn_executable
    character(len=1000) :: nn_files
    character(len=1) :: dummy
    integer, parameter :: nx_nn=1001
    real :: nn_vec(nx_nn,3) 
    real :: n_p(nx_nn)
    real :: t_p(nx_nn)
    real :: w_ped

    integer, parameter :: print_flag=1
    integer, parameter :: test_flag=1

    if (tgyro_ped_model == 1) return

    !-------------------------------------------------------------------------
    ! 1. Initializations
    !
    ! ** temporary **
    neped_in   = 3.62
    zeffped_in = 2.07
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 2. Call EPED1NN 
    !
    if (i_proc_global == 0) then
       !
       ! Write input file for the NN
       open(unit=1,file='input.dat',status='replace')
       write(1,'(a)') '1'
       write(1,'(10(f6.3,1x))') &
            a_in       ,&
            betan_in   ,&
            bt_in      ,&
            delta_in   ,&
            ip_in      ,&
            kappa_in   ,&
            m_in       ,&
            neped_in   ,&
            r_in       ,&
            zeffped_in
       close(1)
       !
       ! Execute the NN
       if (test_flag == 0) then
          !call get_environment_variable('BRAINFUSE_RUN',nn_executable)
          !call get_environment_variable('EPED1NN',nn_files)
          !call execute_command_line(trim(nn_executable)//' '//trim(nn_files)//' input.dat')
       endif

       ! Read neuped_vec=[psi_norm,ne,ptot=2 ne T]
       open(unit=1,file='neuped.profiles',status='old')
       read(1,*) dummy
       read(1,*) nn_vec
       close(1)

       !open (unit=1, file="output.avg", action="read")
       !read(1,*) dummy
       !close(1)
       w_ped   = 0.05
       psi_top(1) = 1.0-1.5*w_ped

    endif

    ! Communicated needed data from output.avg
    call MPI_BCAST(nn_vec,size(nn_vec),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(psi_top,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! n_top: convert to 1/cm^3 from 1/m^3
    nn_vec(:,2) = nn_vec(:,2)*1e-6
    call cub_spline(nn_vec(:,1),nn_vec(:,2),nx_nn,psi_top,n_top,1)
    ! p_top: convert to Ba from Pa 
    nn_vec(:,3) = nn_vec(:,3)*10.0
    call cub_spline(nn_vec(:,1),nn_vec(:,3),nx_nn,psi_top,p_top,1) 

    ! FORMULA: P = 2nkT 

    ! t_top [eV]
    t_top = p_top/(2*n_top*k) 

    ! n', T':
    call bound_deriv(n_p,nn_vec(:,2),nn_vec(:,1),nx_nn)
    call bound_deriv(t_p,nn_vec(:,3)/(2*nn_vec(:,2)*k),nn_vec(:,1),nx_nn)

    ! n_top', t_top'
    call cub_spline(nn_vec(:,1),n_p,nx_nn,psi_top,n_p_top,1)
    call cub_spline(nn_vec(:,1),t_p,nx_nn,psi_top,t_p_top,1) 
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 3. Compute z parameters at pedestal top
    ! 
    ! zn_top = -(1/n) dn/dr 
    ! zt_top = -(1/T) dT/dr 
    !
    call cub_spline(psi_exp,rmin_exp    ,n_exp,psi_top,r_top,1)
    call cub_spline(psi_exp,dpsidr_exp,n_exp,psi_top,dpsidr_top,1)

    zn_top = -n_p_top(1)/n_top(1)*dpsidr_top(1)
    zt_top = -t_p_top(1)/t_top(1)*dpsidr_top(1)
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 4. Integrate to obtain TGYRO pivot

    ! Integration backward from r_top to r_star
    ti(:,n_r) = t_top(1)*exp(0.5*(dlntidr(:,n_r)+zt_top)*(r_top(1)-r(n_r)))
    te(n_r)   = t_top(1)*exp(0.5*(dlntedr(n_r)  +zt_top)*(r_top(1)-r(n_r)))
    ne(n_r)   = n_top(1)*exp(0.5*(dlnnedr(n_r)  +zn_top)*(r_top(1)-r(n_r)))
    !-------------------------------------------------------------------------

    if (print_flag == 1 .and. i_proc_global == 0) then
       print *
       print 10,'n_top [1/cm^3]',n_top
       print 10,'t_top     [eV]',t_top
       print 10,'p_top     [Ba]',p_top
       print 10,'zn_top  [1/cm]',zn_top
       print 10,'zt_top  [1/cm]',zt_top
       print 10,'dlntedr [1/cm]',dlnnedr(n_r)
       print *
       print 10,'psi_top [-]',psi_top(1)
       print 10,' r_top [cm]',r_top(1)
       print 10,'r(n_r) [cm]',r(n_r)
    endif

    stop

10  format(a,1pe12.5)

  end subroutine tgyro_pedestal

end module tgyro_ped
