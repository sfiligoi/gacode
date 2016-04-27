!-----------------------------------------------------------------
! tgyro_ped.f90
!
! PURPOSE:
!  Manage dynamic pedestal interface including infamous
!  no-man's-land (NML).
!-----------------------------------------------------------------

module tgyro_ped

  implicit none

  ! EPED_NN inputs
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

  ! EPED_NN outputs
  integer, parameter :: nx_nn=1001
  real :: nn_vec(nx_nn,3) 
  real :: nn_w_ped

  ! Derivatives wrt psi
  real :: t_vec(nx_nn)
  real :: n_p(nx_nn)
  real :: t_p(nx_nn)

  integer :: n_exp
  real, dimension(:), allocatable :: rmin_exp
  real, dimension(:), allocatable :: psi_exp
  real, dimension(:), allocatable :: dpsidr_exp

  ! Pedestal top scale lengths
  real :: zn_top,zt_top
  real :: r_top(1)

contains

  subroutine tgyro_pedestal

    use mpi
    use tgyro_globals

    implicit none

    ! Parameters interpolated at top of pedestal
    real, dimension(1) :: psi_top
    real, dimension(1) :: n_top,t_top,p_top
    real, dimension(1) :: n_p_top,t_p_top
    real, dimension(1) :: dpsidr_top

    integer, parameter :: print_flag=1

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

       ! Inputs stored in interface
       call tgyro_eped_nn

       psi_top(1) = 1.0-1.5*nn_w_ped

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
    t_vec = nn_vec(:,3)/(2*nn_vec(:,2)*k)
    call bound_deriv(t_p,t_vec,nn_vec(:,1),nx_nn)

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
       !stop
    endif

10  format(a,1pe12.5)

  end subroutine tgyro_pedestal

  subroutine tgyro_pedestal_map

    use tgyro_globals
    use EXPRO_interface

    implicit none

    integer :: i_exp,i0
    real :: r_exp(n_exp),z_exp(n_exp),zt_exp(n_exp)
    real :: ne_exp(n_exp)
    real :: x0

    r_exp(:) = 100.0*EXPRO_rmin(:)

    do i_exp=2,n_exp
       x0 = r_exp(i_exp)
       if (x0 > r(n_r) .and. x0 <= r_top(1)) then 
          z_exp(i_exp) = (dlnnedr(n_r)*(r_top(1)-x0)+zn_top*(x0-r(n_r)))/(r_top(1)-r(n_r))
          zt_exp(i_exp) = (dlntedr(n_r)*(r_top(1)-x0)+zt_top*(x0-r(n_r)))/(r_top(1)-r(n_r))
          !print *,i_exp,r_exp(i_exp),z_exp(i_exp),zt_exp(i_exp)
       endif
       if (x0 > r_top(1)) then
          i0 = i_exp
          exit
       endif
    enddo

    call cub_spline(nn_vec(:,1),-n_p/nn_vec(:,2),nx_nn,psi_exp(i0:n_exp),z_exp(i0:n_exp),n_exp-i0+1) 
    call cub_spline(nn_vec(:,1),-t_p/t_vec(:),nx_nn,psi_exp(i0:n_exp),zt_exp(i0:n_exp),n_exp-i0+1) 

    do i_exp=i0,n_exp
       z_exp(i_exp) = z_exp(i_exp)*dpsidr_exp(i_exp)
       zt_exp(i_exp) = zt_exp(i_exp)*dpsidr_exp(i_exp)
       !print *,i_exp,r_exp(i_exp),z_exp(i_exp)*dpsidr_exp(i_exp),zt_exp(i_exp)*dpsidr_exp(i_exp)
    enddo

    do i_exp=i0,1,-1
          ne_exp(i_exp-1) = ne_exp(i_exp)*exp(0.5*(z_exp(i_exp)+z_exp(i_exp-1))* &
               (r_exp(i_exp)-r_exp(i_exp-1)))
    enddo

  end subroutine tgyro_pedestal_map

end module tgyro_ped
