!-------------------------------------------------------

module tgyro_ped

  ! NEUPED inputs

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

    use tgyro_globals

    implicit none

    real :: r_star,n_star,t_star
    real :: t_ped,n_ped
    real :: zt_ped,zn_ped
    real :: zn_star,zt_star
    real :: width_NP,t_ped_NP,n_ped_NP,n_edge_NP,t_edge_NP
    real :: psi_ped(1),r_ped(1),psip_ped(1)

    integer, parameter :: print_flag=1

    if (tgyro_ped_model == 1) return

    !-------------------------------------------------------------------------
    ! 1. Initializations
    !
    r_star = r(n_r)
    n_star = ne(n_r)

    ! ** temporary **
    neped_in   = 3.62
    zeffped_in = 2.07

    ! All *_in variables set in tgyro_init_profiles

    if (print_flag == 1 .and. i_proc_global == 0) then
       print 10,'    a_in [m]',a_in
       print 10,'betan_in [-]',betan_in
       print 10,'   bt_in [T]',bt_in
       print 10,'delta_in [-]',delta_in
       print 10,'   ip_in [-]',ip_in
       print 10,'kappa_in [-]',kappa_in
       print 10,'    m_in [-]',m_in
       print 10,'neped_in [-]',neped_in
       print 10,'    r_in [m]',r_in
       print 10,'ze-ed_in [m]',zeffped_in
       stop
    endif
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 2. Call neuped
    !
    ! NEUPED OUTPUTS: width_NP, t_ped_NP, n_ped_NP
    width_NP = 0.1
    t_ped_NP = 1e3
    n_ped_NP = 1e6 
    n_edge_NP = 0.0
    t_edge_NP = 0.0
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 3. Map NEUPED outputs to TGYRO pedestal parameters
    !
    nr = size(rmin_eped)
    psi_ped(1) = 1-2.0*width_NP
    call cub_spline(polflux_eped,rmin_eped,nr,psi_ped,r_ped,1)
    call cub_spline(polflux_eped,polfluxp_eped,nr,psi_ped,psip_ped,1)

    t_ped  = t_ped_NP*1e-3 ! Convert from eV to keV
    n_ped  = n_ped_NP*1e-6 ! Convert from 10^13/m^3 to 10^19/m^3

    zt_ped = (t_ped-t_edge_NP)/(2*tanh(1.0))*(-1/width_NP)*(1-tanh(1.0)**2)*psip_ped(1)
    zn_ped = (n_ped-n_edge_NP)/(2*tanh(1.0))*(-1/width_NP)*(1-tanh(1.0)**2)*psip_ped(1)
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 4. Integrate to obtain TGYRO pivot
    
    ! Integration backward from r_ped to r_star
    t_star = t_ped*exp(0.5*(zt_star+zt_ped)*(r_ped(1)-r_star))
    n_star = n_ped*exp(0.5*(zn_star+zn_ped)*(r_ped(1)-r_star))
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 5. Map back to TGYRO
    !
    ni(:,n_r) = n_ped*exp(0.5*(dlnnidr(:,n_r)+zn_ped)*(r_ped(1)-r_star))
    ti(:,n_r) = t_ped*exp(0.5*(dlntidr(:,n_r)+zt_ped)*(r_ped(1)-r_star))
    ne(n_r)   = n_ped*exp(0.5*(dlnnedr(n_r)+zn_ped)*(r_ped(1)-r_star))
    te(n_r)   = t_ped*exp(0.5*(dlntedr(n_r)+zt_ped)*(r_ped(1)-r_star))
    !-------------------------------------------------------------------------

10  format(a,1pe12.5)

  end subroutine tgyro_pedestal

end module tgyro_ped
