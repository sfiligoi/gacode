!-------------------------------------------------------

module tgyro_ped

  ! NEUPED inputs

  implicit none

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

contains

  subroutine tgyro_pedestal

    use tgyro_globals

    implicit none

    real :: r_star
    real :: n_star
    integer, parameter :: print_flag=1

    if (tgyro_ped_model == 1) return

    !-------------------------------------------------------------------------
    ! Initializations
    r_star = r(n_r)
    n_star = ne(n_r)

    ! All *_in variables set in tgyro_init_profiles

    if (print_flag == 1 .and. i_proc_global == 0) then
       print 10,'    a_in [m]',a_in
       print 10,'   bt_in [T]',bt_in
       print 10,'delta_in [-]',delta_in
       print 10,'kappa_in [-]',kappa_in
       print 10,'    m_in [-]',m_in
       print 10,'    r_in [m]',r_in
       stop
    endif

    ! 1. call neuped
    !call neuped

    ! 2. Map NEUPED outputs to TGYRO pedestal parameters
    !    [width_NP, t_ped_NP, n_ped_NP]
    !  r_ped  = psinorm_to_r(1.0-2*width_NP)
    !  t_ped  = t_ped_NP
    !  n_ped  = n_ped_NP
    !  zt_ped =  
    !  zn_ped = 

    !-------------------------------------------------------------------------
    ! 3. integrate to obtain TGYRO pivot

    ! Integration backward from r_ped to r_star
    !  t_star = t_ped*exp(0.5*(zt_star+zt_ped)*(r_ped-r_star))
    !  n_star = n_ped*exp(0.5*(zn_star+zn_ped)*(r_ped-r_star))

    ! Configuration variables for core TGYRO iteration
    !  ni(:,n_r) = n_star
    !  ti(:,n_r) = t_star  
    !  ne(n_r)   = n_star
    !  te(n_r)   = t_star  
    !-------------------------------------------------------------------------

10  format(a,1pe12.5)

  end subroutine tgyro_pedestal

end module tgyro_ped
