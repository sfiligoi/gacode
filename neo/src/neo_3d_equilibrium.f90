module neo_3d_equilibrium

  implicit none

  public :: ThreeD_EQUIL_alloc, ThreeD_EQUIL_do
  
  ! equilibrium parameters (theta,varphi)
  ! (b dot grad)-- d/dth part
  real, dimension(:), allocatable   :: k_par_t_0
  real, dimension(:,:), allocatable :: k_par_t_1          
  ! (b dot grad)-- d/dphi part
  real, dimension(:), allocatable   :: k_par_p_0
  real, dimension(:,:), allocatable :: k_par_p_1
  ! (b cross grad B) dot grad r
  real, dimension(:), allocatable   :: v_drift_x_overB2_0
  real, dimension(:,:), allocatable :: v_drift_x_overB2_1 
  ! B 
  real, dimension(:), allocatable   :: Bmag_0
  real, dimension(:,:), allocatable :: Bmag_1
  ! (b dot grad B)/B
  real, dimension(:), allocatable   :: gradpar_Bmag_overB_0
  real, dimension(:,:), allocatable :: gradpar_Bmag_overB_1 
  ! flux surface avg weights
  real, dimension(:), allocatable   :: w_theta_0
  real, dimension(:,:), allocatable :: w_theta_1
  real :: sum_w_theta_0
  real :: sum_w_theta_1

  real :: d_theta, d_varphi
  
  ! private parameters
  logical, private :: initialized = .false.
  real, dimension(:), allocatable, private :: x
  real, private :: dx
  integer, parameter, private :: nx=7
  real, dimension(:), allocatable, private :: y
  real, private :: dy
  integer, parameter, private :: ny=7
  real, dimension(:), allocatable, private :: R0, Z0
  real, dimension(:,:,:), allocatable, private :: R1, Z1
  real, dimension(:), allocatable, private :: dR0dt, dZ0dt, dR0dr, dZ0dr, &
       dR0dtdt, dZ0dtdt, dR0drdt, dZ0drdt
  real, dimension(:,:,:), allocatable, private :: dZ1dt, dZ1dr, &
       dZ1dtdt, dZ1drdt, dZ1dp, dZ1drdp, dZ1dtdp, dZ1dpdp
  real, private :: c0
  real, dimension(:), allocatable, private :: Jfac, dJfacdt 
  integer, dimension(-2:2), private :: c2deriv
  integer, dimension(:), allocatable, private :: xcyc
  integer, dimension(:), allocatable, private :: ycyc

contains
  
  subroutine ThreeD_EQUIL_alloc(flag)
    use neo_globals, only: theta, n_theta, pi, n_varphi
    use neo_3d_globals, only: varphi
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it, ip
    
    if(flag == 1) then
       if(initialized) return
       
       allocate(theta(n_theta))
       allocate(varphi(n_varphi))

       allocate(k_par_t_0(n_theta))
       allocate(k_par_p_0(n_theta))
       allocate(v_drift_x_overB2_0(n_theta))
       allocate(Bmag_0(n_theta))
       allocate(gradpar_Bmag_overB_0(n_theta))
       allocate(w_theta_0(n_theta))

       allocate(k_par_t_1(n_theta,n_varphi))
       allocate(k_par_p_1(n_theta,n_varphi))
       allocate(v_drift_x_overB2_1(n_theta,n_varphi))
       allocate(Bmag_1(n_theta,n_varphi))
       allocate(gradpar_Bmag_overB_1(n_theta,n_varphi))
       allocate(w_theta_1(n_theta,n_varphi))
       
       d_theta = 2*pi/n_theta
       do it=1,n_theta
          theta(it) = -pi+(it-1)*d_theta
       enddo
       
       d_varphi = 2*pi/n_varphi
       do ip=1,n_varphi
          varphi(ip) = -pi+(ip-1)*d_varphi
       enddo
       
       initialized = .true.
       
    else
       if(.NOT. initialized) return
       
       deallocate(theta)
       deallocate(varphi)
       deallocate(k_par_t_0)
       deallocate(k_par_p_0)
       deallocate(v_drift_x_overB2_0)
       deallocate(Bmag_0)
       deallocate(gradpar_Bmag_overB_0)
       deallocate(w_theta_0)
       deallocate(k_par_t_1)
       deallocate(k_par_p_1)
       deallocate(v_drift_x_overB2_1)
       deallocate(Bmag_1)
       deallocate(gradpar_Bmag_overB_1)
       deallocate(w_theta_1)

       initialized = .false.
       
    endif
    
  end subroutine ThreeD_EQUIL_alloc
 
  subroutine ThreeD_EQUIL_do(ir)
    use neo_globals
    use neo_3d_globals
    implicit none
    integer, intent(in) :: ir
    integer :: it, ip, jt, jp, id, k
    !
    real, parameter :: zs=1.0
    real, parameter :: zc=0.0
    integer, parameter :: M_tor=1
    !
    real, dimension(:), allocatable :: g0,gtt0,gpp0
    real, dimension(:,:), allocatable :: g1,gtt1,gpp1,gtp1
    real, dimension(:,:,:), allocatable :: dR1dt, dR1dr, &
         dR1dtdt, dR1drdt, dR1dp, dR1drdp, dR1dtdp, dR1dpdp
    real :: db0dt, db1dp, db1dt
    real :: dg0dt, dgpp0dt, dgtt0dt
    real :: dg1dp, dg1dt, dgpp1dp, dgpp1dt, dgtt1dp, dgtt1dt, dgtp1dp, dgtp1dt
    !
    real, dimension(:,:), allocatable :: sum1, sum2
    real :: max_sum, s1, s2
    ! local equilibrium parameters (x,y)
  real, dimension(:), allocatable   :: k_par_t_0_loc, k_par_p_0_loc, v_drift_x_overB2_0_loc, &
       Bmag_0_loc, gradpar_Bmag_overB_0_loc, w_theta_0_loc
  real, dimension(:,:), allocatable :: k_par_t_1_loc, k_par_p_1_loc, v_drift_x_overB2_1_loc, &
       Bmag_1_loc, gradpar_Bmag_overB_1_loc, w_theta_1_loc    
  real, dimension(:,:), allocatable :: vec_xp

    ! theta grid for computing the equilibrium
    allocate(x(nx))
    dx = 2*pi/nx
    do it=1,nx
       x(it) = -pi+(it-1)*dx
    enddo
    
    allocate(xcyc(1-nx:2*nx))
    do it=1,nx
       xcyc(it-nx) = it
       xcyc(it) = it
       xcyc(it+nx) = it
    enddo

    ! varphi grid for computing the equilibrium
    allocate(y(ny))
    dy = 2*pi/ny
    do ip=1,ny
       y(ip) = -pi+(ip-1)*dy
    enddo
    
    allocate(ycyc(1-ny:2*ny))
    do ip=1,ny
       ycyc(ip-ny) = ip
       ycyc(ip) = ip
       ycyc(ip+ny) = ip
    enddo
    

    ! coefficients for 4th order centered 2nd-derivative
    c2deriv(-2) = -1
    c2deriv(-1) =  16
    c2deriv(0)  = -30
    c2deriv(1)  =  16
    c2deriv(2)  = -1

    ! specify the axisymmetric components
    allocate(R0(nx))
    allocate(Z0(nx))
    allocate(dR0dt(nx))
    allocate(dZ0dt(nx))
    allocate(dR0dr(nx))
    allocate(dZ0dr(nx))
    allocate(dR0dtdt(nx))
    allocate(dZ0dtdt(nx))
    allocate(dR0drdt(nx))
    allocate(dZ0drdt(nx))
    do it=1,nx
       R0(it)      = rmaj(ir) + r(ir)*cos(x(it))
       dR0dt(it)   = -r(ir)*sin(x(it))
       dR0dr(it)   = cos(x(it))
       dR0dtdt(it) = -r(ir)*cos(x(it))
       dR0drdt(it) = -sin(x(it))

       Z0(it)      = r(ir)*sin(x(it))
       dZ0dt(it)   = r(ir)*cos(x(it))
       dZ0dr(it)   = sin(x(it))
       dZ0dtdt(it) = -r(ir)*sin(x(it))
       dZ0drdt(it) = cos(x(it))
    enddo

    ! specify the non-axisymmetric components of Z = Z1
    ! Z1(:,:,1) = Z1(:,:,2)*cos(N*varphi) + Z1(:,:,3)*sin(N*varphi)
    allocate(Z1(nx,ny,3))
    allocate(dZ1dr(nx,ny,3))
    allocate(dZ1dt(nx,ny,3))
    allocate(dZ1dp(nx,ny,3))
    allocate(dZ1drdt(nx,ny,3))
    allocate(dZ1drdp(nx,ny,3))
    allocate(dZ1dtdt(nx,ny,3))
    allocate(dZ1dtdp(nx,ny,3))
    allocate(dZ1dpdp(nx,ny,3))
    do it=1,nx
       do ip=1,ny
          !
          Z1(it,ip,2)      = z1_mag * r(ir) * (zs * sin(M_tor*x(it))  + zc * cos(M_tor*x(it)))
          Z1(it,ip,3)      = z1_mag * r(ir) * (-zs * cos(M_tor*x(it)) + zc * sin(M_tor*x(it)))
          Z1(it,ip,1)      = Z1(it,ip,2) * cos(N_tor*y(ip)) + Z1(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1dr(it,ip,2)   = Z1(it,ip,2)/r(ir)
          dZ1dr(it,ip,3)   = Z1(it,ip,3)/r(ir)
          dZ1dr(it,ip,1)   = dZ1dr(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1dr(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1dp(it,ip,2)   =  N_tor * Z1(it,ip,3)
          dZ1dp(it,ip,3)   = -N_tor * Z1(it,ip,2)
          dZ1dp(it,ip,1)   = dZ1dp(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1dp(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1dt(it,ip,2)   = -M_tor * Z1(it,ip,3)
          dZ1dt(it,ip,3)   =  M_tor * Z1(it,ip,2)
          dZ1dt(it,ip,1)   = dZ1dt(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1dt(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1drdp(it,ip,2) = dZ1dp(it,ip,2)/r(ir)
          dZ1drdp(it,ip,3) = dZ1dp(it,ip,3)/r(ir)
          dZ1drdp(it,ip,1) = dZ1drdp(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1drdp(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1drdt(it,ip,2) = dZ1dt(it,ip,2)/r(ir)
          dZ1drdt(it,ip,3) = dZ1dt(it,ip,3)/r(ir)
          dZ1drdt(it,ip,1) = dZ1drdt(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1drdt(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1dpdp(it,ip,2) = -N_tor**2 * Z1(it,ip,2)
          dZ1dpdp(it,ip,3) = -N_tor**2 * Z1(it,ip,3)
          dZ1dpdp(it,ip,1) = dZ1dpdp(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1dpdp(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1dtdt(it,ip,2) = -M_tor**2 * Z1(it,ip,2)
          dZ1dtdt(it,ip,3) = -M_tor**2 * Z1(it,ip,3)
          dZ1dtdt(it,ip,1) = dZ1dtdt(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1dtdt(it,ip,3) * sin(N_tor*y(ip))
          !
          dZ1dtdp(it,ip,2) = N_tor * M_tor * Z1(it,ip,2)
          dZ1dtdp(it,ip,3) = N_tor * M_tor * Z1(it,ip,3)
          dZ1dtdp(it,ip,1) = dZ1dtdp(it,ip,2) * cos(N_tor*y(ip)) &
               + dZ1dtdp(it,ip,3) * sin(N_tor*y(ip))
       enddo
    enddo

    allocate(Jfac(nx))
    allocate(dJfacdt(nx))
    do it=1,nx
       Jfac(it) = dR0dr(it)*dZ0dt(it)-dR0dt(it)*dZ0dr(it)
       dJfacdt(it) = dR0drdt(it)*dZ0dt(it) + dR0dr(it)*dZ0dtdt(it) &
            - dR0dtdt(it)*dZ0dr(it) - dR0dt(it)*dZ0drdt(it)
    enddo

    c0 = 0.0
    do it=1,nx
       c0 = c0 + Jfac(it) / R0(it) /nx
    enddo

    ! compute the non-axisymmetric components of R = R1
    allocate(R1(nx,ny,3))
    call threed_solve(ir)

    allocate(dR1dr(nx,ny,3))
    allocate(dR1dt(nx,ny,3))
    allocate(dR1dp(nx,ny,3))
    allocate(dR1drdt(nx,ny,3))
    allocate(dR1drdp(nx,ny,3))
    allocate(dR1dtdt(nx,ny,3))
    allocate(dR1dtdp(nx,ny,3))
    allocate(dR1dpdp(nx,ny,3))
    dR1dr(:,:,:)   = 0.0
    dR1dt(:,:,:)   = 0.0
    dR1dp(:,:,:)   = 0.0
    dR1drdt(:,:,:) = 0.0
    dR1drdp(:,:,:) = 0.0
    dR1dtdt(:,:,:) = 0.0
    dR1dtdp(:,:,:) = 0.0
    dR1dpdp(:,:,:) = 0.0
    do it=1,nx
       do ip=1,ny
          !
          dR1dr(it,ip,2) = R1(it,ip,2)/r(ir) 
          dR1dr(it,ip,3) = R1(it,ip,3)/r(ir)
          dR1dr(it,ip,1) = dR1dr(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1dr(it,ip,3) * sin(N_tor * y(ip))
          !
          dR1dp(it,ip,2) = R1(it,ip,3) * N_tor
          dR1dp(it,ip,3) = R1(it,ip,2) * (-N_tor)
          dR1dp(it,ip,1) = dR1dp(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1dp(it,ip,3) * sin(N_tor * y(ip))
          !
          dR1dpdp(it,ip,2) = R1(it,ip,2) * (-N_tor**2)
          dR1dpdp(it,ip,3) = R1(it,ip,3) * (-N_tor**2)
          dR1dpdp(it,ip,1) = dR1dpdp(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1dpdp(it,ip,3) * sin(N_tor * y(ip))
          !
          dR1drdp(it,ip,2) = R1(it,ip,3)/r(ir) * N_tor
          dR1drdp(it,ip,3) = R1(it,ip,2)/r(ir) * (-N_tor)
          dR1drdp(it,ip,1) = dR1drdp(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1drdp(it,ip,3) * sin(N_tor * y(ip))
          !
          do id=-2,2
             jt = xcyc(it+id)
             dR1dt(it,ip,2)   = dR1dt(it,ip,2)   + cderiv(id)  / (12.0*dx) * R1(jt,ip,2)
             dR1dt(it,ip,3)   = dR1dt(it,ip,3)   + cderiv(id)  / (12.0*dx) * R1(jt,ip,3)
             dR1dtdt(it,ip,2) = dR1dtdt(it,ip,2) + c2deriv(id) / (12.0*dx*dx) * R1(jt,ip,2)
             dR1dtdt(it,ip,3) = dR1dtdt(it,ip,3) + c2deriv(id) / (12.0*dx*dx) * R1(jt,ip,3)
          enddo
          dR1dt(it,ip,1)   = dR1dt(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1dt(it,ip,3) * sin(N_tor * y(ip))
          dR1dtdt(it,ip,1) = dR1dtdt(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1dtdt(it,ip,3) * sin(N_tor * y(ip))
          !
          dR1drdt(it,ip,2) = dR1dt(it,ip,2)/r(ir)
          dR1drdt(it,ip,3) = dR1dt(it,ip,3)/r(ir)
          dR1drdt(it,ip,1) = dR1drdt(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1drdt(it,ip,3) * sin(N_tor * y(ip))
       enddo
    enddo
    do it=1,nx
       do ip=1,ny
          do id=-2,2
             jt = xcyc(it+id)
             dR1dtdp(it,ip,2)   = dR1dtdp(it,ip,2) + cderiv(id)  / (12.0*dx) * dR1dp(jt,ip,2)
             dR1dtdp(it,ip,3)   = dR1dtdp(it,ip,3) + cderiv(id)  / (12.0*dx) * dR1dp(jt,ip,3)
          enddo
          dR1dtdp(it,ip,1) = dR1dtdp(it,ip,2) * cos(N_tor * y(ip)) &
               + dR1dtdp(it,ip,3) * sin(N_tor * y(ip))
       enddo
    enddo

    ! axisymmetric metrics
    allocate(g0(nx))
    allocate(gtt0(nx))
    allocate(gpp0(nx))
    do it=1,nx
       g0(it)   = c0 * R0(it)**2 / r(ir)
       gtt0(it) = (c0*R0(it)/Jfac(it))**2 * (dR0dt(it)**2 + dZ0dt(it)**2) 
       gpp0(it) = R0(it)**2
    enddo
    
    ! non-axisymmetric metrics
    allocate(g1(nx,ny))
    allocate(gpp1(nx,ny))
    allocate(gtt1(nx,ny))
    allocate(gtp1(nx,ny))
    do it=1,nx
       do ip=1,ny
          g1(it,ip) = c0 * R0(it) / r(ir) * R1(it,ip,1) &
               + c0 * R0(it)**2 / (r(ir) *Jfac(it)) * (dZ0dt(it)*dR1dr(it,ip,1) &
               - dZ0dr(it) * dR1dt(it,ip,1) + dR0dr(it)*dZ1dt(it,ip,1) &
               - dR0dt(it)*dZ1dr(it,ip,1))
          gpp1(it,ip) = 2 * R0(it) * R1(it,ip,1)
          gtt1(it,ip) = 2 * (c0*R0(it)/Jfac(it))**2 * (dR0dt(it) * dR1dt(it,ip,1) &
               + dZ0dt(it)*dZ1dt(it,ip,1))
          gtp1(it,ip) = c0*R0(it)/Jfac(it) * (dR0dt(it)*dR1dp(it,ip,1) &
               + dZ0dt(it)*dZ1dp(it,ip,1))
       enddo
    enddo

    ! Compute the equilibrium quantities on the local theta grid
    allocate(k_par_t_0_loc(nx))
    allocate(k_par_p_0_loc(nx))
    allocate(v_drift_x_overB2_0_loc(nx))
    allocate(Bmag_0_loc(nx))
    allocate(gradpar_Bmag_overB_0_loc(nx))
    allocate(w_theta_0_loc(nx))
    allocate(k_par_t_1_loc(nx,ny))
    allocate(k_par_p_1_loc(nx,ny))
    allocate(v_drift_x_overB2_1_loc(nx,ny))
    allocate(Bmag_1_loc(nx,ny))
    allocate(gradpar_Bmag_overB_1_loc(nx,ny))
    allocate(w_theta_1_loc(nx,ny))

    sum_w_theta_0 = 0.0
    sum_w_theta_1 = 0.0
    do it=1,nx

       ! B
       Bmag_0_loc(it) = 1/g0(it) * sqrt(gpp0(it) + gtt0(it)/q(ir)**2)

       ! bhat dot grad
       k_par_p_0_loc(it) = 1.0/(Bmag_0_loc(it)*g0(it))

       k_par_t_0_loc(it) = k_par_p_0_loc(it) * c0 * R0(it) / Jfac(it) / q(ir)

       dg0dt   = 2.0 * c0 * dR0dt(it) * R0(it) / r(ir)
       dgtt0dt = 2.0 * (c0*R0(it)/Jfac(it))**2 &
            * (dR0dt(it) * dR0dtdt(it) + dZ0dt(it) * dZ0dtdt(it) &
            + (1/R0(it) * dR0dt(it) - 1/Jfac(it) * dJfacdt(it)) &
            * (dR0dt(it)**2 + dZ0dt(it)**2) )
       dgpp0dt = 2.0 * dR0dt(it) * R0(it)
       
       db0dt = -1/g0(it)**2 * dg0dt * sqrt(gpp0(it) + gtt0(it)/q(ir)**2) &
            + 1/g0(it) * 0.5 / sqrt(gpp0(it) + gtt0(it)/q(ir)**2) &
            * (dgpp0dt + dgtt0dt / q(ir)**2)

       ! bhat dot grad B / B
       gradpar_Bmag_overB_0_loc(it) = 1.0/(Bmag_0_loc(it)*g0(it)) &
            *c0/q(ir)*R0(it)/Jfac(it)*db0dt * (1.0/Bmag_0_loc(it))
       
       ! bhat X grad B dot grad r / B^2
       v_drift_x_overB2_0_loc(it) = -rho(ir)/(r(ir) &
            * Bmag_0_loc(it)**3 * g0(it)**2) &
            * gpp0(it)*c0/Jfac(it)*R0(it)*db0dt
       
       ! flux-surf avg weights
       w_theta_0_loc(it) =  Jfac(it) / (c0*R0(it)) * g0(it)

       do ip=1,ny

          ! B
          Bmag_1_loc(it,ip) = -g1(it,ip)/g0(it) * Bmag_0_loc(it) &
               + 0.5 / g0(it)**2 / Bmag_0_loc(it) &
               * (gpp1(it,ip) + gtt1(it,ip)/q(ir)**2 + gtp1(it,ip)/q(ir))

          ! bhat dot grad
          k_par_p_1_loc(it,ip) = -(Bmag_1_loc(it,ip)*g0(it) &
               + Bmag_0_loc(it)*g1(it,ip)) / (Bmag_0_loc(it)*g0(it))**2

          k_par_t_1_loc(it,ip) = k_par_p_1_loc(it,ip) * c0 * R0(it) &
               / Jfac(it) / q(ir)

          
          dg1dp = c0/r(ir) * (R0(it)*dR1dp(it,ip,1) &
               + 1/Jfac(it)*R0(it)**2 * (dR1drdp(it,ip,1) * dZ0dt(it) &
               + dZ1dtdp(it,ip,1) * dR0dr(it) - dZ1drdp(it,ip,1) * dR0dt(it) &
               - dR1dtdp(it,ip,1) * dZ0dr(it)))
          
          dg1dt = c0/r(ir) * (R1(it,ip,1) * dR0dt(it) + R0(it)*dR1dt(it,ip,1) &
               + (2.0/Jfac(it) * R0(it) * dR0dt(it) &
               - dJfacdt(it)/Jfac(it)**2 *R0(it)**2) &
               * (dR1dr(it,ip,1)*dZ0dt(it) + dZ1dt(it,ip,1)*dR0dr(it) &
               - dZ1dr(it,ip,1)*dR0dt(it) - dR1dt(it,ip,1)*dZ0dr(it)) &
               + 1/Jfac(it)*R0(it)**2 * (dR1dr(it,ip,1)*dZ0dtdt(it) &
               + dR1drdt(it,ip,1)*dZ0dt(it) + dZ1dtdt(it,ip,1)*dR0dr(it) &
               + dZ1dt(it,ip,1)*dR0drdt(it) - dZ1drdt(it,ip,1)*dR0dt(it) &
               - dZ1dr(it,ip,1)*dR0dtdt(it) - dR1dtdt(it,ip,1)*dZ0dr(it) &
               - dR1dt(it,ip,1)*dZ0drdt(it)))

          dgpp1dp = 2.0 * R0(it) * dR1dp(it,ip,1)
          dgpp1dt = 2.0 * (R0(it) * dR1dt(it,ip,1) + R1(it,ip,1) * dR0dt(it))

          dgtt1dp  = 2.0 * (c0*R0(it)/Jfac(it))**2 &
               * (dR0dt(it)*dR1dtdp(it,ip,1) + dZ0dt(it)*dZ1dtdp(it,ip,1))

          dgtt1dt = 2.0 * (c0*R0(it)/Jfac(it))**2 &
               * (dR0dt(it)*dR1dtdt(it,ip,1) + dR0dtdt(it)*dR1dt(it,ip,1) &
               + dZ0dt(it)*dZ1dtdt(it,ip,1) + dZ0dtdt(it)*dZ1dt(it,ip,1)) &
               + 4.0 * (c0*R0(it)/Jfac(it))**2 * (1/R0(it)*dR0dt(it) &
               - 1/Jfac(it)*dJfacdt(it)) * (dR0dt(it)*dR1dt(it,ip,1) &
               + dZ0dt(it)*dZ1dt(it,ip,1))

          dgtp1dp = (c0*R0(it)/Jfac(it)) * (dR0dt(it)*dR1dpdp(it,ip,1) &
               + dZ0dt(it) * dZ1dpdp(it,ip,1))

          dgtp1dt = (c0/Jfac(it) * dR0dt(it) &
               - dJfacdt(it)/Jfac(it)**2 * R0(it)) &
               * (dR0dt(it)*dR1dp(it,ip,1) + dZ0dt(it)*dZ1dp(it,ip,1)) &
               + (c0*R0(it)/Jfac(it)) * (dR0dtdt(it)*dR1dp(it,ip,1) &
               + dR0dt(it)*dR1dtdp(it,ip,1) + dZ0dtdt(it)*dZ1dp(it,ip,1) &
               + dZ0dt(it)*dZ1dtdp(it,ip,1))

          db1dp = -dg1dp * Bmag_0_loc(it)/ g0(it) &
               + 0.5/(g0(it)**2 * Bmag_0_loc(it)) &
               * (dgpp1dp * 1.0/q(ir)**2 * dgtt1dp + 2.0/q(ir) * dgtp1dp)

          db1dt = -dg1dt * Bmag_0_loc(it) / g0(it) &
               + g1(it,ip)/g0(it)**2 * dg1dt * Bmag_0_loc(it) &
               - g1(it,ip)/g0(it) * db0dt + 0.5/(g0(it)**2 * Bmag_0_loc(it)) &
               * (dgpp1dt * 1.0/q(ir)**2 * dgtt1dt + 2.0/q(ir) * dgtp1dt) &
               - (dg0dt / (g0(it)**3 * Bmag_0_loc(it)) &
               + db0dt/(g0(it)**2 * Bmag_0_loc(it)**2)) &
               * (gpp1(it,ip) + gtt1(it,ip)/q(ir)**2 + 2.0/q(ir) * gtp1(it,ip))


          ! bhat dot grad B / B
          gradpar_Bmag_overB_1_loc(it,ip) = 1.0/(Bmag_0_loc(it)*g0(it)) &
               *c0/q(ir)*R0(it)/Jfac(it)*db0dt &
               * (- Bmag_1_loc(it,ip)/Bmag_0_loc(it)**2) &
               + 1.0/(Bmag_0_loc(it)**2 *g0(it)) &
               * (c0/q(ir)*R0(it)/Jfac(it)*db1dt + db1dp) &
               - 1.0/Bmag_0_loc(it) * 1.0/(Bmag_0_loc(it)*g0(it))**2 &
               * (Bmag_1_loc(it,ip)*g0(it) + Bmag_0_loc(it)*g1(it,ip)) &
               * c0/q(ir)*R0(it)/Jfac(it)*db0dt

          ! bhat X grad B dot grad r / B^2
          v_drift_x_overB2_1_loc(it,ip) = rho(ir) &
               /(r(ir) * Bmag_0_loc(it)**3 * g0(it)**2) &
               * (-gpp0(it)*c0/Jfac(it)*R0(it)*db0dt &
               * (- 3.0*Bmag_1_loc(it,ip)/Bmag_0_loc(it) &
               - 2.0*g1(it,ip)/g0(it)) &
               - (gpp1(it,ip) + gtp1(it,ip)/q(ir)) * c0/Jfac(it)*R0(it)*db0dt &
               - gpp0(it) * c0/Jfac(it)*R0(it)*db1dt &
               + gtt0(it)/q(ir)*db1dp)

          ! flux-surface avg integration weights
          w_theta_1_loc(it,ip) =  Jfac(it) / (c0*R0(it)) * g1(it,ip)

       enddo
    enddo

    ! map the equilibrium paramters from the local theta grid to the 
    ! computational grid
    call cub_spline(x,k_par_t_0_loc,nx,theta,k_par_t_0,n_theta)
    call cub_spline(x,k_par_p_0_loc,nx,theta,k_par_p_0_loc,n_theta)
    call cub_spline(x,Bmag_0_loc,nx,theta,Bmag_0_loc,n_theta)
    call cub_spline(x,gradpar_Bmag_overB_0_loc,nx,theta,gradpar_Bmag_overB_0,n_theta)
    call cub_spline(x,w_theta_0_loc,nx,theta,w_theta_0,n_theta)
    call cub_spline(x,v_drift_x_overB2_0_loc,nx,theta,v_drift_x_overB2_0,n_theta)
    allocate(vec_xp(nx,n_varphi))
    vec_xp(:,:) = 0.0
    do it=1,nx
       call cub_spline(y,k_par_t_1_loc(it,:),ny,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(x,vec_xp(:,ip),nx,theta,k_par_t_1(:,ip),n_theta)
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nx
       call cub_spline(y,k_par_p_1_loc(it,:),ny,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(x,vec_xp(:,ip),nx,theta,k_par_p_1(:,ip),n_theta)
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nx
       call cub_spline(y,Bmag_1_loc(it,:),ny,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(x,vec_xp(:,ip),nx,theta,Bmag_1(:,ip),n_theta)
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nx
       call cub_spline(y,gradpar_Bmag_overB_1_loc(it,:),ny,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(x,vec_xp(:,ip),nx,theta,gradpar_Bmag_overB_1(:,ip),n_theta)
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nx
       call cub_spline(y,w_theta_1_loc(it,:),ny,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(x,vec_xp(:,ip),nx,theta,w_theta_1(:,ip),n_theta)
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nx
       call cub_spline(y,v_drift_x_overB2_1_loc(it,:),ny,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(x,vec_xp(:,ip),nx,theta,v_drift_x_overB2_1(:,ip),n_theta)
    enddo
    deallocate(vec_xp)

    do it=1,n_theta
       sum_w_theta_0 = sum_w_theta_0 + w_theta_0(it)
       do ip=1,n_varphi
          sum_w_theta_1 = sum_w_theta_1 + w_theta_1(it,ip)
       enddo
    enddo

    ! EAB check
    allocate(sum1(nx,ny))
    allocate(sum2(nx,ny))
    sum1(:,:) = 0.0
    sum2(:,:) = 0.0
    do it=1,nx
       do ip=1,ny
          do id=-2,2
             jt = xcyc(it+id)
             jp = ycyc(ip+id)
             sum1(it,ip) = sum1(it,ip) &
                  + cderiv(id) / (12.0*dx) &
                  * c0*R0(it)/Jfac(it) &
                  * ( (gpp1(jt,ip)+gtp1(jt,ip)/q(ir))/g0(jt) &
                  -   (gpp0(jt))*g1(jt,ip)/g0(jt)**2)
             sum2(it,ip) = sum2(it,ip) &
                  + cderiv(id) / (12.0*dy) &
                  * ( (gtp1(it,jp)+gtt1(it,jp)/q(ir))/g0(it) &
                  - gtt0(it)/q(ir)*g1(it,jp)/g0(it)**2 )
          enddo
       enddo
    enddo 
    max_sum = 1.0e-10
    s1 = 0.0
    s2 = 0.0
    do it=1,nx
       do ip=1,ny
          if(abs(sum1(it,ip)-sum2(it,ip)) > abs(max_sum)) then
             max_sum = sum1(it,ip)-sum2(it,ip)
             s1 = sum1(it,ip)
             s2 = sum2(it,ip)
          endif
          !print *, sum1(it,ip), sum2(it,ip)
       enddo
    enddo
    print *, 'g1 eqn max_sum = ', max_sum
    print *, s1, s2       

    sum1(:,:) = 0.0
    sum2(:,:) = 0.0
    do it=1,nx
       do ip=1,ny
          !sum1(it,ip) = sum1(it,ip) &
          !     - 2.0*dR0dt(it)*(gpp1(it,ip) + gtp1(it,ip)/q(ir)) &
          !     + 2.0 * dR0dt(it) * g1(it,ip)/(c0/r(ir)) 
          sum1(it,ip) = sum1(it,ip) &
               + dR1dtdt(it,ip,1)
          do id=-2,2
             jt = xcyc(it+id)
             jp = ycyc(ip+id)
             !sum1(it,ip) = sum1(it,ip) &
             !     + cderiv(id) / (12.0*dx) &
             !     * c0*R0(it)/Jfac(it) / (g0(jt) + g1(jt,ip)) &
             !     * (gpp0(jt)+gpp1(jt,ip)+gtp1(jt,ip)/q(ir))
             !sum2(it,ip) = sum2(it,ip) &
             !     + cderiv(id) / (12.0*dy) &
             !     * 1.0 / (g0(it) + g1(it,jp)) &
             !     * (gtp1(it,jp)+gtt0(it)/q(ir)+gtt1(it,jp)/q(ir))
             !
             !sum1(it,ip) = sum1(it,ip) &
             !     + R0(it)*cderiv(id) / (12.0*dx)*(gpp1(jt,ip) + gtp1(jt,ip)/q(ir)) &
             !     - R0(it)/(c0/r(ir))*cderiv(id) / (12.0*dx)*g1(jt,ip)
             !sum2(it,ip) = sum2(it,ip) &
             !     + Jfac(it)/c0 * cderiv(id) / (12.0*dy)*(gtp1(it,jp)+gtt1(it,jp)/q(ir)) &
             !     - c0/Jfac(it)/q(ir) * (dR0dt(it)**2 + dZ0dt(it)**2) / (c0/r(ir)) &
             !     * cderiv(id) / (12.0*dy) * g1(it,jp)
             sum2(it,ip) = sum2(it,ip) &
                  + cderiv(id) / (12.0*dx)*dR1dt(jt,ip,1)
                 
          enddo
       enddo
    enddo
    max_sum = 1.0e-10
    s1 = 0.0
    s2 = 0.0
    do it=1,nx
       do ip=1,1
          if(abs(sum1(it,ip)-sum2(it,ip)) > abs(max_sum)) then
             max_sum = sum1(it,ip)-sum2(it,ip)
             s1 = sum1(it,ip)
             s2 = sum2(it,ip)
          endif
          !print *, x(it),sum1(it,ip), sum2(it,ip)
       enddo
    enddo
    print *, 'g interm eqn max_sum = ', max_sum
    print *, s1, s2

    sum1(:,:) = 0.0
    sum2(:,:) = 0.0
    do it=1,nx
       do ip=1,ny
          k = 1
          sum1(it,ip) = dZ1dt(it,ip,k) * R0(it)**2 * &
               (-dJfacdt(it)*dR0dr(it)/Jfac(it)**2 + dR0drdt(it)/Jfac(it)) &
               + dZ1dr(it,ip,k) * R0(it)**2 * &
               (dJfacdt(it)*dR0dt(it)/Jfac(it)**2 - dR0dtdt(it)/Jfac(it)) &
               + dZ1dpdp(it,ip,k) * dZ0dt(it) &
               + dZ1dtdt(it,ip,k) * R0(it)**2 * dR0dr(it)/Jfac(it) &
               - dZ1drdt(it,ip,k) * R0(it)**2 * dR0dt(it)/Jfac(it) &
               - 1.0/q(ir)*dZ1dp(it,ip,k) * (-c0/Jfac(it)*dR0dt(it)*dZ0dt(it) &
               - 1/Jfac(it)**2 * dJfacdt(it) * R0(it)*dZ0dt(it) &
               + c0/Jfac(it)*R0(it)*dZ0dtdt(it)) &
               + 1.0/q(ir) * dZ1drdp(it,ip,k) * c0/Jfac(it)**2 * R0(it) &
               * dR0dt(it) &
               * (dR0dt(it)**2 + dZ0dt(it)**2) &
               - 1.0/q(ir) * dZ1dtdp(it,ip,k) * c0/Jfac(it) * R0(it) &
               * (-dZ0dt(it) + 1/Jfac(it)*dR0dr(it)*dR0dt(it)**2 &
               + 1/Jfac(it)*dR0dr(it)*dZ0dt(it)**2)

          !sum2(it,ip) = (-N_tor**2 * r(ir) * cos(x(it)) &
          !     + (1-M_tor**2)*R0(it)**2/r(ir)*cos(x(it))) &
          !     * z1_mag*r(ir) * (zs*sin(M_tor*x(it)) + zc*cos(M_tor*x(it))) &
          !     * cos(N_tor*y(ip)) &
          !     - 1.0/q(ir)*c0*r(ir)*N_tor*sin(x(it))*cos(x(it)) &
          !     * z1_mag*r(ir) * (-zs*cos(M_tor*x(it)) + zc*sin(M_tor*x(it))) &
          !     * cos(N_tor*y(ip)) &
          !     + (-N_tor**2 * r(ir) * cos(x(it)) &
          !     + (1-M_tor**2)*R0(it)**2/r(ir)*cos(x(it))) &
          !     * z1_mag*r(ir) * (-zs*cos(M_tor*x(it)) + zc*sin(M_tor*x(it))) &
          !     * sin(N_tor*y(ip)) &
          !     + 1.0/q(ir)*c0*r(ir)*N_tor*sin(x(it))*cos(x(it)) &
          !     * z1_mag*r(ir) * (zs*sin(M_tor*x(it)) + zc*cos(M_tor*x(it))) &
          !     * sin(N_tor*y(ip))

          sum2(it,ip) = sum2(it,ip) - dR0dt(it)*R1(it,ip,1) &
               + R0(it)*dR1dt(it,ip,1)*(1.0+R0(it)/r(ir)*dZ0drdt(it)) &
               - dR1dr(it,ip,1)*R0(it)*R0(it)/r(ir)*dZ0dtdt(it) &
               - dR1drdt(it,ip,1)*R0(it)*R0(it)/r(ir)*dZ0dt(it) &
               + dR1dtdt(it,ip,1)*R0(it)*R0(it)/r(ir)*dZ0dr(it) &
               + dR1dp(it,ip,1)*c0/q(ir)*(-dR0dt(it)**2/r(ir) &
               + R0(it)/r(ir)*dR0dtdt(it) + r(ir)) &
               - dR1dpdp(it,ip,1)*dR0dt(it) &
               - dR1dtdp(it,ip,1)*c0*R0(it)/q(ir)*(dR0dt(it)/r(ir)+dZ0dr(it)) &
               + dR1drdp(it,ip,1)*c0*R0(it)/q(ir)*dZ0dt(it)
    
       enddo
    enddo
    max_sum = 1.0e-10
    s1 = 0.0
    s2 = 0.0
    do it=1,nx
       do ip=1,ny
          if(abs(sum1(it,ip)-sum2(it,ip)) > abs(max_sum)) then
             max_sum = sum1(it,ip)-sum2(it,ip)
             s1 = sum1(it,ip)
             s2 = sum2(it,ip)
          endif
          if(ip==5) then
             print *, x(it),sum1(it,ip),sum2(it,ip)
          endif
       enddo
    enddo
    print *, 'R eqn max_sum = ', max_sum
    print *, s1, s2
    stop
    

    ! clean-up
    deallocate(k_par_t_0_loc)
    deallocate(k_par_p_0_loc)
    deallocate(v_drift_x_overB2_0_loc)
    deallocate(Bmag_0_loc)
    deallocate(gradpar_Bmag_overB_0_loc)
    deallocate(w_theta_0_loc)
    deallocate(k_par_t_1_loc)
    deallocate(k_par_p_1_loc)
    deallocate(v_drift_x_overB2_1_loc)
    deallocate(Bmag_1_loc)
    deallocate(gradpar_Bmag_overB_1_loc)
    deallocate(w_theta_1_loc)
    deallocate(R0)
    deallocate(Z0)
    deallocate(dR0dt)
    deallocate(dZ0dt)
    deallocate(dR0dr)
    deallocate(dZ0dr)
    deallocate(dR0dtdt)
    deallocate(dZ0dtdt)
    deallocate(dR0drdt)
    deallocate(dZ0drdt)
    deallocate(R1)
    deallocate(dR1dr)
    deallocate(dR1dt)
    deallocate(dR1dp)
    deallocate(dR1drdt)
    deallocate(dR1drdp)
    deallocate(dR1dtdt)
    deallocate(dR1dtdp)
    deallocate(dR1dpdp)
    deallocate(Z1)
    deallocate(dZ1dr)
    deallocate(dZ1dt)
    deallocate(dZ1dp)
    deallocate(dZ1drdt)
    deallocate(dZ1drdp)
    deallocate(dZ1dtdt)
    deallocate(dZ1dtdp)
    deallocate(dZ1dpdp)
    deallocate(Jfac)
    deallocate(dJfacdt)
    deallocate(g0)
    deallocate(gtt0)
    deallocate(gpp0)
    deallocate(g1)
    deallocate(gtt1)
    deallocate(gpp1)
    deallocate(gtp1)
    deallocate(x)
    deallocate(xcyc)
    deallocate(y)
    deallocate(ycyc)

    open(unit=1,file=trim(path)//'out.neo.btest_2',status='replace')
    do it=1, n_theta
       do ip=1,n_varphi
          write(1,'(e16.8)',advance='no') theta(it)
          write(1,'(e16.8)',advance='no') varphi(ip)
          write(1,'(e16.8)',advance='no') Bmag_1(it,ip)
          write(1,'(e16.8)',advance='no') v_drift_x_overB2_1(it,ip)
          write(1,'(e16.8)',advance='no') k_par_t_1(it,ip)
          write(1,'(e16.8)',advance='no') k_par_p_1(it,ip)
          write(1,'(e16.8)',advance='no') gradpar_Bmag_overB_1(it,ip)
          write(1,'(e16.8)',advance='no') w_theta_1(it,ip)/sum_w_theta_1
          write(1,*)
       enddo
    enddo
    close(1)

  end subroutine ThreeD_EQUIL_do

  subroutine threed_solve(ir)
    use neo_globals
    use neo_3d_globals
    implicit none
    integer, intent(in) :: ir
    integer :: it,jt,id,ip,k
    real, dimension(:,:), allocatable :: cmat
    real, dimension(:), allocatable   :: bvec, dvec
    real :: fac
    ! parameters for matrix solve
    integer :: msize
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    real, dimension(:), allocatable    :: work

    allocate(cmat(2*nx,2*nx))
    allocate(bvec(2*nx))

    cmat(:,:) = 0.0
    do it=1,nx
       jt=it
       fac = - dR0dt(it) - (R0(it)/r(ir))**2*dZ0dtdt(it)
       cmat(it,jt)       = cmat(it,jt) + fac
       cmat(it+nx,jt+nx) = cmat(it+nx,jt+nx) + fac
       !
       fac = N_tor*N_tor * dR0dt(it)
       cmat(it,jt)       = cmat(it,jt) + fac
       cmat(it+nx,jt+nx) = cmat(it+nx,jt+nx) + fac
       !
       fac = c0/q(ir)*dZ0dt(it)*R0(it)/r(ir)
       cmat(it,jt+nx)    = cmat(it,jt+nx) + fac*N_tor
       cmat(it+nx,jt)    = cmat(it+nx,jt) + fac*(-N_tor)
       !
       fac = c0/q(ir)/r(ir)*(-dR0dt(it)**2 + R0(it)*dR0dtdt(it) + r(ir)**2)
       cmat(it,jt+nx)    = cmat(it,jt+nx) + fac*N_tor
       cmat(it+nx,jt)    = cmat(it+nx,jt) + fac*(-N_tor)
       do id=-2,2
          jt = xcyc(it+id)
          fac = cderiv(id) / (12.0*dx) * R0(it) &
               * (1.0 + R0(it)/r(ir)*dZ0drdt(it)) &
               - cderiv(id) / (12.0*dx) * (R0(it)/r(ir))**2 * dZ0dt(it) &
               + c2deriv(id) / (12.0*dx*dx) * R0(it)**2/r(ir) * dZ0dr(it)
          cmat(it,jt)    = cmat(it,jt) + fac
          cmat(it+nx,jt+nx) = cmat(it+nx,jt+nx) + fac
       enddo
    enddo

    ! RHS
    bvec(:) = 0.0
    do it=1,nx
       do k=2,3
          if(k==2) then
             ! cos(N*x) part
             jt=0
          else
             ! sin(N*x) part
             jt=nx
          endif
          bvec(it+jt)    = dZ1dt(it,ip,k) * R0(it)**2 * &
               (-dJfacdt(it)*dR0dr(it)/Jfac(it)**2 + dR0drdt(it)/Jfac(it)) &
               + dZ1dr(it,ip,k) * R0(it)**2 * &
               (dJfacdt(it)*dR0dt(it)/Jfac(it)**2 - dR0dtdt(it)/Jfac(it)) &
               + dZ1dpdp(it,ip,k) * dZ0dt(it) &
               + dZ1dtdt(it,ip,k) * R0(it)**2 * dR0dr(it)/Jfac(it) &
               - dZ1drdt(it,ip,k) * R0(it)**2 * dR0dt(it)/Jfac(it) &
               - 1.0/q(ir)*dZ1dp(it,ip,k) * (-c0/Jfac(it)*dR0dt(it)*dZ0dt(it) &
               - 1/Jfac(it)**2 * dJfacdt(it) * R0(it)*dZ0dt(it) &
               + c0/Jfac(it)*R0(it)*dZ0dtdt(it)) &
               + 1.0/q(ir) * dZ1drdp(it,ip,k) * c0/Jfac(it)**2 * R0(it) * dR0dt(it) &
               * (dR0dt(it)**2 + dZ0dt(it)**2) &
               - 1.0/q(ir) * dZ1dtdp(it,ip,k) * c0/Jfac(it) * R0(it) &
               * (-dZ0dt(it) + 1/Jfac(it)*dR0dr(it)*dR0dt(it)**2 &
               + 1/Jfac(it)*dR0dr(it)*dZ0dt(it)**2)
       enddo
    enddo

    !do jt=1,nx
    !   print *, cmat(nx,jt), cmat(nx,jt+nx)
    !enddo
    !stop

    ! enforce the constraints
    cmat(nx,:)=0.0
    cmat(2*nx,:)=0.0
    do jt=1,nx
       cmat(nx,jt)      = 1.0
       cmat(2*nx,jt+nx) = 1.0
    enddo
    bvec(nx)   = 0.0
    bvec(2*nx) = 0.0

    ! matrix solve
    msize=2*nx
    allocate(work(msize))
    allocate(i_piv(msize))
    call DGETRF(msize,msize,cmat(:,:),msize,i_piv,info)

    call DGETRI(msize,cmat(:,:),msize,i_piv,work,msize,info)
    do it=1,2*nx
       print '(50(1pe12.5,1x))', cmat(it,:)
    enddo
    stop

    call DGETRS('N',msize,1,cmat(:,:),msize,i_piv,bvec(:),msize,info)
    deallocate(work)
    deallocate(i_piv)

    allocate(dvec(2*nx))
    dvec(:)=0.0
    do it=1,nx
       do id=-2,2
          jt = xcyc(it+id)
          dvec(it)    = dvec(it) + cderiv(id)  / (12.0*dx) * bvec(jt)
          dvec(it+nx) = dvec(it+nx) + cderiv(id)  / (12.0*dx) * bvec(jt+nx)
       enddo
    enddo
    open(unit=1,file=trim(path)//'out.neo.btest',status='replace')
    do it=1, nx
       write(1,'(e16.8)',advance='no') x(it)
       write(1,'(e16.8)',advance='no') bvec(it)
       write(1,'(e16.8)',advance='no') bvec(it+nx)
       write(1,'(e16.8)',advance='no') dvec(it)
       write(1,'(e16.8)',advance='no') dvec(it+nx)
       write(1,*)
    enddo
    close(1)
    deallocate(dvec)
    
    do it=1,nx
       do ip=1,ny
          R1(it,ip,2) = bvec(it)
          R1(it,ip,3) = bvec(it+nx)
          R1(it,ip,1) = bvec(it) * cos(N_tor * y(ip)) &
               + bvec(it+nx) * sin(N_tor * y(ip))
       enddo
    enddo


    ! clean-up
    deallocate(cmat)
    deallocate(bvec)

  end subroutine threed_solve

end module neo_3d_equilibrium
