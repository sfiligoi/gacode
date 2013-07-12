module neo_3d_equilibrium

  implicit none

  public :: ThreeD_EQUIL_alloc, ThreeD_EQUIL_do
  
  ! equilibrium parameters (theta,varphi)
  ! B 
  real, dimension(:,:), allocatable :: Bmag
  ! (bhat dot grad)-- d/dth part and d/dvarphi parts
  real, dimension(:,:), allocatable :: k_par_t
  real, dimension(:,:), allocatable :: k_par_p
  ! (bhat dot grad B)/B
  real, dimension(:,:), allocatable :: gradpar_Bmag_overB
  ! (bhat cross grad B) dot grad r
  real, dimension(:,:), allocatable :: v_drift_x_overB2
  ! flux surface avg weights
  real, dimension(:,:), allocatable :: w_theta
  real :: sum_w_theta

  real :: d_theta, d_varphi
  
  ! private parameters
  logical, private :: initialized = .false.

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

       allocate(Bmag(n_theta,n_varphi))
       allocate(k_par_t(n_theta,n_varphi))
       allocate(k_par_p(n_theta,n_varphi))
       allocate(gradpar_Bmag_overB(n_theta,n_varphi))
       allocate(v_drift_x_overB2(n_theta,n_varphi))
       allocate(w_theta(n_theta,n_varphi))
       
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
       deallocate(Bmag)
       deallocate(gradpar_Bmag_overB)
       deallocate(k_par_t)
       deallocate(k_par_p)
       deallocate(v_drift_x_overB2)
       deallocate(w_theta)

       initialized = .false.
       
    endif
    
  end subroutine ThreeD_EQUIL_alloc
 
  subroutine ThreeD_EQUIL_do(ir)
    use neo_globals
    use neo_3d_globals
    implicit none
    integer, intent(in) :: ir
    integer :: it,ip
    real, dimension(:,:), allocatable :: floc
    real, dimension(:,:), allocatable :: bloc, bloc_tderiv, bloc_pderiv

    ! Compute equil quantities on le3 straight field-line grid and 
    ! map to neo computational grid(ntheta,nvarphi)
    
    allocate(floc(nt,np))
    allocate(bloc(nt,np))
    allocate(bloc_tderiv(nt,np))
    allocate(bloc_pderiv(nt,np))
    allocate(vec_xp(nt,n_varphi))

    ! B
    do it=1,nt
       do ip=1,np
          bloc(it,ip) = 1.0/sqg(it,ip) &
               * sqrt(gpp(it,ip) + 2.0/q(ir) *gpt(it,ip) &
               + 1.0/q(ir)**2 * gtt(it,ip))
       enddo
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nt
       call cub_spline(p,bloc(it,:),np,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(t,vec_xp(:,ip),nt,theta,Bmag(:,ip),n_theta)
    enddo

    ! deriv of B on local grid
    tderiv = 0.0
    do id=-2,2
       k = tcyc(it+id)
       tderiv = tderiv + bloc(k,ip) * cderiv(id)/dt
    enddo
    pderiv = 0.0
    do id=-2,2
       k = pcyc(ip+id)
       pderiv = pderiv + bloc(it,k) * cderiv(id)/dp
    enddo

    ! (bhat dot grad)-- d/dth part
    do it=1,nt
       do ip=1,np
          floc(it,ip) = 1.0/sqg(it,ip)/bloc(it,ip) &
               * (1.0/q(ir) * dtbdt(it,ip) + dtbdp(it,ip))
       enddo
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nt
       call cub_spline(p,floc(it,:),np,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(t,vec_xp(:,ip),nt,theta,k_par_t(:,ip),n_theta)
    enddo

    ! (bhat dot grad)-- d/dvarphi part
    do it=1,nt
       do ip=1,np
          floc(it,ip) = 1.0/sqg(it,ip)/bloc(it,ip)
       enddo
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nt
       call cub_spline(p,floc(it,:),np,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(t,vec_xp(:,ip),nt,theta,k_par_p(:,ip),n_theta)
    enddo

    ! (bhat dot grad B)/B
    do it=1,nt
       do ip=1,np
          floc(it,ip) = 1.0/bloc(it,ip)**2 / sqg(it,ip) &
               * (1.0/q(ir) * bloc_tderiv(it,ip) + bloc_pderiv(it,ip))
       enddo
    enddo
    vec_xp(:,:) = 0.0
    do it=1,nt
       call cub_spline(p,floc(it,:),np,varphi,vec_xp(it,:),n_varphi)
    enddo
    do ip=1,n_varphi
       call cub_spline(t,vec_xp(:,ip),nt,theta,gradpar_Bmag_overB(:,ip),n_theta)
    enddo

    deallocate(floc)
    deallocate(bloc)
    deallocate(bloc_tderiv)
    deallocate(bloc_pderiv)
    deallocate(vec_xp)

  end subroutine ThreeD_EQUIL_do

end module neo_3d_equilibrium
