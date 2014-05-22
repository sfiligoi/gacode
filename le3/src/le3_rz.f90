subroutine le3_rz(t0,p0,r0,z0,drdt0,drdp0,dzdt0,dzdp0,jac0,rc0)

  use le3_globals

  implicit none

  real, intent(in) :: t0,p0
  real, intent(inout) :: r0,z0,jac0
  real, intent(inout) :: drdt0,dzdt0
  real, intent(inout) :: drdp0,dzdp0
  real, intent(inout) :: rc0
  real :: drdr0,dzdr0
  real :: drdtt0,dzdtt0
  real :: ang,x,a,a_t,a_tt
  integer :: it, ip
  real, dimension(:), allocatable :: cost_geo, cosp_geo, sint_geo, sinp_geo

  select case(equilibrium_model)
     
  case (0)
     ! Miller 2D with simple model for 3D

     ang = m*t0+n*p0
     
     x    = asin(delta)
     a    = t0+x*sin(t0)
     a_t  = 1.0+x*cos(t0)
     a_tt = -x*sin(t0)

     ! R
     r0 = rmaj + rmin*cos(a) + hmin*cos(ang)
     
     ! dR/dr
     drdr0 = shift+cos(a)-s_delta/cos(x)*sin(t0)*sin(a) + dhmindr*cos(ang)
     
     ! dR/d(tb)
     drdt0 = -rmin*sin(a)*a_t - m*hmin*sin(ang) 
     
     ! dR/d(pb)
     drdp0 = -n*hmin*sin(ang)                 
     
     ! d^2R/d(tb)^2
     drdtt0 = -rmin*a_t**2*cos(a)-rmin*a_tt*sin(a)-m**2*hmin*cos(ang)
     
     a    = t0+zeta*sin(2*t0)
     a_t  = 1.0+2.0*zeta*cos(2*t0)
     a_tt = -4*zeta*sin(2*t0)
     
     ! Z
     z0 = zmag+kappa*rmin*sin(a) + hmin*sin(ang)
     
     ! dZ/dr
     dzdr0 = dzmag+kappa*(1.0+s_kappa)*sin(a)+kappa*s_zeta*cos(a)*sin(2*t0) &
          + dhmindr*sin(ang)
     
     ! dZ/d(tb)
     dzdt0 =  kappa*rmin*cos(a)*a_t + m*hmin*cos(ang) 
     
     ! dZ/d(pb)
     dzdp0 = n*hmin*cos(ang)
     
     ! d^2Z/d(tb)^2
     dzdtt0 = -kappa*rmin*sin(a)*a_t**2+kappa*rmin*cos(a)*a_tt&
          -m**2*hmin*sin(ang)
     
  case(1)
     ! numerical 2D + 3D as sum of fourier coeffs

     allocate(cost_geo(0:nts_geo))
     allocate(sint_geo(0:nts_geo))
     allocate(cosp_geo(0:nps_geo))
     allocate(sinp_geo(0:nps_geo))
     do it=0,nts_geo
        cost_geo(it) = cos(it*t0)
        sint_geo(it) = sin(it*t0)
     enddo
     do ip=0,nps_geo
        cosp_geo(ip) = cos(ip*p0)
        sinp_geo(ip) = sin(ip*p0)
     enddo
     

     r0    = 0.0
     z0    = 0.0
     drdr0 = 0.0
     dzdr0 = 0.0
     do it=0,nts_geo
        do ip=0,nps_geo
           r0 = r0 + r_geo(it,ip,1)*sint_geo(it)* cosp_geo(ip) &
                + r_geo(it,ip,2)*sint_geo(it)* sinp_geo(ip) &
                + r_geo(it,ip,3)*cost_geo(it)* cosp_geo(ip) &
                + r_geo(it,ip,4)*cost_geo(it)* sinp_geo(ip) 

           z0 = z0 + z_geo(it,ip,1)*sint_geo(it)* cosp_geo(ip) &
                + z_geo(it,ip,2)*sint_geo(it)* sinp_geo(ip) &
                + z_geo(it,ip,3)*cost_geo(it)* cosp_geo(ip) &
                + z_geo(it,ip,4)*cost_geo(it)* sinp_geo(ip) 

           drdr0 = drdr0 + rd_geo(it,ip,1)*sint_geo(it)* cosp_geo(ip) &
                + rd_geo(it,ip,2)*sint_geo(it)* sinp_geo(ip) &
                + rd_geo(it,ip,3)*cost_geo(it)* cosp_geo(ip) &
                + rd_geo(it,ip,4)*cost_geo(it)* sinp_geo(ip) 

           dzdr0 = dzdr0 + zd_geo(it,ip,1)*sint_geo(it)* cosp_geo(ip) &
                + zd_geo(it,ip,2)*sint_geo(it)* sinp_geo(ip) &
                + zd_geo(it,ip,3)*cost_geo(it)* cosp_geo(ip) &
                + zd_geo(it,ip,4)*cost_geo(it)* sinp_geo(ip) 
        enddo
     enddo

     drdt0  = 0.0
     dzdt0  = 0.0
     do it=0,nts_geo
        do ip=0,nps_geo
           drdt0 = drdt0 + r_geo(it,ip,1)*it*cost_geo(it)* cosp_geo(ip) &
                + r_geo(it,ip,2)*it*cost_geo(it)* sinp_geo(ip) &
                - r_geo(it,ip,3)*it*sint_geo(it)* cosp_geo(ip) &
                - r_geo(it,ip,4)*it*sint_geo(it)* sinp_geo(ip) 

           dzdt0 = dzdt0 + z_geo(it,ip,1)*it*cost_geo(it)* cosp_geo(ip) &
                + z_geo(it,ip,2)*it*cost_geo(it)* sinp_geo(ip) &
                - z_geo(it,ip,3)*it*sint_geo(it)* cosp_geo(ip) &
                - z_geo(it,ip,4)*it*sint_geo(it)* sinp_geo(ip)
        enddo
     enddo

     drdp0  = 0.0
     dzdp0  = 0.0
     do it=0,nts_geo
        do ip=0,nps_geo
           drdp0 = drdp0 - r_geo(it,ip,1)*ip*sint_geo(it)* sinp_geo(ip) &
                + r_geo(it,ip,2)*ip*sint_geo(it)* cosp_geo(ip) &
                - r_geo(it,ip,3)*ip*cost_geo(it)* sinp_geo(ip) &
                + r_geo(it,ip,4)*ip*cost_geo(it)* cosp_geo(ip) 
           
           dzdp0 = dzdp0 - z_geo(it,ip,1)*ip*sint_geo(it)* sinp_geo(ip) &
                + z_geo(it,ip,2)*ip*sint_geo(it)* cosp_geo(ip) &
                - z_geo(it,ip,3)*ip*cost_geo(it)* sinp_geo(ip) &
                + z_geo(it,ip,4)*ip*cost_geo(it)* cosp_geo(ip) 
        enddo
     enddo

     drdtt0 = 0.0
     dzdtt0 = 0.0
     do it=0,nts_geo
        do ip=0,nps_geo
           drdtt0 = drdtt0 - r_geo(it,ip,1)*it**2*sint_geo(it)* cosp_geo(ip) &
                - r_geo(it,ip,2)*it**2*sint_geo(it)* sinp_geo(ip) &
                - r_geo(it,ip,3)*it**2*cost_geo(it)* cosp_geo(ip) &
                - r_geo(it,ip,4)*it**2*cost_geo(it)* sinp_geo(ip) 
           
           dzdtt0 = dzdtt0 - z_geo(it,ip,1)*it**2*sint_geo(it)* cosp_geo(ip) &
                - z_geo(it,ip,2)*it**2*sint_geo(it)* sinp_geo(ip) &
                - z_geo(it,ip,3)*it**2*cost_geo(it)* cosp_geo(ip) &
                - z_geo(it,ip,4)*it**2*cost_geo(it)* sinp_geo(ip) 
        enddo
     enddo

     deallocate(cost_geo)
     deallocate(sint_geo)
     deallocate(cosp_geo)
     deallocate(sinp_geo)

  case default
     print *, 'ERROR: (LE3) Invalid equilibrium_model'
     stop

  end select

  ! J [dtb/dr terms vanish]
  jac0 = r0*(drdr0*dzdt0-drdt0*dzdr0)
  
  ! radius of curvature r_c (this is independent of theta coordinates).
  rc0 = (drdt0**2+dzdt0**2)**1.5/(drdt0*dzdtt0-dzdt0*drdtt0)  

end subroutine le3_rz
