module cgyro_gk

  implicit none

  public :: GK_alloc, GK_init, GK_do

  logical, private :: initialized = .false.

  complex, dimension(:,:,:), allocatable, private :: rhs
  complex, dimension(:,:), allocatable, private :: h0_x

  ! theta derivative variables
  integer, dimension(:), allocatable, private :: thcyc
  integer, dimension(:,:), allocatable, private :: rcyc
  real, dimension(-2:2), private :: uderiv
  real, dimension(-2:2), private :: cderiv
  complex, dimension(:,:), allocatable, private :: dtheta
  complex, dimension(:,:), allocatable, private :: dtheta_up

  real, dimension(:), allocatable :: vec_in, vec_outr, vec_outi

contains

  subroutine GK_alloc(flag)

    use cgyro_globals
    use cgyro_equilibrium, only : d_theta

    implicit none

    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it, jt, id, ir, jr
    complex :: thfac

    if (flag == 1) then

       if (initialized) return

       allocate(rhs(4,nc,nv_loc))
       allocate(h0_x(nc,nv_loc))

       ! cyclic index (for theta-periodicity)
       allocate(thcyc(1-n_theta:2*n_theta))
       do it=1,n_theta
          thcyc(it-n_theta) = it
          thcyc(it) = it
          thcyc(it+n_theta) = it
       enddo
       ! coefficients for 4th order centered derivative
       cderiv(-2) =  1.0 / (12.0 * d_theta)
       cderiv(-1) = -8.0 / (12.0 * d_theta)
       cderiv(0)  =  0.0 / (12.0 * d_theta)
       cderiv(1)  =  8.0 / (12.0 * d_theta)
       cderiv(2)  = -1.0 / (12.0 * d_theta)
       ! coefficients for 4th order filter for 3rd order upwinded derivative
       uderiv(-2) =  1.0 / (12.0 * d_theta)
       uderiv(-1) = -4.0 / (12.0 * d_theta)
       uderiv(0)  =  6.0 / (12.0 * d_theta)
       uderiv(1)  = -4.0 / (12.0 * d_theta)
       uderiv(2)  =  1.0 / (12.0 * d_theta)

       allocate(vec_in(n_xi))
       allocate(vec_outr(n_xi))
       allocate(vec_outi(n_xi))

       ! Indices for parallel streaming with upwinding
       allocate(rcyc(n_radial,1-n_theta:2*n_theta))
       allocate(dtheta(n_theta,-2:2))
       allocate(dtheta_up(n_theta,-2:2))
       do ir=1,n_radial
          do it=1,n_theta
             do id=-2,2
                jt = thcyc(it+id)
                if (it+id < 1) then
                   thfac = cos(2*pi*k_theta*rmin) &
                        + i_c * sin(2*pi*k_theta*rmin)
                   if (ir-1 >= 1) then
                      jr = ir-1
                   else
                      jr = n_radial
                   endif
                else if (it+id > n_theta) then
                   thfac = cos(2*pi*k_theta*rmin) &
                        - i_c * sin(2*pi*k_theta*rmin)
                   if (ir+1 <= n_radial) then
                      jr = ir+1
                   else
                      jr = 1
                   endif
                else
                   thfac = 1.0
                   jr = ir
                endif
                dtheta(it,id)    = cderiv(id)*thfac
                dtheta_up(it,id) = uderiv(id)*thfac
                rcyc(ir,it+id)   = jr

             enddo
          enddo
       enddo

       initialized = .true.

    else

       if (.not. initialized) return

       deallocate(rhs)
       deallocate(h0_x)
       deallocate(thcyc)
       deallocate(rcyc)
       deallocate(dtheta)
       deallocate(dtheta_up)
       deallocate(vec_in)
       deallocate(vec_outr)
       deallocate(vec_outi)

       initialized = .false.

    endif

  end subroutine GK_alloc
  
  subroutine GK_init

    use cgyro_globals
    use cgyro_field

    implicit none

    real :: ang

    h_x(:,:) = (0.0,0.0)


    if (restart_mode == 1) then

       open(unit=io_run,file=trim(path)//runfile_restart,status='old')
       read(io_run,*) h_x
       close(io_run)

    else

       if (zf_test_flag == 1) then

          ! Zonal-flow initial condition

          iv_loc = 0
          do iv=nv1,nv2
             iv_loc = iv_loc+1
             if (is_v(iv) == 1) then
                h_x(:,iv_loc) = 1.0
             else
                h_x(:,iv_loc) = 0.0
             endif
          enddo

       else

          ! Initial condition exponential in ballooning angle.

          iv_loc = 0
          do iv=nv1,nv2
             iv_loc = iv_loc+1
             if (is_v(iv) == 1) then
                do ic=1,nc
                   ang = theta(it_c(ic))+2*pi*(ir_c(ic)-n_radial/2-1)
                   h_x(ic,iv_loc) = rho*exp(-(ang/2)**2) 
                enddo
             endif
          enddo

       endif
    endif

    call FIELDx_do

    field_old = field

  end subroutine GK_init

  subroutine GK_do

    use cgyro_globals
    use cgyro_field

    implicit none

    ! compute new collisionless little_h: h = H - ze/T (G phi - G vpar/c apar)
    ! assumes have h_x, cap_h_x, and fields
    ! RK4

    h0_x = h_x

    ! Stage 1
    call get_gkRHS(1)
    h_x = h0_x + 0.5 * delta_t * rhs(1,:,:)
    call FIELDx_do

    ! Stage 2
    call get_gkRHS(2)
    h_x = h0_x + 0.5 * delta_t * rhs(2,:,:)
    call FIELDx_do

    ! Stage 3
    call get_gkRHS(3)
    h_x = h0_x + delta_t * rhs(3,:,:)
    call FIELDx_do

    ! Stage 4
    call get_gkRHS(4)
    h_x = h0_x + delta_t/6.0 * &
         (rhs(1,:,:)+2.0*rhs(2,:,:)+2.0*rhs(3,:,:)+rhs(4,:,:))  
    call FIELDx_do

  end subroutine GK_do
  
  subroutine get_gkRHS(ij)

    use parallel_lib
    use timer_lib

    use cgyro_globals
    use cgyro_equilibrium

    implicit none

    integer, intent(in) :: ij
    integer :: is, ir, it, ie, ix, jx
    integer :: id, jt, jr, jc
    real    :: rval
    complex :: val

    call timer_lib_in('gkrhs')

    ! Get the RHS collisionless GK little_h: h = H-ze/T (G phi - G vpar/c apar)
    rhs(ij,:,:) = (0.0,0.0)

    iv_loc = 0
    do iv=nv1,nv2

       iv_loc = iv_loc+1

       is = is_v(iv)
       ix = ix_v(iv)
       ie = ie_v(iv)

       do ic=1,nc

          ir = ir_c(ic) 
          it = it_c(ic)

          ! parallel streaming with upwinding
          do id=-2,2
             jt = thcyc(it+id)
             jr = rcyc(ir,it+id)
             jc = ic_c(jr,jt)
             rval = omega_stream(it,is) * sqrt(energy(ie)) * xi(ix) 
             rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
                  - rval*dtheta(it,id)*cap_h_c(jc,iv_loc)  &
                  - abs(rval)*dtheta_up(it,id)*h_x(jc,iv_loc)
          enddo

          ! omega_rdrift
          val = omega_rdrift(it,is) &
               * energy(ie) * (1.0 + xi(ix)**2) &
               * (2.0*pi*i_c*indx_r(ir)*r_length_inv) 
          rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
               - val * cap_h_c(ic,iv_loc)

          ! radial upwind
          val = omega_rdrift(it,is) &
               * energy(ie) * (1.0 + xi(ix)**2)
          val = abs(val) * up_radial & 
               * (2.0*indx_r(ir)/(1.0*n_radial))**(up_radial_n-1.0) &
               * (2.0*pi*indx_r(ir)*r_length_inv)
          rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
               - val * h_x(ic,iv_loc)

          ! omega_dalpha
          val = omega_adrift(it,is) &
               * energy(ie) * (1.0 + xi(ix)**2) &
               * (i_c * k_theta)
          rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
               - val * cap_h_c(ic,iv_loc)

          ! omega_dalpha - pressure component
          val = omega_aprdrift(it,is) &
               * energy(ie) * xi(ix)**2 &
               * (i_c * k_theta)
          rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
               - val * cap_h_c(ic,iv_loc)

          ! omega_star
          val = i_c * k_theta &
               * rho *sqrt(temp(is)*mass(is))/(1.0*z(is)) * vth(is) &
               * (dlnndr(is) + dlntdr(is) * (energy(ie)-1.5))
          rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
               - val * z(is)/temp(is) &
               * gyrox_J0(is,ir,it,ie,ix) * field(ir,it,1)
          if (n_field > 1) then
             rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
                  - val * z(is)/temp(is) &
                  * gyrox_J0(is,ir,it,ie,ix) * field(ir,it,2) &
                  * (-xi(ix) * sqrt(2.0*energy(ie)) * vth(is))
          endif
          

       enddo
    enddo

    ! TRAPPING TERM

    if (collision_model == 0) then

       call parallel_lib_r(transpose(cap_h_c),cap_h_v)
       cap_h_v_prime(:,:) = (0.0,0.0)
       ic_loc = 0
       do ic=nc1,nc2
          ic_loc = ic_loc+1
          it = it_c(ic)
          ir = ir_c(ic)
          do iv=1,nv
             is = is_v(iv)
             ix = ix_v(iv)
             ie = ie_v(iv)
             
             do jx=1,n_xi
                cap_h_v_prime(ic_loc,iv) = cap_h_v_prime(ic_loc,iv) &
                     + xi_deriv_mat(ix,jx) * cap_h_v(ic_loc,iv_v(ie,jx,is))
             enddo
          enddo
       enddo
       
       ! Now have cap_h_v(ic_loc,iv)   
       
       call parallel_lib_f(cap_h_v_prime,cap_h_ct)
       cap_h_c = transpose(cap_h_ct)
       
       iv_loc = 0
       do iv=nv1,nv2
       
          iv_loc = iv_loc+1
          
          is = is_v(iv)
          ix = ix_v(iv)
          ie = ie_v(iv)
          
          do ic=1,nc
             
             ir = ir_c(ic) 
             it = it_c(ic)
             val = omega_trap(it,is) &
                  * sqrt(energy(ie)) &
                  * (1.0 - xi(ix)**2) 
             rhs(ij,ic,iv_loc) = rhs(ij,ic,iv_loc) &
                  - val * cap_h_c(ic,iv_loc)
          enddo
       enddo
    endif

    call timer_lib_out('gkrhs')

  end subroutine get_gkRHS
   
 end module cgyro_gk
