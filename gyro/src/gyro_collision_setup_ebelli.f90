!---------------------------------------------------
! gyro_collision_setup_ebelli.f90 
!
! PURPOSE:
!  Make time-advance matrices, combining RBF
!  expansion in space with Crank-Nicholson
!  scheme in time.  
! 
! NOTES:
!  The parameter ord_rbf=3 controls the order
!  of the RBF expansion.  
!
!  The irregular gridset in the (xi,theta) plane
!  is expanded in an RBF spline series.
!
!  Near the boundaries at xi=-1 and xi=1, the SNaK
!  method is used to improve accuracy.
!
!  The pitch-angle scattering operator is
!
!                d       2   d
!  L = nu_coll_d --- (1-xi ) ---
!                dxi         dxi
!
!
!  L[v_par] = -2 nu_coll_d v_par
!
!---------------------------------------------------

subroutine gyro_collision_setup_ebelli
  
  use gyro_globals
  use gyro_collision_private
  use gyro_pointers
  use math_constants
  
  !-----------------------------------------------------------
  implicit none
  !
  integer :: js, ks, ia, ib, ja, jb
  integer :: p, pp
  !
  real :: h_coll, x_coll, xa, xb, sum_nu
  real :: r1, r2, nu_fac
  real, dimension(n_rbf) :: xp
  real, dimension(n_rbf) :: yp
  real, dimension(n_rbf) :: xc
  real, dimension(n_rbf,n_rbf) :: a1, a_inv, a_lor, a_res
  real, dimension(n_rbf,n_rbf) :: d1, d2
  ! LAPACK
  integer :: info1, info2
  integer, dimension(:), allocatable :: i_piv1, i_piv2
  real, dimension(:), allocatable :: work1, work2
  !
  real, external :: DERF
  real, external :: rbf_ebelli

  !-----------------------------------------------------------

  ! nu_total is not used
  nu_total(:,:,:) = 0.0

  if(.not. (collision_method==3 .or. collision_method==4)) then
     call catch_error('ERROR: Collision_method must be 3 (Connor) or 4 (HS0)')
  endif
  
  if(i_proc==0) then
     print *, '******* collision_method = ', collision_method
  endif

  ! Compute the collision frequencies
  do i=1, n_x
     do is=1,n_kinetic
        do js=1,n_kinetic
           rs_coll_const(i,is,js) = 0.0
           
           do ie=1,n_energy
              
              xa = sqrt(energy(ie,1))
              xb = xa * mu(is)/mu(js) * sqrt(tem_s(is,i)/tem_s(js,i))
              nu_fac = nu_s(is,i) * den_s(js,i) / den_s(is,i) * z(js)**2

              if(collision_method == 3) then
                 ! Connor model: nu_d form dependent on (is,js)
                 if(is == js .or. &  
                      (abs(mu(is) - mu(js)) < epsilon(0.))) then
                    ! case 1: like-species/same mass collisions
                    h_coll = exp(-xb*xb)/(xb*sqrt(pi)) &
                         + (1.0-1.0/(2.0*xb*xb)) * DERF(xb)
                    x_coll = 1.0/xa**3
                    
                 else if(mu(is) > mu(js)) then
                    ! case 2: ele-ion and ion-imp(heavy) collisions
                    h_coll = 1.0
                    x_coll = 1.0/xa**3
                    
                 else
                    ! case 3: ion-ele and imp(heavy)-ion collisions
                    h_coll = 1.0
                    x_coll = 1.0
                    nu_fac = nu_fac * 4.0/(3.0*sqrt(pi)) &
                         * mu(is) / mu(js) * (tem_s(is,i) / tem_s(js,i))**1.5

                 endif
                 
                 nu_coll_d(i,ie,is,js) = nu_fac * h_coll * x_coll
                 
              else 
                 ! Connor-like model with symmetric nu_d
                 
                 ! deflection frequency
                 h_coll = exp(-xb*xb)/(xb*sqrt(pi)) &
                      + (1.0-1.0/(2.0*xb*xb)) * DERF(xb)
                 x_coll = 1.0/xa**3
                 nu_coll_d(i,ie,is,js) = nu_fac * h_coll * x_coll
                 
              endif
              
              ! Rs restoring term constant
              rs_coll_const(i,is,js) = rs_coll_const(i,is,js) &
                   + w_energy(ie,1) *  energy(ie,1) * nu_coll_d(i,ie,is,js) 
           enddo ! ie

           if(abs(rs_coll_const(i,is,js)) > 0.0) then
              rs_coll_const(i,is,js) = 1.5 * sqrt(tem_s(js,i)/tem_s(is,i)) &
                   * mu(is)/mu(js) * den_s(js,i)/den_s(is,i) &
                   / rs_coll_const(i,is,js)
           endif
           
        enddo ! js
        
        do js=1,n_kinetic
           do ks=1,n_kinetic
              rs_nunu_const(i,is,js,ks) = 0.0
              do ie=1,n_energy
                 rs_nunu_const(i,is,js,ks) = rs_nunu_const(i,is,js,ks) &
                      + w_energy(ie,1) *  energy(ie,1) &
                      * nu_coll_d(i,ie,is,js) * nu_coll_d(i,ie,is,ks)
              enddo
           enddo
        enddo
        
     enddo ! is
  enddo    ! i

  !-----------------------------------------------------------

  if (electron_method == 3) then
     call catch_error('ERROR: Collisions are not valid for electron_method=3')
  endif
  
  !-----------------------------------------------------------------
  ! Remove n=0 collisions?
 
  if (kill_coll_flag == 1 .and. n_1(in_1) == 0) then
     nu_coll_d(:,:,:,:) = 0.0
     rs_coll_const(:,:,:) = 0.0
     rs_nunu_const(:,:,:,:) = 0.0
  endif
 
  !-----------------------------------------------------------------

  !------------------------------------------------
  ! Define xi = v_parallel/v, including sign of 
  ! velocity.
  !
  do k=1,n_lambda
     ck = class(k)
     do m=1,n_stack
        if (m <= n_theta(ck)) then
           xi(:,k,m) = -sqrt(abs(1.0-lambda(:,k)*b0_t(:,k,m))) 
        else
           xi(:,k,m) = sqrt(abs(1.0-lambda(:,k)*b0_t(:,k,m))) 
        endif
     enddo
  enddo
  !
  !--------------------------------------------------

  p=0
  do is=1,n_kinetic
     do js=1,n_kinetic
        p = p + 1
        indx_coll(is,js) = p
     end do
  end do

  allocate(i_piv1(n_rbf))
  allocate(work1(n_rbf))
  allocate(i_piv2(n_kinetic**2))
  allocate(work2(n_kinetic**2))

  p_ine_loc = 0
  
  do p_ine = 1+i_proc_1,n_ine_1,n_proc_1

     p_ine_loc = p_ine_loc+1
     
     i  = ine_i(p_ine)
     ie = ine_e(p_ine)
     
     ! Pack single-index RBF array
     ! -1 < xp < 1
     ! -1 <= yp < 1
     p = 0
     do k=1,n_lambda
        do m=1,n_stack

           p = p+1
           xp(p) = xi(i,k,m)
           xc(p) = xi(i,k,m)
           yp(p) = theta_t(i,k,m)/pi

           ! Shift centers near boundary according to SNaK 
           ! method for improved boundary accuracy and 
           ! significantly better time-advance stability. 
           if (k <= 2) then
              if (xp(p) < 0.0) then
                 xc(p) = -1.0-0.25*k
              else
                 xc(p) = 1.0+0.25*k
              endif
           endif

        enddo
     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Build RBF expansion matrix
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do p=1,n_rbf
        do pp=1,n_rbf
           a_inv(p,pp) = rbf_ebelli(pi*(yp(p)-yp(pp)),xp(p)-xc(pp),ord_rbf)
        enddo
     enddo

     ! LAPACK factorization and inverse
     call DGETRF(n_rbf,n_rbf,a_inv,n_rbf,i_piv1,info1)
     call DGETRI(n_rbf,a_inv,n_rbf,i_piv1,work1,n_rbf,info1)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Lorentz derivative operator 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do p=1,n_rbf
        do pp=1,n_rbf
           call drbf_ebelli(pi*(yp(p)-yp(pp)),xp(p)-xc(pp),ord_rbf,r1,r2)
           a1(p,pp) = (1-xp(p)**2)*r1 - 2*xp(p)*r2
        enddo
     enddo
     call DGEMM('N','N',n_rbf,n_rbf,n_rbf,1.0,a1,n_rbf,a_inv,n_rbf,0.0,&
          a_lor,n_rbf)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Momentum restoring xi integral operator
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do p=1,n_rbf
        do pp=1,n_rbf
           call irbf_ebelli(pi*(yp(p)-yp(pp)),xc(pp),ord_rbf,r1)
           a1(p,pp) = r1
        enddo
     enddo
     call DGEMM('N','N',n_rbf,n_rbf,n_rbf,1.0,a1,n_rbf,a_inv,n_rbf,0.0,&
          a_res,n_rbf)

     do is=1,n_kinetic
        d1(:,:) = 0.0
        d2(:,:) = 0.0

        sum_nu=0.0
        do js=1,n_kinetic
           sum_nu = sum_nu + nu_coll_d(i,ie,is,js)
        enddo


        ! Identity matrix
        do p=1,n_rbf
           d1(p,p) = 1.0
           d2(p,p) = 1.0
        enddo

        ! "Implicit" part (inverse)
        d1(:,:) = d1(:,:) - 0.5 * dt * sum_nu * 0.5 * a_lor(:,:)
        call DGETRF(n_rbf,n_rbf,d1,n_rbf,i_piv1,info1)
        call DGETRI(n_rbf,d1,n_rbf,i_piv1,work1,n_rbf,info1)

        ! EAB test -- only R
        !a_lor(:,:) = 0.0
        ! EAB test -- only L
        !a_res(:,:) = 0.0

        ! "Explicit" part 
        d2(:,:) = d2(:,:) + 0.5 * dt * sum_nu * 0.5 * a_lor(:,:)

        call DGEMM('N','N',n_rbf,n_rbf,n_rbf,1.0,d2,n_rbf,d1,n_rbf,0.0,&
             a1,n_rbf)
        d_rbf_lorentz(:,:,is,p_ine_loc) = a1(:,:)

        ! collisional R_s -- fast matrix multiply
        do pp=1,n_rbf
           a1(:,pp) = xp(:) * a_res(:,pp)
        enddo
        call DGEMM('N','N',n_rbf,n_rbf,n_rbf,1.0,&
             d_rbf_lorentz(:,:,is,p_ine_loc),n_rbf,a1,n_rbf,0.0,d1,n_rbf)
        do js=1,n_kinetic
           d_rbf_rs(:,:,is,js,p_ine_loc) = d1(:,:) &
                * sqrt(energy(ie,1)) * nu_coll_d(i,ie,is,js) &
                * rs_coll_const(i,is,js) * 0.5 * dt 
        enddo
                
        ! Matrices for Rs int d^3 sqrt(ene)*nu*xi equation
        ! -- int d_xi xi (...)

        call DGEMM('N','N',n_rbf,n_rbf,n_rbf,1.0,a_res(:,:),n_rbf,&
             d_rbf_lorentz(:,:,is,p_ine_loc),n_rbf,0.0,&
             d_rbf_lorentz_int(:,:,is,p_ine_loc),n_rbf)

        do js=1,n_kinetic
           call DGEMM('N','N',n_rbf,n_rbf,n_rbf,1.0,a_res(:,:),n_rbf,&
                d_rbf_rs(:,:,is,js,p_ine_loc),n_rbf,0.0,&
                d_rbf_rs_int(:,:,is,js,p_ine_loc),n_rbf)
        enddo
        
     enddo

     ! Matrix for Rs int d^3 sqrt(ene)*nu*xi equation
     d_rbf_velint(:,:,p_ine_loc) = 0.0
     do ia=1,n_kinetic
        do ib=1,n_kinetic
           p = indx_coll(ia,ib)
           do ja=1,n_kinetic
              do jb=1,n_kinetic
                 pp = indx_coll(ja,jb)

                 if(jb == ia) then
                    d_rbf_velint(p,pp,p_ine_loc) =  &
                         d_rbf_velint(p,pp,p_ine_loc) &
                         -0.5 * dt  * (2.0/3.0) &
                         * rs_coll_const(i,ia,ja) * rs_nunu_const(i,ia,ja,ib)
                 endif

                 if(p == pp) then
                    d_rbf_velint(p,pp,p_ine_loc) = &
                         d_rbf_velint(p,pp,p_ine_loc) + 1.0
                 endif

              enddo
           enddo
        enddo
     enddo

     ! transpose matrices for optimal index order
     do ks=1,n_kinetic
        d_rbf_lorentz(:,:,ks,p_ine_loc) = &
             transpose(d_rbf_lorentz(:,:,ks,p_ine_loc))
        d_rbf_lorentz_int(:,:,ks,p_ine_loc) = &
             transpose(d_rbf_lorentz_int(:,:,ks,p_ine_loc))
     enddo
     do js=1,n_kinetic
        do ks=1,n_kinetic
           d_rbf_rs(:,:,ks,js,p_ine_loc) = &
                transpose(d_rbf_rs(:,:,ks,js,p_ine_loc))
           d_rbf_rs_int(:,:,ks,js,p_ine_loc) = &
                transpose(d_rbf_rs_int(:,:,ks,js,p_ine_loc))
        enddo
     enddo

     d_rbf_velint(:,:,p_ine_loc) = transpose(d_rbf_velint(:,:,p_ine_loc))


     ! LAPACK factorization and inverse
     call DGETRF(n_kinetic**2,n_kinetic**2,d_rbf_velint(:,:,p_ine_loc),&
          n_kinetic**2,i_piv2,info2)
     call DGETRI(n_kinetic**2,d_rbf_velint(:,:,p_ine_loc),n_kinetic**2,i_piv2,&
          work2,n_kinetic**2,info2)

  enddo ! p_ine

  deallocate(i_piv1)
  deallocate(work1)
  deallocate(i_piv2)
  deallocate(work2)

  if (i_proc == 0) call gyro_collision_grid_write(trim(path)//'gyro_collision_grid.out',1)

  if (i_proc == 0 .and. debug_flag == 1) &
       print *,'[make_collision_stencil called]'

end subroutine gyro_collision_setup_ebelli


real function rbf_ebelli(dt,dy,s)
  implicit none
  real, intent(in) :: dt,dy
  integer, intent(in) :: s
  real :: rl

  rl  = sqrt(2.0-2.0*cos(dt)+dy**2)
  rbf_ebelli = rl**s

end function rbf_ebelli


subroutine drbf_ebelli(dt,dy,s,r1,r2)
  implicit none
  real, intent(in) :: dt,dy
  integer, intent(in) :: s
  real, intent(inout) :: r1,r2
  real :: rL

  rL = sqrt(2.0-2.0*cos(dt)+dy**2)
  r2 = s*rL**(s-2)*dy

  if (rL > 0.0) then
     r1 = (s*(s-2)*dy**2+s*rL**2)*rL**(s-4)
  else
     r1 = 0.0
  endif

end subroutine drbf_ebelli

subroutine irbf_ebelli(dt,y,s,r1)
  use gyro_globals, only : i_proc
  implicit none
  real, intent(in) :: dt, y
  integer, intent(in) :: s
  real, intent(inout) :: r1
  real :: rL, tL, logfac, sgn
  integer :: iy
  real, dimension(2) :: dy
  
  if(.not.(s==3)) then
     call catch_error('ERROR: Collisions are only valid for ord_rbf=3')
  endif
  
  ! integral of d_xi * xi * rbf
  dy(1) =  1.0-y
  dy(2) =  -1.0-y
  tL = 2.0-2.0*cos(dt)

  sgn=1.0
  r1 = 0.0
  do iy=1,2
     rL = sqrt(tL+dy(iy)**2)

     if(dy(iy) + rL <= 0.0) then
        logfac = 0.0
     else
        logfac = tL**2 * log(dy(iy) + rL)
     endif

        r1 = r1 + sgn * ( 1.0/5.0 * rL**5 &
             + y*(1.0/4.0 * dy(iy)**3 * rL &
             + 5.0/8.0 * tL * dy(iy) * rL &
             + 3.0/8.0 * logfac))

        sgn = -1.0
     enddo
  
end subroutine irbf_ebelli
