!-----------------------------------------------------------------
! gyro_fieldeigen.f90
!
! PURPOSE:
!  Solves the Maxwell dispersion matrix eigenvalue problem (in 
!  other word, the field equations) for the real frequency
!  and growth rate and computes the corresponding eigenfunctions.
!  We make heavy use of BLAS routines to ensure maximum performance
!  of most all matrix operations.
!
! INPUT PARAMETERS:
!  FIELDEIGEN_ROOT_METHOD
!  FIELDEIGEN_WR
!  FIELDEIGEN_WI
!  FIELDEIGEN_TOL
!    
! OUTPUT:
!  - fieldeigen.out
!  - eigenfunctions are written in balloon_*.out as usual.
!------------------------------------------------------------------

subroutine gyro_fieldeigen

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants
  use gyro_collision_private, only : nu_total
  use gyro_fieldeigen_private

  !---------------------------------------------------------------
  implicit none
  !
  ! Variables for calculation of parallel motion
  !
  real    :: v0
  integer :: mff,mf
  integer :: mc
  integer :: mb,mbb
  integer :: p  
  complex :: pff,pf
  complex :: pc
  complex :: pb,pbb
  !
  ! Additional variables
  !
  complex :: f1,f2,f3
  complex :: z1,z2,z3
  complex :: dz
  integer :: iw
  complex :: a(3,3)
  complex :: fv(3)
  integer :: ipiv(3)
  integer :: ierr
  !---------------------------------------------------

  !-------------------------------------------------------------------------
  ! Validity checks
  !
  if (boundary_method == 2) then
     call catch_error('ERROR: BOUNDARY_METHOD=2 not supported in fieldeigen.')
  endif

  if (n_1(in_1) == 0) then
     call catch_error('ERROR: n=0 not supported in fieldeigen.')
  endif

  if (electron_method == 2) then
     call catch_error('ERROR: ELECTRON_METHOD=2 not supported in fieldeigen.')
  endif
  !-------------------------------------------------------------------------

  if (i_proc == 0 .and. output_flag == 1) then
     open(unit=1,file='fieldeigen.out',status='replace')
     close(1)
  endif
  if (i_proc == 0 .and. silent_flag == 0) then
     print *,' Re(omega)      Im(omega)      |det|          error'
  endif

  !-------------------------------------------------------------
  ! Collisionless allocations
  !-------------------------------------------------------------

  n_eigen = n_x* n_blend*n_field
  n_im    = n_x*n_stack

  allocate(vdotgrad(n_im,n_im,n_lambda,n_energy/n_proc,n_kinetic))
  allocate(v_omegastar(n_x,n_lambda*n_stack))
  allocate(fg(n_blend,-m_gyro:m_gyro-i_gyro,n_x,n_lambda*n_stack,n_field))

  allocate(a_eigen(n_eigen,n_eigen))
  allocate(a_eigen_loc(n_eigen,n_eigen))
  allocate(b_eigen(n_eigen))
  allocate(i_piv_eigen(n_eigen))

  allocate(i_piv_im(n_im))
  allocate(propinv(n_im,n_im))
  allocate(prod(n_eigen,n_im))
  allocate(work_im(n_im))
  allocate(gk_left(n_eigen,n_im))
  allocate(gk_right(n_im,n_eigen))

  allocate(diag_scale(n_eigen))

  diag_scale(:) = -1.0

  !-----------------------------------------------------------------------
  ! Extra setup required for collisions
  !
  if (collision_flag == 1) then

     n_imk = n_x*n_stack*n_lambda

     allocate(cmk(n_lambda*n_stack,n_x))
     allocate(nu_op(n_lambda*n_stack,n_lambda*n_stack,n_x,n_energy/n_proc))

     allocate(i_piv_imk(n_imk))
     allocate(propinvk(n_imk,n_imk))
     allocate(prodk(n_eigen,n_imk))
     allocate(work_imk(n_imk))
     allocate(gk_leftk(n_eigen,n_imk))
     allocate(gk_rightk(n_imk,n_eigen))

     ! Store collision operator in canonical storage format 
     ! for fieldeigen.

     ! (1) define phase function

     do i=1,n_x
        do m=1,n_stack
           do k=1,n_lambda 
              mk = m+(k-1)*n_stack
              cmk(mk,i) = exp(i_c*angp(i)*theta_t(i,k,m))  
           enddo
        enddo
     enddo

     ! (2) apply phase function and collision frequency to RBF derivative:

     do ie2=1,n_energy/n_proc
        ie = ie2+i_proc_1*n_energy/n_proc
        do i=1,n_x
           do mk=1,n_lambda*n_stack 
              do mkp=1,n_lambda*n_stack

                 nu_op(mkp,mk,i,ie2) = nu_total(i,ie,indx_e)*&
                      cmk(mkp,i)*d1_rbf(mkp,mk)/cmk(mk,i) 

              enddo ! mkp
           enddo ! mk
        enddo ! i
     enddo ! ie2

     deallocate(cmk)

  endif
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! H-independent field equation pieces:
  !
  !  [-lambda^2 grad^2] phi
  !  [-2*rho^2/beta grad^2 ] Apar
  !  [1/beta_e] Bpar 
  !
  ! Note: do not need electron current pieces since em effects are
  ! not valid for electron_method=1 and eigensolve not valid for
  ! electron_method=2.
  !
  allocate(ap_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
  call gyro_blend_poisson(1)

  if (n_field > 1) then
     allocate(aa_mm(n_x,-mg_dx:mg_dx-ig_dx,n_blend,n_blend))
     call gyro_blend_ampere
  endif

  if (n_field > 2) then
     allocate(ab_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
     allocate(abp_mm(n_x,-m_gyro:m_gyro-i_gyro,n_blend,n_blend))
     call gyro_blend_ampereperp
  endif
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! vdotgrad = (vpar b dot grad + v_d dot grad)
  !
  vdotgrad(:,:,:,:,:) = (0.0,0.0)

  p_nek_loc = 0

  do ie2=1,n_energy/n_proc
     ie =  ie2+i_proc_1*n_energy/n_proc
     do k=1,n_lambda

        p_nek_loc = p_nek_loc+1

        ck = class(k)

        do is=1,n_kinetic
           do i=1,n_x
              do m=1,n_stack

                 m0 = m_phys(ck,m)
                 im = i + (m-1)*n_x

                 ! d/dr part of the drift + upwind
                 do i_diff=-m_dx,m_dx-i_dx

                    ip = i_cyc(i+i_diff)

                    impr = ip + (m-1)*n_x

                    vdotgrad(im,impr,k,ie2,is) = &
                         vdotgrad(im,impr,k,ie2,is) &
                         - w_d1(i_diff)*omega_dr(m,i,p_nek_loc,is) &
                         - s_d1(i_diff)*abs(omega_dr(m,i,p_nek_loc,is))

                 enddo


                 vdotgrad(im,im,k,ie2,is) = &
                      vdotgrad(im,im,k,ie2,is) &
                      - i_c * omega_d1(m,i,p_nek_loc,is)

                 ! Sign flip for LHS
                 v0 = -mu(is)*sqrt(tem_s(is,i))*v_theta(i,ie,k,is)/d_tau(ck)

                 p = p_phys(ck,m)

                 mff = m_cyc(ck,m+2,p)
                 mf  = m_cyc(ck,m+1,p)
                 mc  = m_cyc(ck,m,p)
                 mb  = m_cyc(ck,m-1,p)
                 mbb = m_cyc(ck,m-2,p)

                 pff = p_cyc(ck,i,m+2,p)
                 pf  = p_cyc(ck,i,m+1,p)
                 pc  = p_cyc(ck,i,m,p)
                 pb  = p_cyc(ck,i,m-1,p)
                 pbb = p_cyc(ck,i,m-2,p)

                 im  = i + (m-1)*n_x

                 impr = i + (mff-1)*n_x
                 vdotgrad(im,impr,k,ie2,is) = &
                      vdotgrad(im,impr,k,ie2,is) &
                      + (-1.0/12.0)*pff*sigma_tau(ck,p)*v0

                 impr = i + (mf-1)*n_x
                 vdotgrad(im,impr,k,ie2,is) = &
                      vdotgrad(im,impr,k,ie2,is) &
                      + (8.0/12.0)*pf*sigma_tau(ck,p)*v0

                 impr = i + (mb-1)*n_x   
                 vdotgrad(im,impr,k,ie2,is) = &
                      vdotgrad(im,impr,k,ie2,is) &
                      + (-8.0/12.0)*pb*sigma_tau(ck,p)*v0

                 impr = i + (mbb-1)*n_x
                 vdotgrad(im,impr,k,ie2,is) = &
                      vdotgrad(im,impr,k,ie2,is) &
                      + (1.0/12.0)*pbb*sigma_tau(ck,p)*v0

              enddo ! m
           enddo ! i
        enddo ! is
     enddo ! k
  enddo ! ie2
  !-----------------------------------------------------------------------

  select case (fieldeigen_root_method)

  case (1) 

     !-------------------------------------------------------------
     ! Secant method for calculation of root:
     !-------------------------------------------------------------

     ! Initial guess for omega

     z1 = fieldeigen_wr + i_c*fieldeigen_wi
     if (abs(fieldeigen_wr) < 1e-10 .and. abs(fieldeigen_wi) < 1e-10) then
        z2 = (fieldeigen_wr+0.001) +  i_c * (fieldeigen_wi + 0.001)
     else
        z2 = 1.001*z1
     endif


     do iw=-1,iwmax

        if (iw == -1) then
           dz = (1.0,0.0)
           omega_eigen = z1
        else if (iw == 0) then
           dz = (1.0,0.0)
           omega_eigen = z2
        else
           dz = (z1-z2)/(f1-f2)*f1
           if (abs(dz)/abs(omega_eigen) < fieldeigen_tol) then
              exit
           endif
           z2 = z1
           z1 = z1-dz
           omega_eigen = z1
        endif

        !  INPUT: omega_eigen
        ! OUTPUT: det
        call gyro_fieldeigen_kernel

        error_eigen = abs(dz)/abs(omega_eigen) 

        if (i_proc == 0 .and. silent_flag == 0) then

           print '(t2,5(1pe14.7,1x))',&
                real(omega_eigen),&
                aimag(omega_eigen), &
                abs(det),&
                error_eigen

        endif
        if (i_proc == 0 .and. output_flag == 1) then

           open(unit=1,file='fieldeigen.out',status='old',position='append')

           write(1,'(t2,5(1pe14.7,1x))') &
                real(omega_eigen), &
                aimag(omega_eigen), &
                abs(det), &
                error_eigen

           close(1)

        endif

        if (iw == -1) then
           f1 = det
        else if (iw == 0) then
           f2 = det
        else
           f2 = f1
           f1 = det
        endif

     enddo ! iw

  case (2)

     !-------------------------------------------------------------
     ! Rational function method for calculation of root:
     !-------------------------------------------------------------

     ! Initial guess for omega

     z1 = fieldeigen_wr + i_c*fieldeigen_wi
     z2 = z1+(0.01,0.01)
     z3 = z1-(0.01,0.01)

     do iw=-2,iwmax

        if (iw == -2) then
           dz = (1.0,0.0)
           omega_eigen = z1
        else if (iw == -1) then
           dz = (1.0,0.0)
           omega_eigen = z2
        else if (iw == 0) then
           dz = (1.0,0.0)
           omega_eigen = z3
        else

           fv(1) = f1
           fv(2) = f2
           fv(3) = f3

           a(1,1) = 1.0
           a(1,2) = z1
           a(1,3) = -z1*f1

           a(2,1) = 1.0
           a(2,2) = z2
           a(2,3) = -z2*f2

           a(3,1) = 1.0
           a(3,2) = z3
           a(3,3) = -z3*f3

           ! LAPACK matrix factorization into L/U components
           call ZGETRF(3,3,a,3,ipiv,ierr)

           ! LAPACK matrix solve
           call ZGETRS('N',3,1,a,3,ipiv,fv,3,ierr)

           z3 = z2
           z2 = z1
           z1 = -fv(1)/fv(2)

           if (abs(z2-z1)/abs(z1) < fieldeigen_tol) then
              exit
           endif

           omega_eigen = z1

        endif

        !  INPUT: omega_eigen
        ! OUTPUT: det
        call gyro_fieldeigen_kernel

        error_eigen = abs(z2-z1)/abs(omega_eigen) 

        if (i_proc == 0 .and. silent_flag == 0) then

           print '(t2,5(1pe14.7,1x))',&
                real(omega_eigen),&
                aimag(omega_eigen), &
                abs(det), &
                error_eigen

        endif
        if (i_proc == 0 .and. output_flag == 1) then

           open(unit=1,file='fieldeigen.out',status='old',position='append')

           write(1,'(t2,5(1pe14.7,1x))') &
                real(omega_eigen), &
                aimag(omega_eigen), &
                abs(det), &
                error_eigen

           close(1)

        endif

        if (iw == -2) then
           f1 = det
        else if (iw == -1) then
           f2 = det
        else if (iw == 0) then
           f3 = det
        else
           f3 = f2
           f2 = f1
           f1 = det
        endif

     enddo ! iw

  end select

  !  if (i_proc == 0 .and. output_flag == 1) then
  !     open(unit=1,file='fieldeigen_INPUT.out',status='replace')
  !     write(1,'(a,1pe12.4)') 'FIELDEIGEN_WR=',real(omega_eigen)
  !     write(1,'(a,1pe12.4)') 'FIELDEIGEN_WI=',aimag(omega_eigen)
  !     close(1)
  !  endif

  !----------------------------------------------------------------------------
  ! Compute the eigenfunction
  !
  if (iw < iwmax) then

     ! Re-create a_eigen with the last (converged) eigenvalue

     a_eigen(:,:) = a_eigen_loc(:,:)

     ! Impose an inhomogeneous constraint: c_phi(i_end,j_end) = 1.0
     i  = n_x
     j  = n_blend
     ix = 1
     ij = i + (j-1)*n_x
     do ip=1,n_x
        do jp=1,n_blend
           do ixp=1,n_field

              ijp = ip + (jp-1)*n_x + (ixp-1)*n_x*n_blend

              if (ip == n_x .and. jp == n_blend .and. ixp == ix) then
                 a_eigen(ij,ijp) = (1.0,0.0)
              else
                 a_eigen(ij,ijp) = (0.0,0.0)
              endif
           enddo
        enddo
     enddo
     b_eigen(:) = (0.0,0.0)
     b_eigen(n_x*n_blend) = (1.0,0.0)

     ! Factorize and solve
     call ZGETRF(n_eigen,n_eigen,a_eigen,n_eigen,i_piv_eigen,info_eigen)
     call ZGETRS('N',n_eigen,1,a_eigen,n_eigen,i_piv_eigen, &
          b_eigen,n_eigen,info_eigen)

     do i=1,n_x
        do j=1,n_blend
           ij  = i + (j-1)*n_x 
           field_blend(j,i,1) = b_eigen(ij)
        enddo
     enddo
     if (n_field > 1) then
        do i=1,n_x
           do j=1,n_blend
              ij  = i + (j-1)*n_x + n_x*n_blend
              field_blend(j,i,2) = b_eigen(ij)
           enddo
        enddo
     endif
     if (n_field > 2) then
        do i=1,n_x
           do j=1,n_blend
              ij  = i + (j-1)*n_x + 2*n_x*n_blend
              field_blend(j,i,3) = b_eigen(ij)
           enddo
        enddo
     endif

  else

     error_eigen = 1.0

  endif

  !-------------------------------------------------------------------
  ! Call plotting and other IO functions:
  !
  call gyro_fieldeigen_df
  call gyro_field_fluxave
  call gyro_field_time_derivative
  call get_field_plot
  call gyro_moments_plot
  if (io_method == 1) then 
     call gyro_write_timedata
  else
     call gyro_write_timedata_hdf5
  endif
  step = 1
  call gyro_write_precision(10,abs(omega_eigen))
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Clean-up

  deallocate(ap_mm)

  if (n_field > 1) then
     deallocate(aa_mm)
  endif

  if(n_field == 3) then
     deallocate(ab_mm)
     deallocate(abp_mm)
  endif

  deallocate(vdotgrad)
  deallocate(v_omegastar)
  deallocate(fg)

  deallocate(a_eigen)
  deallocate(a_eigen_loc)
  deallocate(b_eigen)
  deallocate(i_piv_eigen)

  deallocate(i_piv_im)
  deallocate(propinv)
  deallocate(prod)
  deallocate(work_im)
  deallocate(gk_left)
  deallocate(gk_right)

  deallocate(diag_scale)

  if (collision_flag == 1) then
     deallocate(nu_op)
     deallocate(i_piv_imk)
     deallocate(propinvk)
     deallocate(prodk)
     deallocate(work_imk)
     deallocate(gk_leftk)
     deallocate(gk_rightk)
  endif
  !-------------------------------------------------------------------

end subroutine gyro_fieldeigen
