!------------------------------------------------------------------
! gyro_diffusivity.f90
!
! PURPOSE:
!  Renormalize fluxes computed in gyro_nonlinear_flux to 
!  obtain diffusivities.  Also compute averaged diffusivities.
!------------------------------------------------------------------

subroutine gyro_diffusivity

  use mpi
  use gyro_globals
  use math_constants

  !-------------------------------------------------------------
  implicit none
  !
  real, dimension(n_kinetic,n_field,n_moment,n_x) :: diff_i_loc
  real, dimension(n_kinetic,n_field,n_moment,n_x) :: diff_i_loc_trapped
  real, dimension(n_kinetic,n_field,n_moment) :: diff_loc
  real, dimension(n_kinetic,n_field,n_moment) :: diff_loc_trapped
  real :: diff_norm
  !-------------------------------------------------------------


  !-----------------------------------------------------------------
  ! Define factors useful for quasilinear normalizations
  !
  call gyro_phi_kp_squared
  call get_g_squared

  ! phi_squared_QL_n = <|phi_n(x)|^2>/rho_s^2   
  ! g_squared_QL_n(1) = <|g_n(x)|^2>/rho_s^2  
  ! g_squared_QL_n(2) = <|phi_n(x)|^2>/rho_s^2     
  ! g_squared_QL_n(3) = sqrt[<|g_n(x)|^2><|phi_n(x)|^2>]/rho_s^2  
  !
  ! where <.> is the radial/flux-surface average.

  phi_squared_QL_n(:) = 0.0
  g_squared_QL_n(:) = 0.0

  if (n_1(in_1) /= 0) then
     do i=1,n_x

        ! Factor of 2 for both signs of n, and remove 
        ! krho_i(in_1,i) factor from phi_squared.

        phi_squared_QL_n(1) = phi_squared_QL_n(1) + &
             2.0*phi_squared(i)/krho_i(in_1,i)/rhos_norm**2

        g_squared_QL_n(1:2) = g_squared_QL_n(1:2) + &
             2.0*g_squared(1:2,i)/rhos_norm**2

     enddo
     phi_squared_QL_n(1) = phi_squared_QL_n(1)/n_x
     g_squared_QL_n(1:2) = g_squared_QL_n(1:2)/n_x

  endif

  g_squared_QL_n(3) = sqrt(g_squared_QL_n(1)*g_squared_QL_n(2)) 
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  ! Renomalize linear diffusivities
  ! 
  if (nonlinear_flag == 0 .and. n_1(in_1) /= 0 .and. lindiff_method > 2) then

     select case (lindiff_method)

     case (3)

        ! Gamma_n/(k_theta rho_s |phi_n|^2)
        ! Q_n/(k_theta rho_s |phi_n|^2)

        do i=1,n_x

           diff_i_loc(:,:,:,i) = (nonlinear_flux_passing(i,:,:,1:n_moment)+ &
                nonlinear_flux_trapped(i,:,:,1:n_moment))/sum(phi_squared(:)/n_x)

           diff_i_loc_trapped(:,:,:,i) = nonlinear_flux_trapped(i,:,:,1:n_moment)/ &
                sum(phi_squared(:)/n_x)

        enddo ! i

     case (4)

        ! rho_s^2 D_n/<|phi_n|^2>
        ! rho_s^2 Chi_n/<|phi_n|^2>

        diff_i_loc(:,:,:,:) = &
             diff_i_loc(:,:,:,:)/phi_squared_QL_n(1)
        diff_i_loc_trapped(:,:,:,:) = &
             diff_i_loc_trapped(:,:,:,:)/phi_squared_QL_n(1)

     case (5)

        ! rho_s^2 D_n/<|phi_n|^2>
        ! rho_s^2 Chi_n/<|phi_n|^2>

        do i=1,n_x
           diff_i_loc(:,:,:,i) = diff_i_loc(:,:,:,i)/ &
                (2.0*phi_squared(i)/krho_i(in_1,i)/rhos_norm**2)
           diff_i_loc_trapped(:,:,:,i) = diff_i_loc_trapped(:,:,:,i)/ &
                (2.0*phi_squared(i)/krho_i(in_1,i)/rhos_norm**2)
        enddo

     case (6)

        ! rho_s^2 D_n/<|g_n|^2>
        ! rho_s^2 Chi_n/<|g_n|^2>

        diff_i_loc(:,:,:,:) = diff_i_loc(:,:,:,:)/g_squared_QL_n(1)
        diff_i_loc_trapped(:,:,:,:) = diff_i_loc_trapped(:,:,:,:)/g_squared_QL_n(1)

     case (7)

        ! rho_s^2 D_n/<|phi_n|^2>
        ! rho_s^2 Chi_n/<|phi_n|^2>

        diff_i_loc(:,:,:,:) = diff_i_loc(:,:,:,:)/g_squared_QL_n(2)
        diff_i_loc_trapped(:,:,:,:) = diff_i_loc_trapped(:,:,:,:)/g_squared_QL_n(2)

     case (8)

        ! rho_s^2 D_n/<|phi_n||g_n|>
        ! rho_s^2 Chi_n/<|phi_n||g_n|>

        diff_i_loc(:,:,:,:) = diff_i_loc(:,:,:,:)/g_squared_QL_n(3)
        diff_i_loc_trapped(:,:,:,:) = diff_i_loc_trapped(:,:,:,:)/g_squared_QL_n(3)

     end select

  else

     !----------------------------------------------------------
     ! Compute particle/momentum/energy diffusivities:
     !
     ! Below, note that nonlinear_flux = 0 if n=0, and 
     ! also that the factor of 2.0 comes from 2 signs of n.
     !-----------------------------------------------------

     do ix=1,n_field
        do is=1,n_kinetic
           do i=1,n_x

              ! Particle diffusivity (gyroBohm units):
              diff_i_loc(is,ix,1,i) = 2.0*(nonlinear_flux_passing(i,is,ix,1)+ &
                   nonlinear_flux_trapped(i,is,ix,1))/ &
                   (den_s(is,i)*dlnndr_s(is,i))/rhos_norm**2

              ! Trapped particle diffusivity (gyroBohm units):
              diff_i_loc_trapped(is,ix,1,i) = 2.0*nonlinear_flux_trapped(i,is,ix,1)/ &
                   (den_s(is,i)*dlnndr_s(is,i))/rhos_norm**2

              ! Energy diffusivity (gyroBohm units):
              diff_i_loc(is,ix,2,i) = 2.0*(nonlinear_flux_passing(i,is,ix,2)+ &
                   nonlinear_flux_trapped(i,is,ix,2))/ &
                   (den_s(is,i)*tem_s(is,i)*dlntdr_s(is,i))/rhos_norm**2

              ! Trapped particle energy diffusivity (gyroBohm units):
              diff_i_loc_trapped(is,ix,2,i) = 2.0*nonlinear_flux_trapped(i,is,ix,2)/ &
                   (den_s(is,i)*tem_s(is,i)*dlntdr_s(is,i))/rhos_norm**2

           enddo ! i
        enddo ! is
     enddo ! ix
     !----------------------------------------------------------

  endif
  !-----------------------------------------------------------------

  !----------------------------------------------------------
  ! Get (crude) averages over radius.  Do not include the 
  ! region inside the explicit damper.
  !
  diff_loc(:,:,:)         = 0.0
  diff_loc_trapped(:,:,:) = 0.0
  !
  do i=1+n_explicit_damp,n_x-n_explicit_damp
     diff_loc(:,:,:) = diff_loc(:,:,:)+diff_i_loc(:,:,:,i)
     diff_loc_trapped(:,:,:) = diff_loc_trapped(:,:,:)+diff_i_loc_trapped(:,:,:,i)
  enddo ! i
  diff_loc(:,:,:) = diff_loc(:,:,:)/(n_x-2*n_explicit_damp)
  diff_loc_trapped(:,:,:) = diff_loc_trapped(:,:,:)/(n_x-2*n_explicit_damp)
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! Store diff_loc as diff_n for n-dependent diffusivity:
  !
  diff_n(:,:,:) = diff_loc(:,:,:)
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! At this point, all "loc" quantities are local to 
  ! subgroups.  We need to sum these over n (i.e.,
  ! over NEW_COMM_2):
  !
  call MPI_ALLREDUCE(diff_loc, &
       diff, &
       size(diff), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call MPI_ALLREDUCE(diff_loc_trapped, &
       diff_trapped, &
       size(diff_trapped), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call MPI_ALLREDUCE(diff_i_loc, &
       diff_i, &
       size(diff_i), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call MPI_ALLREDUCE(diff_i_loc_trapped, &
       diff_i_trapped, &
       size(diff_i_trapped), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)
  !----------------------------------------------------------

  if (nonlinear_flag == 0 .and. n_1(in_1) /= 0) then
     if (lindiff_method == 2) then

        ! [D_a, Chi_a] / Chi_i normalization for linear diffusivities.

        diff_norm = diff(1,1,2)

        diff(:,:,:)             = diff(:,:,:)/diff_norm
        diff_trapped(:,:,:)     = diff_trapped(:,:,:)/diff_norm
        diff_i(:,:,:,:)         = diff_i(:,:,:,:)/diff_norm
        diff_i_trapped(:,:,:,:) = diff_i_trapped(:,:,:,:)/diff_norm

     endif
  endif


  if (i_proc == 0) then
     if (debug_flag == 1) print *,'[gyro_diffusivity called]'
  endif

end subroutine gyro_diffusivity
