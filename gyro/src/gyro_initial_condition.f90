!-----------------------------------------------------
! gyro_initial_condition.f90
!
! PURPOSE:
!  Set initial value of h.
!
! NOTES:
!  This is controlled by input parameters
!  AMP_PHI_0     -> amp_0
!  AMP_PHI_N     -> amp_n
!  AMP_PHI_STUDY -> amp_study (and n_study) 
!-------------------------------------------------------

subroutine gyro_initial_condition

  use gyro_globals
  use gyro_pointers
  use math_constants

  !--------------------------------- 
  implicit none
  !
  real :: x
  real :: theta_0
  real :: scale
  real :: a_spec(n_kinetic)
  real :: hs(n_x)
  real :: as_0
  real :: as_n
  !----------------------------------

  !-------------------------------------------------------------------
  h   = (0.0,0.0)
  h_0 = (0.0,0.0)
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Normalization should keep mean-square potential 
  ! invariant.  The overall normalization has been 
  ! retained for backward compatibility with 16 modes.
  !
  scale = (rhos_norm/0.00357)*sqrt(16.0/n_n)
  !
  as_0 = amp_0*scale
  as_n = amp_n*scale
  !
  ! Next, we need to correct for possible n_study reset
  !
  if (n_study /= 0) then
     if (n_1(in_1) == n(n_study)) then
        if (n_1(in_1) == 0) then
           as_0 = amp_study*scale
        else
           as_n = amp_study*scale
        endif
     endif
  endif
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! We need the following normalization factor to ensure that
  ! multiple identical species (with densities which sum to 
  ! unity) give results identical to a single species.  The 
  ! point here is to weight h by the fractional density.
  ! 
  a_spec(1:n_kinetic) = den_s(1:n_kinetic,ir_norm)/den_s(n_spec,ir_norm)
  !-------------------------------------------------------------------

  select case (boundary_method)

  case (1)

     ! PERIODIC

     do i=1,n_x
        x = (i-1.0)/n_x 
        hs(i) = sin(2*pi*x)
     enddo

  case (2)

     ! NONPERIODIC

     do i=1,n_x
        x = (i-1.0)/(n_x-1) 
        hs(i) = sin(2*pi*x)
     enddo

  end select


  ! Generate ICs

  p_nek_loc = 0

  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     if (n_1(in_1) /= 0) then

        !==========
        ! n > 0
        !==========

        if (amp_n > 0.0) then

           ! Hump in theta

           do m=1,n_stack
              do i=1,n_x
                 theta_0 = theta_t(i,k,m)
                 h(m,i,p_nek_loc,:) = a_spec(:)*as_n*cos(theta_0/2.0)**2
              enddo ! i
           enddo ! m

        else

           ! Sharp hump in theta

           do m=1,n_stack
              do i=1,n_x
                 theta_0 = theta_t(i,k,m)
                 h(m,i,p_nek_loc,:) = a_spec(:)*abs(as_n)*cos(theta_0/2.0)**8
              enddo ! i
           enddo ! m

        endif

     else

        !==========
        ! n = 0
        !==========

        ! n=0 modes need zero radial average

        if (amp_0 > 0.0) then

           ! Sine wave in radius 

           do m=1,n_stack
              do i=1,n_x
                 h(m,i,p_nek_loc,:) = a_spec(:)*as_0*hs(i)
              enddo ! i
           enddo ! m

        else

           ! Sine^9 in radius

           do m=1,n_stack
              do i=1,n_x
                 h(m,i,p_nek_loc,:) = a_spec(:)*as_0*hs(i)**9
              enddo ! i
           enddo ! m

       endif

     endif

  enddo ! p_nek


  ! Compute fields from distribution
  call gyro_field_solve_explicit

  ! Generate interpolation of fields suitable for plotting
  call get_field_plot

  ! Compute gyro_u needed for gyro_moments_plot
  call gyro_field_interpolation

  ! Generate data for moments plot
  call gyro_moments_plot

  if (field_r0_flag == 1) call get_field_r0_plot 

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_initial_condition done]' 
  endif
 
end subroutine gyro_initial_condition
