subroutine expromake_set

  use expromake_globals
  use EXPRO_interface

  implicit none

  integer :: i,j,is
  real :: x
  real, dimension(:), allocatable :: chi_t
  real, dimension(EXPRO_n_exp) :: a

  if(set_exm_b_ref == 1) then
     EXPRO_b_ref = exm_b_ref
  endif

  if(set_exm_arho == 1) then
     EXPRO_arho = exm_arho
  endif

  if(set_exm_rho == 1) then
     do i=1,EXPRO_n_exp
        x = real(i-1)/(EXPRO_n_exp-1) 
        EXPRO_rho(i)  = x
     enddo
  endif

  if(set_exm_rmin == 1) then
     do i=1,EXPRO_n_exp
        x = real(i-1)/(EXPRO_n_exp-1) 
        EXPRO_rmin(i)  = exm_a_meters*x
     enddo
  endif

  if(set_exm_rmaj == 1) then
     EXPRO_rmaj = exm_rmaj
  endif

  if(set_exm_q == 1) then
     EXPRO_q = exm_q
  endif

  if(set_exm_kappa == 1) then
     EXPRO_kappa = exm_kappa
  endif

  if(set_exm_delta == 1) then
     EXPRO_delta = exm_delta
  endif

  if(set_exm_z_eff == 1) then
     EXPRO_z_eff = exm_z_eff
  endif

  if(set_exm_w0 == 1) then
     EXPRO_w0 = exm_w0
  endif

  if(set_exm_flow_mom == 1) then
     EXPRO_flow_mom = exm_flow_mom
  endif

  if(set_exm_pow_e == 1) then
     EXPRO_pow_e = exm_pow_e
  endif

  if(set_exm_pow_i == 1) then
     EXPRO_pow_i = exm_pow_i
  endif

  if(set_exm_pow_ei == 1) then
     EXPRO_pow_ei = exm_pow_ei
  endif

  if(set_exm_zeta == 1) then
     EXPRO_zeta = exm_zeta
  endif

  if(set_exm_flow_beam == 1) then
     EXPRO_flow_beam = exm_flow_beam
  endif

  if(set_exm_flow_wall == 1) then
     EXPRO_flow_wall = exm_flow_wall
  endif

  if(set_exm_zmag == 1) then
     EXPRO_zmag = exm_zmag
  endif

  if(set_exm_ptot == 1) then
     EXPRO_ptot = exm_ptot
  endif

  if(set_exm_polflux == 1) then
     ! d psi = 1/q * d chi_t
     allocate(chi_t(EXPRO_n_exp))
     do i=1,EXPRO_n_exp
        chi_t(i) = 0.5 * EXPRO_b_ref * (EXPRO_arho*EXPRO_rho(i))**2
     enddo
     EXPRO_poloidalfluxover2pi(1) = 0.0
     do i=2,EXPRO_n_exp
        EXPRO_poloidalfluxover2pi(i) = EXPRO_poloidalfluxover2pi(i-1) &
             + 1.0/(0.5*(EXPRO_q(i) + EXPRO_q(i+1))) &
             * (chi_t(i) - chi_t(i-1))
     enddo
     deallocate(chi_t)
  endif

  !-----------------------------------------------------------------------
  ! Set electron and ion densities
  !
  if (set_exm_ne > 0) then
     select case (exm_ne_model)
     case (1)
        ! Constant density
        EXPRO_ne(:) = exm_ne_axis
     case (2)
        ! Fixed gradient length
        EXPRO_ne(:) = exm_ne_axis &
             * exp(-EXPRO_rmin(:)/EXPRO_rmin(EXPRO_n_exp)*exm_alne) 
        case (3)
           ! Set electron density by quasineutrality
           a(:) = 0.0
           do is=1,nions_max
                 a(:) = a(:)-exm_z(is)*EXPRO_ni(is,:)
           enddo
           EXPRO_ne(:) = a(:)
        case (4)
           ! Import data scaled by factor 
           call expromake_importvec(exm_ne_data,EXPRO_rho(:),EXPRO_ne(:),EXPRO_n_exp)
     case default
        print *, 'ERROR: (expromake) NE_MODEL must be 1-4.'
        stop
     end select
  endif

  do i=1,nions_max

     if (set_exm_ni(i) > 0) then

        select case(exm_ni_model(i))
        case (1)
           ! Constant density
           EXPRO_ni(i,:) = exm_ni_axis(i)
        case (2)
           ! Fixed gradient length
           EXPRO_ni(i,:) = exm_ni_axis(i) &
                * exp(-EXPRO_rmin(:)/EXPRO_rmin(EXPRO_n_exp)*exm_alni(i)) 
        case (3)
           ! Set ion-i density by quasineutrality 
           a(:) = EXPRO_ne(:)
           do is=1,nions_max
              if (is /= i) then
                 a(:) = a(:)-exm_z(is)*EXPRO_ni(is,:)
              endif
           enddo
           EXPRO_ni(i,:) = a(:)/exm_z(i)
        case (4)
           ! Import data scaled by factor 
           call expromake_importvec(exm_ni_data(i),EXPRO_rho(:),EXPRO_ni(i,:),EXPRO_n_exp)
        case default
           print *, 'ERROR: (expromake) NI_MODEL must be 1-4.'
           stop
        end select

     endif
     !-----------------------------------------------------------------------

     !-----------------------------------------------------------------------
     ! Set electron and ion temperatures
     !
     if (set_exm_te == 1) then
        if(exm_te_model == 1) then
           EXPRO_te(:) = exm_te_axis
        endif
        if(exm_te_model == 2) then
           do j=1,EXPRO_n_exp
              EXPRO_te(j) = exm_te_axis &
                   * exp(-EXPRO_rmin(j)/EXPRO_rmin(EXPRO_n_exp)*exm_alte) 
           enddo
        endif
     endif

     if(set_exm_ti(i) == 1) then
        if(exm_ti_model(i) == 1) then
           EXPRO_ti(i,:) = exm_ti_axis(i)
        endif
        if(exm_ti_model(i) == 2) then
           do j=1,EXPRO_n_exp
              EXPRO_ti(i,j) = exm_ti_axis(i) &
                   * exp(-EXPRO_rmin(j)/EXPRO_rmin(EXPRO_n_exp)*exm_alti(i)) 
           enddo
        endif
     endif
     !-----------------------------------------------------------------------

     if(set_exm_vpol(i) == 1) then
        EXPRO_vpol(i,:) = exm_vpol(i)
     endif

     if(set_exm_vtor(i) == 1) then
        EXPRO_vtor(i,:) = exm_vtor(i)
     endif

  enddo

end subroutine expromake_set
