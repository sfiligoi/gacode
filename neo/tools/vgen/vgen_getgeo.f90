subroutine vgen_getgeo
  use vgen_globals
  use EXPRO_interface
  use GEO_interface

  implicit none

  integer :: i, j
  real :: r_min
  !integer :: n_theta=2
  !real, dimension(2) :: theta=(/ 0.0, -1.8325957/)
  integer :: n_theta=20
  real, dimension(:), allocatable :: theta

  allocate(theta(n_theta))
  do j=1,n_theta
     theta(j) = -pi + (j-1)*2*pi/n_theta
  enddo

  GEO_nfourier_in = EXPRO_nfourier
  GEO_signb_in    = EXPRO_ctrl_signb
  call GEO_alloc(1)

  open(unit=1,file='vgen_geo.out',status='replace')
  r_min = EXPRO_rmin(EXPRO_n_exp)

  do i=2,EXPRO_n_exp
     
     ! Parameters to be passed to GEO library   
     !
     ! NOTE: dp/dr set to zero without loss of generality.
     ! 
     GEO_rmin_in      = EXPRO_rmin(i)/r_min
     GEO_rmaj_in      = EXPRO_rmaj(i)/r_min
     GEO_drmaj_in     = EXPRO_drmaj(i)
     GEO_zmag_in      = EXPRO_zmag(i)/r_min
     GEO_dzmag_in     = EXPRO_dzmag(i)
     GEO_q_in         = EXPRO_q(i)
     GEO_s_in         = EXPRO_s(i)
     GEO_kappa_in     = EXPRO_kappa(i)
     GEO_s_kappa_in   = EXPRO_skappa(i)
     GEO_delta_in     = EXPRO_delta(i)
     GEO_s_delta_in   = EXPRO_sdelta(i)
     GEO_zeta_in      = EXPRO_zeta(i)
     GEO_s_zeta_in    = EXPRO_szeta(i)
     GEO_beta_star_in = 0.0
     !
     if (EXPRO_ctrl_numeq_flag == 0) then
        ! Call GEO with model shape
        GEO_model_in = 0
        call GEO_do()
     else
        ! Call GEO with general (numerical) shape
        GEO_model_in = 1
        GEO_fourier_in(1:4,:) = EXPRO_geo(:,:,i)/r_min
        GEO_fourier_in(5:8,:) = EXPRO_dgeo(:,:,i)
        call GEO_do()
     endif
    
     write(1,'(e16.8,$)') EXPRO_rho(i)
     do j=1, n_theta
        call GEO_interp(theta(j))
        write(1,'(e16.8,$)') theta(j)
        write(1,'(e16.8,$)') GEO_b  * EXPRO_bunit(i)
        write(1,'(e16.8,$)') GEO_bp * EXPRO_bunit(i)
        write(1,'(e16.8,$)') GEO_bt * EXPRO_bunit(i)
        write(1,'(e16.8,$)') GEO_bigr * r_min
     enddo
     write(1,*)
     
  enddo
  
  ! Clean up
  close(1)
  deallocate(theta)
  call GEO_alloc(0)

end subroutine vgen_getgeo
