!----------------------------------------------------------------
! prgen_read_omfit.f90
!
! PURPOSE:
!  Read data from OMFIT/GACODE flux-surface mapper 
!----------------------------------------------------------------

subroutine prgen_read_omfit

  use prgen_globals
  use expro

  implicit none

  integer :: npsi,nf,i,j,ip
  real :: fa,fb
  real, dimension(:,:), allocatable :: efit_si,efit_ci
  real, dimension(:), allocatable :: efit_rho,efit_psi,efit_q,efit_p
  real, dimension(:), allocatable :: efit_rmin,efit_rmaj,efit_kappa,efit_zmaj  
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho

  !----------------------------------------------------------------
  open(unit=1,file='out.dim',status='old')
  read(1,*) npsi
  read(1,*) nf
  read(1,*) rcentr
  read(1,*) bcentr
  read(1,*) current
  close(1)

  allocate(efit_si(npsi,nf))
  allocate(efit_ci(npsi,nf))
  allocate(efit_rmin(npsi))
  allocate(efit_rmaj(npsi))
  allocate(efit_kappa(npsi))
  allocate(efit_zmaj(npsi))
  allocate(efit_psi(npsi))
  allocate(efit_q(npsi))
  allocate(efit_p(npsi))
  allocate(efit_rho(npsi))

  open(unit=1,file='out.data',status='old',access='stream')
  read(1) efit_psi
  read(1) efit_q
  read(1) efit_p
  read(1) efit_si
  read(1) efit_ci
  read(1) efit_rmin
  read(1) efit_rmaj
  read(1) efit_kappa
  read(1) efit_zmaj
  close(1)

  efit_psi  = efit_psi-efit_psi(1)
  dpsi_data = dpsi(nx)
  dpsi_efit = efit_psi(npsi)

  ! Get rho on efit mesh (so have psi,rho,q)
  call prgen_get_chi(npsi,efit_q,efit_psi,efit_rho,torfluxa)

  if (format_type /= 3) then
     ! We have rho on statefile grid
     call cub_spline(efit_rho,efit_psi,npsi,rho,dpsi,nx)
  else
     ! We have psinorm on statefile grid
     dpsi = dpsi*dpsi_efit
     call cub_spline(efit_psi,efit_q,npsi,dpsi,q,nx)
     call prgen_get_chi(nx,q,dpsi,rho,torfluxa)
  endif

  call cub_spline(efit_rho,efit_rmin,npsi,rho,rmin,nx)
  call cub_spline(efit_psi,efit_q,npsi,dpsi,q,nx)
  call cub_spline(efit_psi,efit_p,npsi,dpsi,p_tot,nx)
  call cub_spline(efit_psi,efit_zmaj,npsi,dpsi,zmag,nx)
  call cub_spline(efit_psi,efit_rmaj,npsi,dpsi,rmaj,nx)
  call cub_spline(efit_psi,efit_kappa,npsi,dpsi,kappa,nx)
  call cub_spline(efit_psi,efit_si(:,2),npsi,dpsi,delta,nx) ; delta = sin(delta)
  call cub_spline(efit_psi,efit_si(:,3),npsi,dpsi,zeta,nx) ; zeta = -zeta

  ! New shape coefficients
  call cub_spline(efit_psi,efit_si(:,4),npsi,dpsi,shape_sin3,nx)
  call cub_spline(efit_psi,efit_ci(:,1),npsi,dpsi,shape_cos0,nx)
  call cub_spline(efit_psi,efit_ci(:,2),npsi,dpsi,shape_cos1,nx)
  call cub_spline(efit_psi,efit_ci(:,3),npsi,dpsi,shape_cos2,nx)
  call cub_spline(efit_psi,efit_ci(:,4),npsi,dpsi,shape_cos3,nx)

  !==============================================================================

  if (nfourier > 0) then

     ! Old Fourier rpresentation
     allocate(g3vec(npsi,0:nfourier,4))
     allocate(g3rho(nx,0:nfourier,4))

     open(unit=1,file='fluxfit.geo',status='old',access='stream')
     do i=1,4
        read(1) g3vec(:,:,i)
     enddo
     close(1)

     ! Explicitly set rmin=0 at origin
     g3vec(1,1:nfourier,:) = 0.0

     ! Extrapolate centers (R0,Z0) to origin
     do i=1,4
        call bound_extrap(fa,fb,g3vec(:,0,i),efit_psi,npsi)
        g3vec(1,0,i) = fa
     enddo

     ! Map results onto poloidal flux (dpsi) grid:
     ! NOTE: need sqrt here to get sensible behaviour as r -> 0.
     do i=1,4
        do ip=0,nfourier
           call cub_spline(sqrt(efit_psi),g3vec(:,ip,i),npsi,sqrt(dpsi),g3rho(:,ip,i),nx)
        enddo
     enddo

     ! Explicitly set rmin=0 at origin
     g3rho(:,1:nfourier,1) = 0.0

     open(unit=1,file='input.gacode.geo',status='replace')
     write(1,'(a)') '# input.gacode.geo'
     write(1,'(a)') '#'
     write(1,'(a)') '# File format:'
     write(1,'(a)') '#-------------------'
     write(1,'(a)') '# nfourier'
     write(1,'(a)') '# a[4,0:nfourier,nx]'  
     write(1,'(a)') '#-------------------'
     write(1,'(a)') '#'
     write(1,'(a)') '# NOTE: nx=EXPRO_n_exp is defined in input.gacode'
     write(1,*) nfourier
     do j=1,nx
        do ip=0,nfourier
           do i=1,4
              write(1,'(1pe20.13)') g3rho(j,ip,i)
           enddo
        enddo
     enddo
     close(1)
     !----------------------------------------------------

     ! Clean up
     deallocate(g3vec)
     deallocate(g3rho)

  endif

end subroutine prgen_read_omfit
