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

  allocate(efit_si(npsi,0:nf-1))
  allocate(efit_ci(npsi,0:nf-1))
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
  dpsi_efit = efit_psi(npsi)
  dpsi_data = dpsi(nx)

  ! Get rho on efit mesh (so have psi,rho,q)
  call prgen_get_chi(npsi,efit_q,efit_psi,efit_rho,torfluxa)

  select case (format_type)

  case (3,4,5)
     ! Statefile (required) supplies psi_norm (dpsi)
     ! PEQDSK=3, CORSICA=5, GENF=9
     dpsi = dpsi*dpsi_efit
     ! Get q
     call cub_spline(efit_psi,efit_q,npsi,dpsi,q,nx)
     ! Get rho
     call prgen_get_chi(nx,q,dpsi,rho,torfluxa)
  case (0,1,2,6,7,8)
     ! Statefile supplies rho
     ! Get dpsi
     call cub_spline(efit_rho,efit_psi,npsi,rho,dpsi,nx)
  end select

  call bound_interp(efit_rho,efit_rmin,npsi,rho,rmin,nx)
  call bound_interp(efit_psi,efit_q,npsi,dpsi,q,nx)
  call bound_interp(efit_psi,efit_p,npsi,dpsi,p_tot,nx)
  call bound_interp(efit_rho,efit_zmaj,npsi,rho,zmag,nx)
  call bound_interp(efit_rho,efit_rmaj,npsi,rho,rmaj,nx)
  call bound_interp(efit_rho,efit_kappa,npsi,rho,kappa,nx)

  ! Miller extended harmonic shape coefficients
  call bound_interp(efit_rho,efit_ci(:,0),npsi,rho,shape_cos(0,:),nx)

  if (nf-1 > 6) then
     print '(a)','ERROR: (prgen_read_omfit) max nharm is 6'
     stop
  endif
  
  do i=1,nf-1
     call bound_interp(efit_rho,efit_si(:,i),npsi,rho,shape_sin(i,:),nx)
     call bound_interp(efit_rho,efit_ci(:,i),npsi,rho,shape_cos(i,:),nx)
  enddo
  
  delta = sin(shape_sin(1,:))
  zeta  = -shape_sin(2,:)
  !==============================================================================

  if (nfourier > 0) then

     ! Legacy direct Fourier representation
     allocate(g3vec(npsi,0:nfourier,4))
     allocate(g3rho(nx,0:nfourier,4))

     open(unit=1,file='fluxfit.geo',status='old',access='stream')
     do i=1,4
        read(1) g3vec(:,:,i)
     enddo
     close(1)

     ! Map results onto poloidal flux (dpsi) grid:
     ! NOTE: need sqrt here to get sensible behaviour as r -> 0.
     do i=1,4
        do ip=0,nfourier
           call bound_interp(efit_rho,g3vec(:,ip,i),npsi,rho,g3rho(:,ip,i),nx)
        enddo
     enddo

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

  deallocate(efit_si)
  deallocate(efit_ci)
  deallocate(efit_rmin)
  deallocate(efit_rmaj)
  deallocate(efit_kappa)
  deallocate(efit_zmaj)
  deallocate(efit_psi)
  deallocate(efit_q)
  deallocate(efit_p)
  deallocate(efit_rho)

end subroutine prgen_read_omfit
