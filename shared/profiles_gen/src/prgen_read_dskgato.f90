!---------------------------------------------------------
! prgen_read_dskgato.f90
!
! PURPOSE:
!  Read DSKGATO format.
!
!----------------------------------------------------------

subroutine prgen_read_dskgato

  use prgen_globals

  implicit none

  character (len=1) :: a
  integer :: i,ip
  integer :: neqtyp
  integer :: ntht,neqsym
  real :: fa,fb
  real :: dummy(4)
  real, dimension(:), allocatable :: psi
  real, dimension(:), allocatable :: pdum
  real, dimension(:), allocatable :: tdum
  real, dimension(:,:), allocatable :: xs
  real, dimension(:,:), allocatable :: zs
  real, dimension(:), allocatable :: sqdpsi

  real, dimension(:,:), allocatable :: gvec
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho

  character (len=70) :: cdum

  !---------------------------------------------------
  ! 
  neqtyp=1
  open(unit=1,file='grid.dat',status='old')
  if (neqtyp == 0) then
     read(1,'(a)') cdum
     read(1,'(a)') cdum
     read(1,*) nsurf,ntht,neqsym
  else
     read(1,*) nsurf,ntht
     neqsym = 1
  endif

  if (neqsym == 0) narc = ntht
  if (neqsym == 1) narc = 2*(ntht-1) + 1

  read (1,10) dummy(1:4)

  if (neqtyp == 0) read(1,10) dummy(1:2)
  if (neqtyp == 1) read(1,10) dummy(1:3)

  allocate(tdum(ntht))
  allocate(psi(nsurf))
  allocate(pdum(nsurf))
  allocate(xs(nsurf,narc))
  allocate(zs(nsurf,narc))

  read(1,10) psi(:)
  read(1,10) pdum(:) ! fval
  read(1,10) pdum(:) ! ffprime
  read(1,10) pdum(:) ! sp
  read(1,10) pdum(:) ! pprime
  read(1,10) pdum(:) ! qsfin
  read(1,10) pdum(:) ! nel

  read(1,10) tdum(:) ! seqdpdr
  read(1,10) tdum(:) ! seqdpdz

  read(1,10) xs(:,1:ntht)
  read(1,10) zs(:,1:ntht)
  close(1)

  zs(:,1)    = 0.0
  zs(:,ntht) = 0.0

  do i=ntht+1,narc
     ip = narc+1-i
     xs(:,i) = +xs(:,ip)
     zs(:,i) = -zs(:,ip)
  enddo
  !----------------------------------------------------

  !----------------------------------------------------------------
  ! Correct flux variation in profile data (ITERDB, etc) 
  !
  if (format_type == 0 .or. format_type == 7) then
     ! Case 1: raw gfile mode ; dpsi is undefined
     dpsi_data = psi(nsurf)
     dpsi_efit = psi(nsurf)
     do i=1,nx
        ! Shrink slightly for safety
        dpsi(i) = (i-1)*dpsi_efit/(nx-1+1e-12)
     enddo
  else
     ! Case 2: typical case 
     dpsi_data = dpsi(nx)
     dpsi_efit = psi(nsurf)

     ! Ensure max(dpsi) = max(gato_psi) 
     dpsi(:)  = dpsi(:)*dpsi_efit/dpsi_data
     ! Extra insurance against roundoff
     dpsi(nx) = dpsi_efit 
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Get Miller-style (model shape) coefficients using 
  ! fluxfit routines:
  !
  ! 1  r (cm)
  ! 2  z_mag (cm)
  ! 3  R0 (cm)
  ! 4  kappa
  ! 5  delta
  ! 6  zeta
  !
  call fluxfit_driver(1,1,nsurf,narc,xs,zs,verbose_flag)
  !
  allocate(gvec(6,nsurf))
  !
  open(unit=1,file='fluxfit.profile',status='old')
  read(1,*) a
  read(1,*) gvec(:,:)
  close(1)

  allocate(sqdpsi(nx))
  sqdpsi = sqrt(dpsi)

  ! Explicitly set rmin=0 at origin
  gvec(1,1) = 0.0

  ! Use extrapolation to get values of shape parameters at origin
  do i=2,6
     call bound_extrap(fa,fb,gvec(i,:),psi,nsurf)
     gvec(i,1) = fa
  enddo

  ! Map shape coefficients onto poloidal flux (dpsi) grid:
  call cub_spline(sqrt(psi),gvec(1,:),nsurf,sqdpsi,rmin,nx)
  call cub_spline(sqrt(psi),gvec(2,:),nsurf,sqdpsi,zmag,nx)
  call cub_spline(sqrt(psi),gvec(3,:),nsurf,sqdpsi,rmaj,nx)
  call cub_spline(sqrt(psi),gvec(4,:),nsurf,sqdpsi,kappa,nx)
  call cub_spline(sqrt(psi),gvec(5,:),nsurf,sqdpsi,delta,nx)
  call cub_spline(sqrt(psi),gvec(6,:),nsurf,sqdpsi,zeta,nx)

  ! Explicitly set rmin=0 at origin
  rmin(1) = 0.0

  deallocate(gvec)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Get general geometry coefficients
  !
  call fluxfit_driver(2,nfourier,nsurf,narc,xs,zs,verbose_flag)

  allocate(g3vec(4,0:nfourier,nsurf))
  allocate(g3rho(4,0:nfourier,nx))

  open(unit=1,file='fluxfit.geo',status='old')
  read(1,*) ip
  read(1,*) g3vec(:,:,:)
  close(1)

  ! Explicitly set rmin=0 at origin
  g3vec(:,1:nfourier,1) = 0.0

  ! Extrapolate centers (R0,Z0) to origin
  do i=1,4
     call bound_extrap(fa,fb,g3vec(i,0,:),psi,nsurf)
     g3vec(i,0,1) = fa
  enddo

  ! Map Fourier coefficients onto poloidal flux (dpsi) grid 
  do i=1,4
     do ip=0,nfourier
        call cub_spline(sqrt(psi),g3vec(i,ip,:),nsurf,sqdpsi,g3rho(i,ip,:),nx)
     enddo
  enddo

  ! Explicitly set rmin=0 at origin
  g3rho(:,1:nfourier,1) = 0.0

  open(unit=1,file='input.profiles.geo',status='replace')
  write(1,'(a)') '# input.profiles.geo'
  write(1,'(a)') '#'
  write(1,'(a)') '# See https://fusion.gat.com/theory/input.profiles.geo for complete documentation.'
  write(1,'(a)') '#'
  write(1,'(a)') '# File format:'
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '# nfourier'
  write(1,'(a)') '# a[4,0:nfourier,nx]'  
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '#'
  write(1,'(a)') '# NOTE: nx=EXPRO_n_exp is defined in input.profiles'
  write(1,*) nfourier
  write(1,'(1pe20.13)') g3rho(:,:,:)
  close(1)
  !----------------------------------------------------------------

  ! Cleanup
  deallocate(g3vec)
  deallocate(g3rho)
  deallocate(xs)
  deallocate(zs)
  deallocate(psi)
  deallocate(tdum)
  deallocate(pdum)

10 format(1p4e19.12)

end subroutine prgen_read_dskgato
