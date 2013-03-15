!---------------------------------------------------------
! prgen_read_omfit.f90
!
! PURPOSE:
!  Read data from OMFIT flux-surface mapper and fit to 
!  geometry coefficients.
!----------------------------------------------------------

subroutine prgen_read_omfit

  use prgen_globals

  implicit none

  integer :: i,i1,i2
  integer :: ip
  integer :: ntot
  character (len=1) :: a
  integer, dimension(:), allocatable :: narcv
  real, dimension(:), allocatable :: r_raw,z_raw,l_raw
  real, dimension(:,:), allocatable :: r2,z2
  real, dimension(:), allocatable :: psi
  real, dimension(:,:), allocatable :: gvec
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho

  ! Get the gfile header
  open(unit=1,file='gfile',status='old')
  read(1,'(a)') efit_header
  close(1) 

  ! Read the OMFIT mapper file
  open(unit=1,file='ffile',status='old')

  do i=1,6
     read(1,*) a
  enddo
  read(1,*) nsurf
  allocate(narcv(nsurf))
  read(1,*) narcv(:)

  ntot = sum(narcv)

  narc = narcv(1)

  allocate(psi(nsurf))
  allocate(r_raw(ntot))
  allocate(z_raw(ntot))
  allocate(l_raw(ntot))

  read(1,*) psi(:) 
  read(1,*) r_raw(:)
  read(1,*) z_raw(:)
  read(1,*) l_raw(:)
  close(1)

  psi(:) = psi(:)-psi(1)

  !----------------------------------------------------
  ! Get Miller-style (model shape) coefficients using 
  ! fluxfit routines:
  !
  ! 1  r (cm)
  ! 2  z_mag (cm)
  ! 3  R0 (cm)
  ! 4  kappa
  ! 5  delta
  ! 6  zeta

  allocate(r2(narc,nsurf))
  allocate(z2(narc,nsurf))
  
  do i=1,nsurf
     if (i == 1) then
        i1 = 1
     else
        i1 = sum(narcv(1:(i-1)))+1
     endif
     i2 = sum(narcv(1:i))
     r2(:,i) = r_raw(i1:i2)
     z2(:,i) = z_raw(i1:i2)
  enddo

  ! When verbose_flag=1, fluxfit will echo lots of data
  call fluxfit_driver(1,1,nsurf,narc,r2,z2,verbose_flag)

  allocate(gvec(6,nsurf))

  open(unit=1,file='fluxfit.profile',status='old')
  read(1,*) a
  read(1,*) gvec(:,:)
  close(1)

  ! Use linear interpolation to get values of 
  ! shape parameters at origin
  gvec(:,1) = (psi(3)*gvec(:,2)-psi(2)*gvec(:,3))/(psi(3)-psi(2))

  ! Manage total flux variation 
  if (format_type == 0 .or. format_type == 7) then
     ! Case 1: raw gfile mode ; dpsi is undefined
     dpsi_data = psi(nsurf)
     dpsi_gato = psi(nsurf)
     do i=1,nx
        dpsi(i) = (i-1)*dpsi_gato/(nx-1)
     enddo
  else
     ! Case 2: typical case 
     dpsi_data = dpsi(nx)
     dpsi_gato = psi(nsurf)

     ! Ensure max(dpsi) = max(gato_psi) 
     dpsi(:)  = dpsi(:)*dpsi_gato/dpsi_data
     ! Extra insurance against roundoff
     dpsi(nx) = dpsi_gato 
  endif

  ! Explicitly set rmin=0 at origin
  gvec(1,1) = 0.0

  ! Map shape coefficients onto poloidal flux (dpsi) grid:
  call cub_spline(psi,gvec(1,:)**2,nsurf,dpsi,rmin,nx)
  ! Want to spline r^2 since its linear in psi
  rmin = sqrt(rmin)
  call cub_spline(psi,gvec(2,:),nsurf,dpsi,zmag,nx)
  call cub_spline(psi,gvec(3,:),nsurf,dpsi,rmaj,nx)
  call cub_spline(psi,gvec(4,:),nsurf,dpsi,kappa,nx)
  call cub_spline(psi,gvec(5,:),nsurf,dpsi,delta,nx)
  call cub_spline(psi,gvec(6,:),nsurf,dpsi,zeta,nx)

  ! Explicitly set rmin=0 at origin
  rmin(1) = 0.0

  deallocate(gvec)
  !----------------------------------------------------

  !----------------------------------------------------
  ! Get general geometry coefficients
  !
  ! When verbose_flag=1, fluxfit will echo lots of data
  call fluxfit_driver(2,nfourier,nsurf,narc,r2,z2,verbose_flag)

  allocate(g3vec(4,0:nfourier,nsurf))
  allocate(g3rho(4,0:nfourier,nx))

  open(unit=1,file='fluxfit.geo',status='old')
  read(1,*) ip
  read(1,*) g3vec(:,:,:)
  close(1)

  ! Use linear interpolation to get values of 
  ! (general) shape parameters at origin
  g3vec(:,:,1) = (psi(3)*g3vec(:,:,2)-psi(2)*g3vec(:,:,3))/(psi(3)-psi(2))

  ! Map results onto poloidal flux (dpsi) grid:
  do i=1,4
     do ip=0,nfourier
        call cub_spline(psi,g3vec(i,ip,:),nsurf,dpsi,g3rho(i,ip,:),nx)
     enddo
  enddo

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
  !----------------------------------------------------

  ! Clean up
  deallocate(g3vec)
  deallocate(g3rho)
  deallocate(r2)
  deallocate(z2)
  deallocate(r_raw)
  deallocate(z_raw)
  deallocate(l_raw)
  deallocate(narcv)
  deallocate(psi)

end subroutine prgen_read_omfit
