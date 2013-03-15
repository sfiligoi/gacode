!---------------------------------------------------------
! prgen_read_gato.f90
!
! PURPOSE:
!  Read data from GATO (eqgta, grid.dat, o1gta):
!
!  - psi (poloidal flux over 2pi)
!  - V(psi) (volume)
!  - dV(psi)/dpsi
!  - R(psi,arclength)
!  - Z(psi,arclength)
!  - q
!----------------------------------------------------------

subroutine prgen_read_gato

  use prgen_globals

  implicit none

  integer :: i
  integer :: ip
  real, dimension(:), allocatable :: gato_psi  
  real, dimension(:), allocatable :: gato_volume
  real, dimension(:), allocatable :: gato_dvoldpsi
  real, dimension(:), allocatable :: gato_q
  real, dimension(:,:), allocatable :: gato_bigr
  real, dimension(:,:), allocatable :: gato_bigz
  real, dimension(:,:), allocatable :: gvec
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho

  character (len=70) :: cdum

  character (len=10) :: a
  character (len=10) :: b
  character (len=10) :: c

  !--------------------------------------------------
  ! Read gfile header info
  !
  open(unit=1,file='eqgta',status='old')
  read(1,'(a)') efit_header
  close(1) 
  !--------------------------------------------------

  !---------------------------------------------------
  ! Read flux-surface data in GATO's grid.dat file
  !
  ! nsurf = number of finite flux surfaces
  ! 
  open(unit=1,file='grid.dat',status='old')
  read(1,'(a)') cdum
  read(1,*) nsurf,narc

  allocate(gato_psi(0:nsurf))
  allocate(gato_volume(0:nsurf))
  allocate(gato_dvoldpsi(0:nsurf))
  allocate(gato_bigr(narc,nsurf))
  allocate(gato_bigz(narc,nsurf))

  read(1,'(a)') cdum
  read(1,'(a)') cdum
  read(1,*) gato_psi(:) 
  read(1,'(a)') cdum
  read(1,*) gato_volume(:)
  read(1,'(a)') cdum
  read(1,*) gato_dvoldpsi(:)
  read(1,'(a)') cdum
  read(1,*) gato_bigr(:,:)
  read(1,'(a)') cdum
  read(1,*) gato_bigz(:,:)

  close(1)
  !----------------------------------------------------

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

  ! When verbose_flag=1, fluxfit will echo lots of data
  call fluxfit_driver(1,1,nsurf,narc,gato_bigr,gato_bigz,verbose_flag)
  allocate(gvec(6,0:nsurf))

  open(unit=1,file='fluxfit.profile',status='old')
  read(1,*) cdum
  read(1,*) gvec(:,1:nsurf)
  close(1)

  ! Use linear interpolation to get values of 
  ! shape parameters at origin
  gvec(:,0) = (gato_psi(2)*gvec(:,1)-&
       gato_psi(1)*gvec(:,2))/&
       (gato_psi(2)-gato_psi(1))

  ! Manage total flux variation 
  if (format_type == 0 .or. format_type == 7) then
     ! Case 1: raw gfile mode ; dpsi is undefined
     dpsi_data = gato_psi(nsurf)
     dpsi_gato = gato_psi(nsurf)
     do i=1,nx
        dpsi(i) = (i-1)*dpsi_gato/(nx-1)
     enddo
  else
     ! Case 2: typical case 
     dpsi_data = dpsi(nx)
     dpsi_gato = gato_psi(nsurf)

     ! Ensure max(dpsi) = max(gato_psi) 
     dpsi(:)  = dpsi(:)*dpsi_gato/dpsi_data
     ! Extra insurance against roundoff
     dpsi(nx) = dpsi_gato 
  endif

  ! Map shape coefficients onto poloidal flux (dpsi) grid:
  call cub_spline(gato_psi,gvec(1,:),nsurf+1,dpsi,rmin,nx)
  call cub_spline(gato_psi,gvec(2,:),nsurf+1,dpsi,zmag,nx)
  call cub_spline(gato_psi,gvec(3,:),nsurf+1,dpsi,rmaj,nx)
  call cub_spline(gato_psi,gvec(4,:),nsurf+1,dpsi,kappa,nx)
  call cub_spline(gato_psi,gvec(5,:),nsurf+1,dpsi,delta,nx)
  call cub_spline(gato_psi,gvec(6,:),nsurf+1,dpsi,zeta,nx)

  ! Explicitly set rmin=0 at origin
  rmin(1) = 0.0

  deallocate(gvec)
  !----------------------------------------------------

  !----------------------------------------------------
  ! Get general geometry coefficients
  !
  ! When verbose_flag=1, fluxfit will echo lots of data
  call fluxfit_driver(2,nfourier,nsurf,narc,gato_bigr,gato_bigz,verbose_flag)

  allocate(g3vec(4,0:nfourier,0:nsurf))
  allocate(g3rho(4,0:nfourier,nx))

  open(unit=1,file='fluxfit.geo',status='old')
  read(1,*) ip
  read(1,*) g3vec(:,:,1:nsurf)
  close(1)

  ! Use linear interpolation to get values of 
  ! (general) shape parameters at origin
  g3vec(:,:,0) = (gato_psi(2)*g3vec(:,:,1)-&
       gato_psi(1)*g3vec(:,:,2))/(gato_psi(2)-gato_psi(1))

  ! Map results onto poloidal flux (dpsi) grid:
  do i=1,4
     do ip=0,nfourier
        call cub_spline(&
             gato_psi,g3vec(i,ip,:),nsurf+1,dpsi,g3rho(i,ip,:),nx)
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

  !-------------------------------------------------------------
  ! Read q profile data in GATO's o1gta file and map to psi grid
  !
  allocate(gato_q(0:nsurf))
  open(unit=1,file='o1gta',status='old')

  do
     read(1,*) a,b,c
     if (a(1:1) == 'q') then
        if (c(1:2) == '(f') then
           read(1,*) gato_q(:)
           exit
        endif
     endif
  enddo

  close(1)

  call cub_spline(gato_psi,gato_q,nsurf+1,dpsi,q_gato,nx)
  if (nogatoq_flag == 0 .or. format_type == 3 .or. format_type == 7) then
     q(:) = q_gato(:)
  endif
  !-------------------------------------------------------------  

  ! Clean up
  deallocate(g3vec)
  deallocate(g3rho)
  deallocate(gato_psi)
  deallocate(gato_volume)
  deallocate(gato_dvoldpsi)
  deallocate(gato_bigr)
  deallocate(gato_bigz)
  deallocate(gato_q)

end subroutine prgen_read_gato
