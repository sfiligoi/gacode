module le3_globals

  real, parameter :: pi=3.141592653589793

  logical :: initialized = .false.

  ! Resolution parameters
  integer :: nt,np
  integer :: nts,nps
  integer :: m,n

  ! Geometry parameters
  real :: rmin,rmaj,hmin,dhmindr,q,s
  real :: kappa, s_kappa, delta, s_delta, zeta, s_zeta
  real :: shift, zmag, dzmag
  real :: iota,iota_p
  real :: dt,dp

  real, dimension(:), allocatable    :: t,p
  real, dimension(:,:), allocatable  :: tb
  real, dimension(:,:), allocatable  :: dtbdp,dtbdt
  real, dimension(:,:), allocatable  :: r,z
  real, dimension(:,:), allocatable  :: drdtb,dzdtb
  real, dimension(:,:), allocatable  :: drdpb,dzdpb
  real, dimension(:,:), allocatable  :: jac
  real, dimension(:,:), allocatable  :: rp,rt
  real, dimension(:,:), allocatable  :: zp,zt
  real, dimension(:,:), allocatable  :: fp,ft
  real, dimension(:,:), allocatable  :: as,bs,cs,ds
  real, dimension(:,:), allocatable  :: sinm,cosm
  real, dimension(:,:), allocatable  :: sinn,cosn

  real, dimension(:,:), allocatable :: gpp,gtt,gpt
  real, dimension(:,:), allocatable :: cosu,rc
  real, dimension(:,:), allocatable :: rs,zs
  real, dimension(:,:), allocatable :: bpol, btor
  real, dimension(:,:), allocatable :: dbdt, dbdp

  !--------------------------------------------------------
  ! MINPACK variables
  integer :: info
  integer :: nwork
  integer :: msize
  real, dimension(:), allocatable :: work
  !
  real :: tol
  real, dimension(:), allocatable :: xfunc
  real, dimension(:), allocatable :: yfunc
  !--------------------------------------------------------

  real, dimension(:,:), allocatable :: bmag
  real, dimension(:,:), allocatable :: bdotgrad
  real, dimension(:,:), allocatable :: bdotgradB_overB
  real, dimension(:,:), allocatable :: vdrift_x,dgdp,vexb_dt,vexb_dp
  real, dimension(:,:), allocatable :: g
  real, dimension(:,:), allocatable :: mat_stream_dt
  real, dimension(:,:), allocatable :: mat_stream_dp
  real, dimension(:,:), allocatable :: mat_trap
  real, dimension(:,:), allocatable :: mat_coll
  real, dimension(:,:), allocatable :: mat_vexb_dt, mat_vexb_dp
  real :: vprime
  integer :: matsize
  integer :: indx_c00
  real, dimension(:,:), allocatable :: basis, basis_prime, &
       basis_dt_prime, basis_dp_prime
  integer, dimension(:), allocatable :: m_indx, n_indx, itype

end module le3_globals
