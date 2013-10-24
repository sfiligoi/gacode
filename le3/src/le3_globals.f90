module le3_globals

  real, parameter :: pi=3.141592653589793

  logical :: initialized = .false.

  ! Resolution parameters
  integer :: nt,np
  integer :: nts,nps
  integer :: m,n

  ! Geometry parameters
  real :: rmin,rmaj,hmin,dhmindr,q
  real :: kappa, s_kappa, delta, s_delta, zeta, s_zeta
  real :: shift, zmag, dzmag
  real :: iota
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

end module le3_globals
