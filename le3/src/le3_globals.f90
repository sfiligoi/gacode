module le3_globals

  real, parameter :: pi=3.141592653589793

  integer :: nt,np
  integer :: m,n
  real :: rmin,rmaj,hmin,q
  real :: dt,dp
  real :: tol
  integer :: restart_flag

  real, dimension(:), allocatable    :: t,p
  real, dimension(:,:), allocatable  :: tb
  real, dimension(:,:), allocatable  :: dtbdp,dtbdt
  real, dimension(:,:), allocatable  :: r,z
  real, dimension(:,:), allocatable  :: drdtb,dzdtb
  real, dimension(:,:), allocatable  :: drdpb,dzdpb
  real, dimension(:,:), allocatable  :: drdr,dzdr
  real, dimension(:,:), allocatable  :: jac
  real, dimension(:,:), allocatable  :: bp,br,bz
  real, dimension(:,:), allocatable  :: rp,rt
  real, dimension(:,:), allocatable  :: zp,zt
  real, dimension(:,:), allocatable  :: fp,ft
  real, dimension(:,:), allocatable  :: fpt,ftp

  integer, dimension(:), allocatable :: tcyc, pcyc
  real, dimension(-2:2) :: cderiv

end module le3_globals
