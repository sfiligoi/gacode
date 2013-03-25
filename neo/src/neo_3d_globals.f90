module neo_3d_globals

  real, dimension(:), allocatable     :: varphi
  integer, dimension(:), allocatable  :: ip_indx
  real, dimension(:,:), allocatable   :: driftx_2d
  real, dimension(:,:,:), allocatable :: driftx_3d
  integer, dimension(:), allocatable  :: vpcyc

end module neo_3d_globals
