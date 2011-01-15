module gyro_fieldeigen_private

  integer, parameter :: iwmax=20
  !
  integer :: n_eigen
  integer :: info_eigen
  integer, dimension(:), allocatable :: i_piv_eigen
  complex, dimension(:,:), allocatable :: a_eigen
  complex, dimension(:,:), allocatable :: a_eigen_loc
  complex, dimension(:), allocatable :: b_eigen
  !
  complex, dimension(:,:,:,:,:), allocatable :: vdotgrad
  complex, dimension(:,:),   allocatable :: v_omegastar
  complex, dimension(:,:,:,:,:), allocatable :: fg
  !
  integer :: n_im
  integer :: info_im
  integer, dimension(:), allocatable :: i_piv_im
  complex, dimension(:,:), allocatable :: propinv
  complex, dimension(:,:), allocatable :: prod
  complex, dimension(:), allocatable :: work_im
  complex, dimension(:,:), allocatable :: gk_left
  complex, dimension(:,:), allocatable :: gk_right
  !
  integer :: n_imk
  integer :: info_imk
  integer, dimension(:), allocatable :: i_piv_imk
  complex, dimension(:,:), allocatable :: propinvk
  complex, dimension(:,:), allocatable :: prodk
  complex, dimension(:), allocatable :: work_imk
  complex, dimension(:,:), allocatable :: gk_leftk
  complex, dimension(:,:), allocatable :: gk_rightk
  !
  real :: sgn
  real, dimension(:), allocatable :: diag_scale
  complex :: det
  complex, dimension(:,:), allocatable :: cmk
  !
  real :: error_eigen
  complex :: omega_eigen
  !
  ! Stacked indices
  integer :: im, impr
  integer :: mp, mk, mkp
  integer :: ie2
  integer :: ipp, ippp, ij, ijp
  integer :: imk, imkp, kp
  integer :: idiff_ipp,idiff_ippp,ixp
  !
  ! Constants for BLAS
  !
  complex :: cmplx_0=(0.0,0.0)
  complex :: cmplx_1=(1.0,0.0)
  !---------------------------------------------------

end module gyro_fieldeigen_private
