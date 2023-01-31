program cgyro_test_fft

  use, intrinsic :: iso_c_binding
  implicit none
#ifdef HIPGPU
  type(C_PTR) :: plan_c2r_many
  type(C_PTR) :: plan_r2c_many
#else
  integer(c_int) :: plan_c2r_many
  integer(c_int) :: plan_r2c_many
#endif
  ! nl03 n=256
  complex, dimension(0:95,0:767,72) :: fxmany
  real, dimension(0:189,0:767,72)  :: uvmany
  real, dimension(0:189,0:767,72) :: exp_uxmany
  complex, dimension(0:95,0:767,72) :: exp_fvmany
  real, dimension(0:189,0:767,72) :: comp_uxmany
  complex, dimension(0:95,0:767,72) :: comp_fvmany

!$acc enter data create(fxmany,uvmany, exp_uxmany, exp_fvmany) &
!$acc&           create(comp_uxmany)

  open(unit=1,file="data/bin.fxmany.raw",form='unformatted',access='direct',recl=size(fxmany)*16)
  read(1,rec=1) fxmany
  close(1)
!$acc update device(fxmany) async
  open(unit=1,file="data/bin.uvmany.raw",form='unformatted',access='direct',recl=size(uvmany)*8)
  read(1,rec=1) uvmany
  close(1)
!$acc update device(uvmany) async
  open(unit=1,file="data/bin.uxmany.raw",form='unformatted',access='direct',recl=size(exp_uxmany)*8)
  read(1,rec=1) exp_uxmany
  close(1)
!$acc update device(exp_uxmany) async
  open(unit=1,file="data/bin.final.raw",form='unformatted',access='direct',recl=size(exp_fvmany)*16)
  read(1,rec=1) exp_fvmany
  close(1)
!$acc update device(exp_fvmany) async

  call cgyro_setup_fft(plan_c2r_many,plan_r2c_many)
!$acc wait
  call cgyro_do_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany,comp_uxmany,comp_fvmany)

contains
  subroutine cgyro_do_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany,uxmany,fvmany)

  use, intrinsic :: iso_c_binding
#ifdef HIPGPU
  use hipfort_hipfft
#else
  use cufft
#endif
  implicit none
  !-----------------------------------
#ifdef HIPGPU
  type(C_PTR), intent(inout) :: plan_c2r_many
  type(C_PTR), intent(inout) :: plan_r2c_many
#else
  integer(c_int), intent(inout) :: plan_c2r_many
  integer(c_int), intent(inout) :: plan_r2c_many
#endif
  complex, dimension(:,:,:), intent(inout) :: fxmany
  real, dimension(:,:,:), intent(inout) :: uvmany
  real, dimension(:,:,:), intent(inout) :: uxmany
     complex, dimension(:,:,:), intent(inout) :: fvmany
     !-----------------------------------
     integer :: rc

     ! --------------------------------------
     ! Forward
     ! --------------------------------------
!$acc wait
!$acc  host_data &
!$acc& use_device(fxmany,uxmany) 

#ifdef HIPGPU
  rc = hipfftExecZ2D(plan_c2r_many,fxmany,uxmany)
#else
  rc = cufftExecZ2D(plan_c2r_many,fxmany,uxmany)
#endif

!$acc wait
!$acc end host_data

     ! --------------------------------------
     ! Backward
     ! --------------------------------------

!$acc wait
!$acc host_data use_device(uvmany,fvmany)
#ifdef HIPGPU
  rc = hipfftExecD2Z(plan_r2c_many,uvmany,fvmany)
#else
  rc = cufftExecD2Z(plan_r2c_many,uvmany,fvmany)
#endif
!$acc wait
!$acc end host_data
!$acc wait

  end subroutine cgyro_do_fft


  subroutine cgyro_setup_fft(plan_c2r_many,plan_r2c_many)

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
#ifdef HIPGPU
  use hipfort_hipfft
#else
  use cufft
#endif
  implicit none
  !-----------------------------------
#ifdef HIPGPU
  type(C_PTR), intent(inout) :: plan_c2r_many
  type(C_PTR), intent(inout) :: plan_r2c_many
#else
  integer(c_int), intent(inout) :: plan_c2r_many
  integer(c_int), intent(inout) :: plan_r2c_many
#endif
  integer :: istatus
  integer, parameter :: irank = 2
  integer, dimension(irank) :: ndim,inembed,onembed
  integer :: idist,odist,istride,ostride,nsplit

     nsplit = 72

     ! nl03
     ndim(1) = 768
     ndim(2) = 190
     idist = 96*768
     odist = 190*768
     istride = 1
     ostride = 1
     inembed = 96
     onembed = 190

#ifdef HIPGPU
     istatus = hipfftCreate(plan_c2r_many)
     istatus = hipfftPlanMany(&
          plan_c2r_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          HIPFFT_Z2D, &
          nsplit)
#else
     istatus = cufftPlanMany(&
          plan_c2r_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          CUFFT_Z2D, &
          nsplit)
#endif
     idist = 190*768
     odist = 96*768
     inembed = 190
     onembed = 96
     istride = 1
     ostride = 1
#ifdef HIPGPU
     istatus = hipfftCreate(plan_r2c_many)
     istatus = hipfftPlanMany(&
          plan_r2c_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          HIPFFT_D2Z, &
          nsplit)
#else
     istatus = cufftPlanMany(&
          plan_r2c_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          CUFFT_D2Z, &
          nsplit)
#endif
  end subroutine cgyro_setup_fft




end program cgyro_test_fft
