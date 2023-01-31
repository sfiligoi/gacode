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
  complex, dimension(*), intent(inout) :: fxmany
  real, dimension(*), intent(inout) :: uvmany
  real, dimension(*), intent(inout) :: uxmany
  complex, dimension(*), intent(inout) :: fvmany
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
!$acc host_data use_device(uvmany,fnmany)
#ifdef HIPGPU
  rc = hipfftExecD2Z(plan_r2c_many,uvmany,fnmany)
#else
  rc = cufftExecD2Z(plan_r2c_many,uvmany,fnmany)
#endif
!$acc wait
!$acc end host_data
!$acc wait

end subroutine cgyro_do_fft


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

!$acc wait
  call cgyro_do_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany,comp_uxmany,comp_fvmany)

end program cgyro_test_fft
