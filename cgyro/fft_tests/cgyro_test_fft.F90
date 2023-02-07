program cgyro_test_fft

  use, intrinsic :: iso_c_binding
  implicit none
#ifdef _OPENACC
  
#ifdef HIPGPU
  ! HIP
  type(C_PTR) :: plan_c2r_many
  type(C_PTR) :: plan_r2c_many
#else
  ! CUDA
  integer(c_int) :: plan_c2r_many
  integer(c_int) :: plan_r2c_many
#endif

#else
  ! FFTW
  type(C_PTR) :: plan_c2r_many
  type(C_PTR) :: plan_r2c_many
#endif
  ! nl03 n=256
  complex, dimension(0:95,0:767,72) :: fxmany
  real, dimension(0:189,0:767,72)  :: uvmany
  real, dimension(0:189,0:767,72) :: exp_uxmany
  complex, dimension(0:95,0:767,72) :: exp_fvmany
  real, dimension(0:189,0:767,72) :: comp_uxmany
  complex, dimension(0:95,0:767,72) :: comp_fvmany
  integer :: i,j
  integer :: start_count, end_count 
  integer :: count_rate, count_max

!$acc enter data create(fxmany,uvmany, exp_uxmany, exp_fvmany) &
!$acc&           create(comp_uxmany, comp_fvmany)

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
  call cgyro_setup_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany,comp_uxmany,comp_fvmany)

  open(unit=1,file="data/bin.fxmany.raw",form='unformatted',access='direct',recl=size(fxmany)*16)
  read(1,rec=1) fxmany
  close(1)
!$acc update device(fxmany) async
  open(unit=1,file="data/bin.uvmany.raw",form='unformatted',access='direct',recl=size(uvmany)*8)
  read(1,rec=1) uvmany
  close(1)
!$acc update device(uvmany) async

!$acc wait
  call cgyro_do_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany,comp_uxmany,comp_fvmany)

  call cgyro_comp_D(190*768*72,comp_uxmany,exp_uxmany);
  call cgyro_comp_Z(96*768*72,comp_fvmany,exp_fvmany);

  open(unit=1,file="data/bin.fxmany.raw",form='unformatted',access='direct',recl=size(fxmany)*16)
  read(1,rec=1) fxmany
  close(1)
!$acc update device(fxmany)

  do j=1,5
   call SYSTEM_CLOCK(start_count, count_rate, count_max)
   do i=1,100
    call cgyro_ident_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany)
   enddo
   call SYSTEM_CLOCK(end_count, count_rate, count_max)

   if ((end_count<start_count).and.(count_max>0)) end_count = end_count + count_max
   write(*,*) "2x100 FFT took ", (1.0*(end_count-start_count))/count_rate, " seconds"
  enddo

contains
  subroutine cgyro_comp_D(nels,comp_many,exp_many)
     implicit none
     !-----------------------------------
     integer, intent(in) :: nels
     real, dimension(*), intent(in) :: comp_many,exp_many
     !-----------------------------------
     integer :: n
     real :: comp_val,exp_val

!$acc parallel loop present(comp_many,exp_many) copyin(nels) private(comp_val,exp_val)
     do n=1,nels
       comp_val=comp_many(n)
       exp_val=exp_many(n)
       if (exp_val<1.e-5) then
          if (abs(comp_val - exp_val) > 1.e-10)  then
             write(*,*) "ERROR: cgyro_comp_D ",n,comp_val,exp_val
          endif
       else
          if (abs(comp_val - exp_val) > 1.e-7)  then
             write(*,*) "ERROR: cgyro_comp_D ",n,comp_val,exp_val
          endif
       endif
     enddo
     write(*,*) "cgyro_comp_D completed"

  end subroutine cgyro_comp_D

  subroutine cgyro_comp_Z(nels,comp_many,exp_many)
     implicit none
     !-----------------------------------
     integer, intent(in) :: nels
     complex, dimension(*), intent(in) :: comp_many,exp_many
     !-----------------------------------
     integer :: n
     real :: comp_val,exp_val

!$acc parallel loop present(comp_many,exp_many) copyin(nels) private(comp_val,exp_val)
     do n=1,nels
       comp_val=real(comp_many(n))
       exp_val=real(exp_many(n))
       if (exp_val<1.e-5) then
          if (abs(comp_val - exp_val) > 1.e-10)  then
             write(*,*) "ERROR: cgyro_comp_Z real ",n,comp_val,exp_val
          endif
       else
          if (abs(comp_val - exp_val) > 1.e-7)  then
             write(*,*) "ERROR: cgyro_comp_Z real ",n,comp_val,exp_val
          endif
       endif
       comp_val=aimag(comp_many(n))
       exp_val=aimag(exp_many(n))
       if (exp_val<1.e-5) then
          if (abs(comp_val - exp_val) > 1.e-10)  then
             write(*,*) "ERROR: cgyro_comp_Z imag ",n,comp_val,exp_val
          endif
       else
          if (abs(comp_val - exp_val) > 1.e-7)  then
             write(*,*) "ERROR: cgyro_comp_Z imag ",n,comp_val,exp_val
          endif
       endif
     enddo
     write(*,*) "cgyro_comp_Z completed"

  end subroutine cgyro_comp_Z

  subroutine cgyro_ident_fft(plan_c2r_many,plan_r2c_many,fxmany,uxmany)

     use, intrinsic :: iso_c_binding
#ifdef _OPENACC
#ifdef HIPGPU
     use hipfort_hipfft
#else
     use cufft
#endif
#endif

     implicit none
     !-----------------------------------
#ifndef _OPENACC
     include 'fftw3.f03'
#endif

#ifdef _OPENACC
  
#ifdef HIPGPU
     type(C_PTR), intent(inout) :: plan_c2r_many
     type(C_PTR), intent(inout) :: plan_r2c_many
#else
     integer(c_int), intent(inout) :: plan_c2r_many
     integer(c_int), intent(inout) :: plan_r2c_many
#endif

#else
     ! FFTW
     type(C_PTR), intent(inout) :: plan_c2r_many
     type(C_PTR), intent(inout) :: plan_r2c_many
#endif
     complex, dimension(:,:,:), intent(inout) :: fxmany
     real, dimension(:,:,:), intent(inout) :: uxmany
     !-----------------------------------
     integer :: rc

     ! --------------------------------------
     ! Forward
     ! --------------------------------------
!$acc wait
!$acc  host_data &
!$acc& use_device(fxmany,uxmany) 

#ifdef _OPENACC
  
#ifdef HIPGPU
     rc = hipfftExecZ2D(plan_c2r_many,c_loc(fxmany),c_loc(uxmany))
#else
     rc = cufftExecZ2D(plan_c2r_many,fxmany,uxmany)
#endif

#else
     ! FFTW
     call fftw_execute_dft_c2r(plan_c2r_many,fxmany,uxmany) 
     rc = 0
#endif
     if (rc/=0) then
        write(*,*) "ERROR: fftExec D2Z failed! ", rc
        call abort
     endif

!$acc wait
!$acc end host_data

     ! --------------------------------------
     ! Backward
     ! --------------------------------------

!$acc wait
!$acc host_data use_device(uxmany,fxmany)
#ifdef _OPENACC
  
#ifdef HIPGPU
     rc = hipfftExecD2Z(plan_r2c_many,c_loc(uxmany),c_loc(fxmany))
#else
     rc = cufftExecD2Z(plan_r2c_many,uxmany,fxmany)
#endif

#else
     ! FFTW
     call fftw_execute_dft_r2c(plan_r2c_many,uxmany,fxmany) 
     rc = 0
#endif

!$acc wait
!$acc end host_data
     if (rc/=0) then
        write(*,*) "ERROR: fftExec Z2D failed! ", rc
        call abort
     endif
!$acc wait

  end subroutine cgyro_ident_fft


  subroutine cgyro_do_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany,uxmany,fvmany)

     use, intrinsic :: iso_c_binding
#ifdef _OPENACC
#ifdef HIPGPU
     use hipfort_hipfft
#else
     use cufft
#endif
#endif

     implicit none
     !-----------------------------------
#ifndef _OPENACC
     include 'fftw3.f03'
#endif

#ifdef _OPENACC
  
#ifdef HIPGPU
     type(C_PTR), intent(inout) :: plan_c2r_many
     type(C_PTR), intent(inout) :: plan_r2c_many
#else
     integer(c_int), intent(inout) :: plan_c2r_many
     integer(c_int), intent(inout) :: plan_r2c_many
#endif

#else
     ! FFTW
     type(C_PTR), intent(inout) :: plan_c2r_many
     type(C_PTR), intent(inout) :: plan_r2c_many
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

#ifdef _OPENACC
  
#ifdef HIPGPU
     rc = hipfftExecZ2D(plan_c2r_many,c_loc(fxmany),c_loc(uxmany))
#else
     rc = cufftExecZ2D(plan_c2r_many,fxmany,uxmany)
#endif

#else
     ! FFTW
     call fftw_execute_dft_c2r(plan_c2r_many,fxmany,uxmany) 
     rc = 0
#endif
     if (rc/=0) then
        write(*,*) "ERROR: fftExec D2Z failed! ", rc
        call abort
     else
        write(*,*) "INFO: fftExec D2Z done."
        call flush(6)
     endif

!$acc wait
!$acc end host_data

     ! --------------------------------------
     ! Backward
     ! --------------------------------------

!$acc wait
!$acc host_data use_device(uvmany,fvmany)
#ifdef _OPENACC
  
#ifdef HIPGPU
     rc = hipfftExecD2Z(plan_r2c_many,c_loc(uvmany),c_loc(fvmany))
#else
     rc = cufftExecD2Z(plan_r2c_many,uvmany,fvmany)
#endif

#else
     ! FFTW
     call fftw_execute_dft_r2c(plan_r2c_many,uvmany,fvmany) 
     rc = 0
#endif

!$acc wait
!$acc end host_data
     if (rc/=0) then
        write(*,*) "ERROR: fftExec Z2D failed! ", rc
        call abort
     else
        write(*,*) "INFO: fftExec Z2D done."
        call flush(6)
     endif
!$acc wait

  end subroutine cgyro_do_fft


  subroutine cgyro_setup_fft(plan_c2r_many,plan_r2c_many,fxmany,uvmany,uxmany,fvmany)

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
#ifdef _OPENACC
#ifdef HIPGPU
  use hipfort_hipfft
#else
  use cufft
#endif
#endif
  implicit none
#ifndef _OPENACC
     include 'fftw3.f03'
     integer, external :: omp_get_max_threads
#endif
  !-----------------------------------
#ifdef _OPENACC
  
#ifdef HIPGPU
  type(C_PTR), intent(inout) :: plan_c2r_many
  type(C_PTR), intent(inout) :: plan_r2c_many
#else
  integer(c_int), intent(inout) :: plan_c2r_many
  integer(c_int), intent(inout) :: plan_r2c_many
#endif

#else
  ! FFTW
  type(C_PTR), intent(inout) :: plan_c2r_many
  type(C_PTR), intent(inout) :: plan_r2c_many
#endif
  complex, dimension(:,:,:), intent(inout) :: fxmany
  real, dimension(:,:,:), intent(inout) :: uvmany
  real, dimension(:,:,:), intent(inout) :: uxmany
  complex, dimension(:,:,:), intent(inout) :: fvmany
  !-----------------------------------
  integer :: istatus
  integer, parameter :: irank = 2
  integer, dimension(irank) :: ndim,inembed,onembed
  integer :: idist,odist,istride,ostride,nsplit

#ifndef _OPENACC
     istatus = fftw_init_threads()
     call fftw_plan_with_nthreads(omp_get_max_threads())
#endif

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

#ifdef _OPENACC
  
#ifdef HIPGPU
     plan_c2r_many = c_null_ptr
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

#else
  ! FFTW
     plan_c2r_many = fftw_plan_many_dft_c2r(&
          irank, &
          ndim, &
          nsplit, &
          fxmany, &
          inembed, &
          istride, &
          idist, &
          uxmany, &
          onembed, &
          ostride, &
          odist, &
          FFTW_PATIENT)
     istatus = 0
#endif
     if (istatus/=0) then
        write(*,*) "ERROR: fftPlanMany Z2D failed! ", istatus
        call abort
     else
        write(*,*) "INFO: fftPlanMany Z2D done."
        call flush(6)
     endif

     idist = 190*768
     odist = 96*768
     inembed = 190
     onembed = 96
     istride = 1
     ostride = 1
#ifdef _OPENACC
  
#ifdef HIPGPU
     plan_r2c_many = c_null_ptr
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

#else
  ! FFTW
     plan_r2c_many = fftw_plan_many_dft_r2c(&
          irank, &
          ndim, &
          nsplit, &
          uvmany, &
          inembed, &
          istride, &
          idist, &
          fvmany, &
          onembed, &
          ostride, &
          odist, &
          FFTW_PATIENT)
     istatus = 0
#endif
     if (istatus/=0) then
        write(*,*) "ERROR: fftPlanMany D2Z failed! ", istatus
        call abort
     else
        write(*,*) "INFO: fftPlanMany D2Z done."
        call flush(6)
     endif

  end subroutine cgyro_setup_fft

end program cgyro_test_fft
