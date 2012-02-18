!\************************************************************
! tglf_init.f90: initialize the parallel TGLF model, it
!     will assign the number of ky spectrums between CPUs,
!     and setup the group communicators, total number of
!     group PEs, and local group PE ids, et al.
!
! REVERSION HISTORY:
! 26 May 2011, xyuan@pppl.gov, Version 1.0
!/************************************************************
  subroutine tglf_mpi_init
    use tglf_mpi
    use tglf_kyspectrum
    implicit none

    ! local variables
    integer, allocatable, dimension(:) :: nkyProc
    integer :: npts,nres
    integer :: i, ierr


    !\--------------------------
    ! exective code starts here
    !/

    ! free memory first
    if(allocated(nkyProc)) deallocate(nkyProc)
    if(allocated(nky0))    deallocate(nky0)
    if(allocated(nky1))    deallocate(nky1)
 
    ! allocate memory
    allocate(nkyProc(nProcTglf),nky0(nProcTglf),nky1(nProcTglf),stat=ierr)
    if(ierr /= 0) then
       write(*,*) ' memory allocation error!'
       return
    end if

    CALL get_ky_spectrum

    ! initialize nkyProc
    nkyProc(:) = 0
    nky0(:)    = 0
    nky1(:)    = 0

    ! setup the number of spectrum for each CPUs
    if(nProcTglf > nky) then
      ! warning message
      write(*,*) ' too many PEs applied'
      nkyProc(1:nky) = 1
      
      do i = 1, nky
        nky0(i) = i
        nky1(i) = i
      end do

    else
      npts = nky/nProcTglf
      nres = mod(nky,nProcTglf)
      
      nkyProc(:) = npts
      if(nres /= 0) then
        do i = 1, nres
          nkyProc(i) = npts + 1
        end do
      end if

      nky0(1) = 1
      nky1(1) = nkyProc(1)

      do i = 2, nProcTglf
        nky0(i) = nky1(i-1)+1
        nky1(i) = nky0(i)+nkyProc(i)-1 
      end do

    end if
    
  end subroutine tglf_mpi_init
  !\**********************************
  ! this is the END OF FILE
  !/**********************************
