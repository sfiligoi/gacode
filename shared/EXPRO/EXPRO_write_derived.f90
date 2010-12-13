!--------------------------------------------------------
! EXPRO_write_derived.f90
!
! PURPOSE:
!  Write all derived quantities to input.profiles.extra
!--------------------------------------------------------

subroutine EXPRO_write_derived(path)

  use EXPRO_interface

  character(len=*) :: path
  integer, parameter :: io=1

  open(unit=io,file=trim(path)//'input.profiles.extra',status='replace')

  write(io,*) EXPRO_bunit(:)    ! 1
  write(io,*) EXPRO_s(:)        ! 2
  write(io,*) EXPRO_drmaj(:)    ! 3
  write(io,*) EXPRO_dzmag(:)    ! 4 
  write(io,*) EXPRO_sdelta(:)   ! 5  
  write(io,*) EXPRO_skappa(:)   ! 6
  write(io,*) EXPRO_szeta(:)    ! 7
  write(io,*) EXPRO_dlnnedr(:)  ! 8
  write(io,*) EXPRO_dlntedr(:)  ! 9 
  write(io,*) transpose(EXPRO_dlnnidr(:,:)) ! 10-14 
  write(io,*) transpose(EXPRO_dlntidr(:,:)) ! 15-19
  write(io,*) EXPRO_w0(:)       ! 20
  write(io,*) EXPRO_w0p(:)      ! 21
  write(io,*) EXPRO_vol(:)      ! 22
  write(io,*) EXPRO_volp(:)     ! 23
  write(io,*) EXPRO_cs(:)       ! 24
  write(io,*) EXPRO_rhos(:)     ! 25

  close(io)

end subroutine EXPRO_write_derived
