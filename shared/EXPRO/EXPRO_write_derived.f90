!--------------------------------------------------------
! EXPRO_write_derived.f90
!
! PURPOSE:
!  Write all derived quantities to input.profiles.extra
!--------------------------------------------------------

subroutine EXPRO_write_derived

  use EXPRO_globals
  use EXPRO_interface

  integer, parameter :: io=1

  open(unit=io,file=trim(path)//'input.profiles.extra',status='replace')

  write(io,10) EXPRO_bunit(:)    ! 1
  write(io,10) EXPRO_s(:)        ! 2
  write(io,10) EXPRO_drmaj(:)    ! 3
  write(io,10) EXPRO_dzmag(:)    ! 4 
  write(io,10) EXPRO_sdelta(:)   ! 5  
  write(io,10) EXPRO_skappa(:)   ! 6
  write(io,10) EXPRO_szeta(:)    ! 7
  write(io,10) EXPRO_dlnnedr(:)  ! 8
  write(io,10) EXPRO_dlntedr(:)  ! 9 
  write(io,10) transpose(EXPRO_dlnnidr(:,:)) ! 10-14 
  write(io,10) transpose(EXPRO_dlntidr(:,:)) ! 15-19
  write(io,10) EXPRO_dlnptotdr(:)! 20 
  write(io,10) EXPRO_drdrho(:)   ! 21
  write(io,10) EXPRO_w0p(:)      ! 22
  write(io,10) EXPRO_vol(:)      ! 23
  write(io,10) EXPRO_volp(:)     ! 24
  write(io,10) EXPRO_cs(:)       ! 25
  write(io,10) EXPRO_rhos(:)     ! 26
  write(io,10) EXPRO_ni_new(:)   ! 27
  write(io,10) EXPRO_dlnnidr_new(:) ! 28 
  write(io,10) EXPRO_grad_r0(:)     ! 29
  write(io,10) EXPRO_ave_grad_r(:)  ! 30 
  write(io,10) EXPRO_bp0(:)      ! 31
  write(io,10) EXPRO_bt0(:)      ! 32
  write(io,10) EXPRO_gamma_e(:)  ! 33
  write(io,10) EXPRO_gamma_p(:)  ! 34
  write(io,10) EXPRO_mach(:)     ! 35

  close(io)

10 format(t2,1pe14.7)

end subroutine EXPRO_write_derived
