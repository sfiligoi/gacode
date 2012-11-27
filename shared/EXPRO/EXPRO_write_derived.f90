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
  write(io,*) EXPRO_dlnptotdr(:)! 20 
  write(io,*) EXPRO_drdrho(:)   ! 21
  write(io,*) EXPRO_w0(:)       ! 22
  write(io,*) EXPRO_w0p(:)      ! 23
  write(io,*) EXPRO_vol(:)      ! 24
  write(io,*) EXPRO_volp(:)     ! 25
  write(io,*) EXPRO_cs(:)       ! 26
  write(io,*) EXPRO_rhos(:)     ! 27
  write(io,*) EXPRO_ni_new(:)   ! 28
  write(io,*) EXPRO_dlnnidr_new(:) ! 29 
  write(io,*) EXPRO_grad_r0(:)     ! 30
  write(io,*) EXPRO_ave_grad_r(:)  ! 31 
  write(io,*) EXPRO_bp0(:)      ! 32
  write(io,*) EXPRO_bt0(:)      ! 33
  write(io,*) EXPRO_poloidalfluxover2pi(:) ! 34
  write(io,*) EXPRO_gamma_e(:)  ! 35
  write(io,*) EXPRO_gamma_p(:)  ! 36
  write(io,*) EXPRO_mach(:)     ! 37

  close(io)

end subroutine EXPRO_write_derived
