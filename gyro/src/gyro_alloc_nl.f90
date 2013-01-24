!----------------------------------------------------------
! gyro_alloc_nl
!
! PURPOSE:
!  Create and destroy nonlinear work arrays.
!
! NOTES:
!  flag=0: deallocate
!  flag=1: allocate
!----------------------------------------------------------

subroutine gyro_alloc_nl(flag)

  use gyro_globals 
  use gyro_nl_private

  integer, intent(in) :: flag

  if (flag == 1 .and. allocated(n_p)) then
     if (i_proc == 0) then
        print *,'WARNING: already allocated arrays in gyro_alloc_nl'
     endif
     return
  endif
  if (flag == 0 .and. .not.allocated(n_p)) then
     if (i_proc == 0) then
        print *,'WARNING: cannot deallocate arrays in gyro_alloc_nl'
     endif
     return
  endif

  if (flag == 1) then

     if (nonlinear_flag == 1) then
        ! Shared allocations
        allocate(n_p(-n_max:n_max))
        allocate(i_p(-n_max:n_max))
        allocate(c_nl_i(n_x))
     endif

  else 

     ! Shared deallocations
     if (allocated(n_p)) deallocate(n_p)
     if (allocated(i_p)) deallocate(i_p)
     if (allocated(c_nl_i)) deallocate(c_nl_i)

     ! Special FFTW2 deallocations
     if (allocated(v_fft)) deallocate(v_fft)
     if (allocated(vt_fft)) deallocate(vt_fft)

     ! Special FFTW3 deallocations
     if (allocated(v_fft3)) deallocate(v_fft3)
     if (allocated(vt_fft3)) deallocate(vt_fft3)

     ! Special ESSL deallocations
     if (allocated(aux1_dcrft)) deallocate(aux1_dcrft)
     if (allocated(aux2_dcrft)) deallocate(aux2_dcrft)
     if (allocated(aux1_drcft)) deallocate(aux1_drcft)
     if (allocated(aux2_drcft)) deallocate(aux2_drcft)

  endif

end subroutine gyro_alloc_nl
