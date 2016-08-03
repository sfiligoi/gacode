!------------------------------------------------------------
! tgyro_quasigrad.f90
!
! PURPOSE:
!  Simple routine to enforce derivative of
!
!             ne = sum_i n_i z_i
!
! according to different rules.
!------------------------------------------------------------

subroutine tgyro_quasigrad

  use tgyro_globals

  implicit none
  integer :: is
  integer :: n_off
  real, dimension(n_r) :: delta_n

  ! Compute density offset that needs to be corrected
  delta_n = 0.0
  if (evo_e(is) > 0) then
     delta_n(:) = delta_n(:) - ne(:)*dlnnedr(:)
  endif
  do is=0,loc_n_ion
     if (evo_e(is) > 0) then
        delta_n(:) = delta_n(:) + zi_vec(is)*ni(is,:)*dlnnidr(is,:)
     endif
  enddo

  ! Number of species to offset
  do is=0,loc_n_ion
     if (evo_e(is) == -1) then
        n_off = n_off+1
     endif
  enddo

  ! Perform offsets
  if (evo_e(0) == -1) then
     dlnnedr(:) = delta_n(:)/(ne(:)*n_off)
  endif
  do is=1,loc_n_ion
     if (evo_e(is) == -1) then
        dlnnidr(is,:) = -delta_n(:)/(zi_vec(is)*ni(is,:)*n_off)
     endif
  enddo

end subroutine tgyro_quasigrad

