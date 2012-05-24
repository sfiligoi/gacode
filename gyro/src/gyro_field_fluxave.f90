!------------------------------------------------------
! gyro_field_fluxave.f90
!
! PURPOSE:
!  Compute fast flux-surface average of fields.
!   
! NOTES:
!  Because c_fluxave is zero for n > 0, the 
!  value of field_fluxave will also be zero 
!  for n > 0.
!------------------------------------------------------

subroutine gyro_field_fluxave

  use gyro_globals
  use gyro_pointers

  !------------------------------------
  implicit none
  !------------------------------------

  field_fluxave(:,:) = 0.0

  do ix=1,n_field
     do i=1,n_x
        do j=1,n_blend
           field_fluxave(i,ix) = field_fluxave(i,ix) &
                + c_fluxave(j,i)*field_blend(j,i,ix)
        enddo ! j
     enddo ! i
  enddo ! ix

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_field_fluxave called]'
  endif

end subroutine gyro_field_fluxave
