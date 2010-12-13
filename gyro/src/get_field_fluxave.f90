!------------------------------------------------------
! get_field_fluxave.f90 [caller gyro_write_master]
!
! PURPOSE:
!  Compute fast flux-surface average of (phi,a)
!   
! NOTES:
!  Because c_fluxave is zero for n > 0, the 
!  quantities phi_fluxave and a_fluxave will 
!  also be zero for n > 0.
!------------------------------------------------------

subroutine get_field_fluxave

  use gyro_globals
  use gyro_pointers

  !------------------------------------
  implicit none
  !------------------------------------

  phi_fluxave(:) = 0.0
  a_fluxave(:)   = 0.0
  aperp_fluxave(:)   = 0.0

  do i=1,n_x
     do j=1,n_blend

        phi_fluxave(i) = phi_fluxave(i)+c_fluxave(j,i)*field_blend(j,i,1)

     enddo ! j
  enddo ! i

  if (n_field > 1) then
     do i=1,n_x
        do j=1,n_blend

           a_fluxave(i)   = a_fluxave(i)+c_fluxave(j,i)*field_blend(j,i,2)

        enddo ! j
     enddo ! i
  endif

  if (n_field > 2) then
     do i=1,n_x
        do j=1,n_blend

           aperp_fluxave(i)   = aperp_fluxave(i)+c_fluxave(j,i)*field_blend(j,i,3)

        enddo ! j
     enddo ! i
  endif

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[get_field_fluxave called]'
  endif

end subroutine get_field_fluxave
