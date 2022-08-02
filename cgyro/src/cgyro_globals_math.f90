!-----------------------------------------------------------------
! cgyro_globals_math.f90
!
! PURPOSE:
!  Helper functions for scalable math operations
!  on matrices with known sizes
!-----------------------------------------------------------------

module cgyro_globals_math

    implicit none 

contains


  !=========================================================
  ! Velocity distributed arrays
  !=========================================================

subroutine cgyro_vel_inplace_fma1(left, c1, r1)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(inout), dimension(nc,nv_loc) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(nc,nv_loc) :: r1
    !-------------------------------------------------------
    call cgyro_cmpl_inplace_fma1(nc*nv_loc, left,c1,r1)
end subroutine cgyro_vel_inplace_fma1

subroutine cgyro_vel_inplace_fma2(left, c1, r1, c2, r2)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(inout), dimension(nc,nv_loc) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(nc,nv_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc) :: r2
    !-------------------------------------------------------
    call cgyro_cmpl_inplace_fma2(nc*nv_loc, left,c1,r1,c2,r2)
end subroutine cgyro_vel_inplace_fma2


subroutine cgyro_vel_inplace_fma3(left, c1, r1, c2, r2, c3, r3)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(inout), dimension(nc,nv_loc) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(nc,nv_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc) :: r3
    !-------------------------------------------------------
    call cgyro_cmpl_inplace_fma3(nc*nv_loc, left,c1,r1,c2,r2,c3,r3)
end subroutine cgyro_vel_inplace_fma3

subroutine cgyro_vel_inplace_fma4(left, c1, r1, c2, r2, c3, r3, c4, r4)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(inout), dimension(nc,nv_loc) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(nc,nv_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(nc,nv_loc) :: r4
    !-------------------------------------------------------
    call cgyro_cmpl_inplace_fma4(nc*nv_loc, left,c1,r1,c2,r2,c3,r3,c4,r4)
end subroutine cgyro_vel_inplace_fma4

subroutine cgyro_vel_copy(left, r1)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc) :: left
    complex, intent(in), dimension(nc,nv_loc) :: r1
    !-------------------------------------------------------
    call cgyro_cmpl_copy(nc*nv_loc, left,r1)
end subroutine cgyro_vel_copy

subroutine cgyro_vel_fma2(left, r1, c2, r2)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc) :: left
    complex, intent(in), dimension(nc,nv_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc) :: r2
    !-------------------------------------------------------
    call cgyro_cmpl_fma2(nc*nv_loc, left,r1,c2,r2)
end subroutine cgyro_vel_fma2

subroutine cgyro_vel_fma3(left, r1, c2, r2, c3, r3)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc) :: left
    complex, intent(in), dimension(nc,nv_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc) :: r3
    !-------------------------------------------------------
    call cgyro_cmpl_fma3(nc*nv_loc, left,r1,c2,r2,c3,r3)
end subroutine cgyro_vel_fma3

subroutine cgyro_vel_fma4(left, r1, c2, r2, c3, r3, c4, r4)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc) :: left
    complex, intent(in), dimension(nc,nv_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(nc,nv_loc) :: r4
    !-------------------------------------------------------
    call cgyro_cmpl_fma4(nc*nv_loc, left,r1,c2,r2,c3,r3,c4,r4)
end subroutine cgyro_vel_fma4

subroutine cgyro_vel_fma5(left, r1, c2, r2, c3, r3, c4, r4, c5, r5)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc) :: left
    complex, intent(in), dimension(nc,nv_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(nc,nv_loc) :: r4
    real, intent(in) :: c5
    complex, intent(in), dimension(nc,nv_loc) :: r5
    !-------------------------------------------------------
    call cgyro_cmpl_fma5(nc*nv_loc, left,r1,c2,r2,c3,r3,c4,r4,c5,r5)
end subroutine cgyro_vel_fma5

end module cgyro_globals_math


