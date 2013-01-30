!------------------------------------------------------------
! expromake.f90
!
! PURPOSE:
!------------------------------------------------------------

program expromake

  use expromake_globals
  use EXPRO_interface

  implicit none

  integer :: ierr


  call expromake_read_input 

  EXPRO_ctrl_density_method = 1
  EXPRO_ctrl_z(1:3) = exm_z(1:3)
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_signq = (-1)*(-1)
  EXPRO_ctrl_signb = -(-1)
  EXPRO_ctrl_rotation_method = 1

  ! We're going to see if this file exists
  open(unit=1,file='input.profiles.gen',status='old',iostat=ierr)

  if (ierr > 0) then

     mode = 1

     ! Starting from scratch (no input.profiles)

     ! Setting EXPRO_n_exp > 0 skips reading input.profiles

     nx = 11
 
     EXPRO_n_exp  = nx
     EXPRO_ncol   = 5
     EXPRO_nblock = 8

     call EXPRO_alloc('./',1) 

     ! Set initial values of parameters
     call expromake_init

  else

     mode = 2

     ! Starting from pre-existing input.profiles

     close(1)
     call EXPRO_alloc('./',1) 
     nx = EXPRO_n_exp
     call EXPRO_read

  endif

  call expromake_write
  call EXPRO_alloc('./',0) 

end program expromake
