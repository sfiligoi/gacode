!------------------------------------------------------------
! expromake.f90
!
! PURPOSE:
!------------------------------------------------------------

program expromake

  use expromake_globals
  use EXPRO_interface

  implicit none

  integer :: i
  integer :: ierr
  integer :: indx
  integer :: n_indx
  integer :: nx


  call expromake_read_input 

  EXPRO_ctrl_density_method = 1
  EXPRO_ctrl_z(1:3) = z(1:3)
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_signq = (-1)*(-1)
  EXPRO_ctrl_signb = -(-1)
  EXPRO_ctrl_rotation_method = 1

  open(unit=1,file='input.profiles.gen',status='old',iostat=ierr)

  if (ierr > 0) then
     ! Setting EXPRO_n_exp > 0 skips reading input.profiles
     EXPRO_n_exp  = 10
     EXPRO_ncol   = 5
     EXPRO_nblock = 8
     call EXPRO_alloc('./',1) 
     call expromake_init
  else
     close(1)
     ! Have input.profiles
     call EXPRO_alloc('./',1) 
     call EXPRO_read
  endif

  nx = EXPRO_n_exp

  open(unit=1,file='input.profiles',status='replace')

  !---------------------------------------------------------------
  ! Basic information
  !
  write(1,20) '# expromake '
  write(1,30) '#    RADIAL GRIDPOINTS : ',nx
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Merged files and extra information
  !
  write(1,20) '# * Additional information:'
  write(1,20) '#'
  write(1,'(a,1pe8.2,a)') '#  PLASMA MAJOR RADIUS : ',EXPRO_rmaj(nx),' m'
  write(1,'(a,1pe8.2,a)') '#  PLASMA MINOR RADIUS : ',EXPRO_rmin(nx),' m'
  write(1,20) '#'
  !---------------------------------------------------------------

  write(1,20) '# '
  if (nx > 99) then
     write(1,30) 'N_EXP=',nx
  else
     write(1,25) 'N_EXP=',nx
  endif

  write(1,'(a,sp,1pe14.7)') 'BT_EXP=',EXPRO_b_ref
  write(1,60) 'ARHO_EXP=',EXPRO_arho

  n_indx = size(EXPRO_tag)

  do indx=1,n_indx,5

     write(1,20) '# '
     write(1,20) '#',EXPRO_tag(indx:indx+4)
     select case (indx)
     case(1)
        do i=1,nx
           write(1,10) EXPRO_rho(i),EXPRO_rmin(i),EXPRO_rmaj(i),EXPRO_q(i),EXPRO_kappa(i) 
        enddo
     case(6)
        do i=1,nx
           write(1,10) EXPRO_delta(i),EXPRO_te(i),EXPRO_ne(i),EXPRO_z_eff(i),EXPRO_w0(i) 
        enddo
     case(11)
        do i=1,nx
           write(1,10) EXPRO_flow_mom(i),EXPRO_pow_e(i),EXPRO_pow_i(i),EXPRO_pow_ei(i),EXPRO_zeta(i) 
        enddo
     case(16)
        do i=1,nx
           write(1,10) EXPRO_flow_beam(i),EXPRO_flow_wall(i),EXPRO_zmag(i),EXPRO_ptot(i),EXPRO_poloidalfluxover2pi(i) 
        enddo
     case(21)
        do i=1,nx
           write(1,10) EXPRO_ni(1,i),EXPRO_ni(2,i),EXPRO_ni(3,i),EXPRO_ni(4,i),EXPRO_ni(5,i)
        enddo
     case(26)
        do i=1,nx
           write(1,10) EXPRO_ti(1,i),EXPRO_ti(2,i),EXPRO_ti(3,i),EXPRO_ti(4,i),EXPRO_ti(5,i)
        enddo
     case(31)
        do i=1,nx
           write(1,10) EXPRO_vtor(1,i),EXPRO_vtor(2,i),EXPRO_vtor(3,i),EXPRO_vtor(4,i),EXPRO_vtor(5,i)
        enddo
     case(36)
        do i=1,nx
           write(1,10) EXPRO_vpol(1,i),EXPRO_vpol(2,i),EXPRO_vpol(3,i),EXPRO_vpol(4,i),EXPRO_vpol(5,i)
        enddo
     end select

  enddo

  close(1)

10 format(5(1pe14.7,2x))
20 format(30(a))
25 format(a,i2)
30 format(a,i3)
40 format(a,i6)
60 format(a,1pe13.7)

  call EXPRO_alloc('./',0) 

end program expromake
