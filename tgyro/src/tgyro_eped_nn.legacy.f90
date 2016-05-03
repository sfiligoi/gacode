subroutine tgyro_eped_nn

  use tgyro_ped
  
  implicit none

  character(len=1000) :: nn_executable
  character(len=1000) :: nn_files
  character(len=1) :: dummy

  ! Write input file for the NN
  open(unit=1,file='input.dat',status='replace')
  write(1,'(a)') '1'
  write(1,'(10(f6.3,1x))') &
       a_in       ,&
       betan_in   ,&
       bt_in      ,&
       delta_in   ,&
       ip_in      ,&
       kappa_in   ,&
       m_in       ,&
       neped_in   ,&
       r_in       ,&
       zeffped_in 
  close(1)
  !
  ! Execute the NN
  call get_environment_variable('EPEDNN_EXEC',nn_executable)
  call get_environment_variable('EPEDNN_MODEL',nn_files)
  !write(*,*)'EPEDNN executable: ',trim(nn_executable)
  !write(*,*)'EPEDNN model: ',trim(nn_files)
  call system(trim(nn_executable)//' '//trim(nn_files)//' input.dat')

  ! Read epednn_vec=[psi_norm,ne,ptot=2 ne T]
  open(unit=1,file='epednn.profiles',status='old')
  read(1,*) dummy
  read(1,*) nn_vec
  close(1)

  ! Temporary -- need to read this
  nn_w_ped = 0.05

end subroutine tgyro_eped_nn
