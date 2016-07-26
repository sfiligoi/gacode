subroutine tgyro_eped_nn

  use tgyro_ped
  
  implicit none

  character(len=1000) :: nn_executable
  character(len=1000) :: nn_files
  character(len=1) :: dummy

  character(len=1000) :: epednn_model
  real*4 :: INPUT_PARAMETERS(10)
  real*4 :: OUTPUT_PARAMETERS(5)
  integer :: ierr
  CHARACTER NUL
  PARAMETER(NUL = CHAR(0))

  include 'brainfuse_lib.inc'

  INPUT_PARAMETERS( 1)= a_in
  INPUT_PARAMETERS( 2)= betan_in
  INPUT_PARAMETERS( 3)= bt_in
  INPUT_PARAMETERS( 4)= delta_in
  INPUT_PARAMETERS( 5)= ip_in
  INPUT_PARAMETERS( 6)= kappa_in
  INPUT_PARAMETERS( 7)= m_in
  INPUT_PARAMETERS( 8)= neped_in
  INPUT_PARAMETERS( 9)= r_in
  INPUT_PARAMETERS(10)= zeffped_in
  
  WRITE(*,*)INPUT_PARAMETERS
  
  call get_environment_variable('TGLFNN_MODEL_DIR',epednn_model)
  write(*,*)TRIM(epednn_model)
  ierr=load_anns(1, TRIM(epednn_model)//NUL,'brainfuse'//NUL)
  ! ierr=load_anns_inputs(INPUT_PARAMETERS)
  ! ierr=run_anns()
  ! ierr=get_anns_avg_array(OUTPUT_PARAMETERS)
  
  ! WRITE(*,*)OUTPUT_PARAMETERS
  
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
       neped_in   ,& ! at ped
       r_in       ,& 
       zeffped_in    ! at ped
  close(1)
  !
  ! Execute the NN
  call get_environment_variable('EPEDNN_EXEC',nn_executable)
  call get_environment_variable('EPEDNN_MODEL',nn_files)
  write(*,*)'EPEDNN_EXEC: ',trim(nn_executable)
  write(*,*)'EPEDNN_MODEL: ',trim(nn_files)
  call gacode_system(trim(nn_executable)//' '//trim(nn_files)//' input.dat')

  ! Read epednn_vec=[psi_norm,ne,ptot=2 ne T]
  open(unit=1,file='epednn.profiles',status='old')
  read(1,*) dummy
  read(1,*) nn_vec
  close(1)
  ! Correct n
  nn_vec(nx_nn,2) = 2*nn_vec(nx_nn-1,2)-nn_vec(nx_nn-2,2)
  
  ! NOTE:  w_ped currently fixed so that psi_top = 0.9 = 1-2.5*w_ped
  !        OR w_ped = 0.1/2.5

  nn_w_ped = 0.1/2.5
  
end subroutine tgyro_eped_nn
