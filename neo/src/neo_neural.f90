module neo_neural

  implicit none

  public :: NEURAL_alloc, NEURAL_do, NEURAL_write
  real :: jpar_nn_neo, jtor_nn_neo, ke_nn_neo, ki1_nn_neo, ki2_nn_neo
  character(len=80),private :: runfile_nn = 'out.neo.transport_nn'
  logical, private :: initialized = .false.

contains

  subroutine NEURAL_alloc(flag)
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
  end subroutine NEURAL_alloc

  subroutine NEURAL_do(ir)
    implicit none
    integer, intent (in) :: ir
    jpar_nn_neo = 0.0
    jtor_nn_neo = 0.0
    ke_nn_neo   = 0.0
    ki1_nn_neo  = 0.0
    ki2_nn_neo  = 0.0
  end subroutine NEURAL_do
  
  subroutine NEURAL_write(ir)
    implicit none
    integer, intent (in) :: ir
  end subroutine NEURAL_write
  
end module neo_neural
