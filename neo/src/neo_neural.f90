module neo_neural

  implicit none

  public :: NEURAL_alloc, NEURAL_do, NEURAL_write
  real :: jpar_nn_neo, jtor_nn_neo
  real, dimension(:), allocatable :: kbig_upar_nn_neo, vpol_th0_nn_neo
  character(len=80),private :: runfile_nn = 'out.neo.transport_nn'
  logical, private :: initialized = .false.

contains

  subroutine NEURAL_alloc(flag)
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    if(flag == 1) then
       if(initialized) return
       allocate(kbig_upar_nn_neo(n_species))
       allocate(vpol_th0_nn_neo(n_species))
       initialized = .true.
    else
       if(.NOT. initialized) return
       deallocate(kbig_upar_nn_neo)
       deallocate(vpol_th0_nn_neo)
       initialized = .false.
    end if
  end subroutine NEURAL_alloc

  subroutine NEURAL_do(ir)
    implicit none
    integer, intent (in) :: ir
    jpar_nn_neo = 0.0
    jtor_nn_neo = 0.0
    kbig_upar_nn_neo(:) = 0.0
    vpol_th0_nn_neo(:)  = 0.0
  end subroutine NEURAL_do
  
  subroutine NEURAL_write(ir)
    implicit none
    integer, intent (in) :: ir
  end subroutine NEURAL_write
  
end module neo_neural
