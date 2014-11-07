module cgyro_gyro
  implicit none

  public :: GYRO_alloc
  logical, private :: initialized = .false.
  real, dimension(:,:,:,:,:), allocatable :: gyrox_J0

contains

  subroutine GYRO_alloc(flag)
    use cgyro_globals
    use cgyro_equilibrium
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    real, external :: BESJ0
    integer :: is,ir,it,ie,ix
    real :: arg

    if (flag == 1) then
       if (initialized) return

       allocate(gyrox_J0(n_species,n_radial,n_theta,n_energy,n_xi))
       do is=1,n_species
          do ir=1,n_radial
             do it=1,n_theta
                do ie=1,n_energy
                   do ix=1,n_xi
                      arg = k_perp(it,ir) * rho * vth(is) &
                           * mass(is) / (z(is) * Bmag(it)) &
                           * sqrt(2.0* energy(ie)) &
                           * sqrt(1-xi(ix)**2)
                      gyrox_J0(is,ir,it,ie,ix) = BESJ0(arg)
                   enddo
                enddo
             enddo
          enddo
       enddo

       initialized = .true.

    else
       if(.NOT. initialized) return
       deallocate(gyrox_J0)
       initialized = .false.
    endif
    
  end subroutine gyro_alloc

end module cgyro_gyro
