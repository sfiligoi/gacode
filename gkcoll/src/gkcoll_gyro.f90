module gkcoll_gyro
  implicit none

  public :: GYRO_alloc
  logical, private :: initialized = .false.
  real, dimension(:,:,:,:,:), allocatable :: gyrop_J0
  real, dimension(:,:,:,:,:), allocatable :: gyrox_J0

contains

  subroutine GYRO_alloc(flag)
    use gkcoll_globals
    use gkcoll_equilibrium
    use gkcoll_legendre
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    real, external :: BESJ0
    integer :: is,ir,it,ie,ix,jx
    real :: arg, val
    real, dimension(:), allocatable :: jxi, w_jxi
    integer, parameter :: n_jxi=100

    if(flag == 1) then
       if(initialized) return

       allocate(gyrop_J0(n_species,n_radial,n_theta,n_energy,n_xi))
       allocate(jxi(n_jxi))
       allocate(w_jxi(n_jxi))
       call gauss_legendre(-1.0,1.0,jxi,w_jxi,n_jxi)

       do is=1,n_species
          do ir=1,n_radial
             do it=1,n_theta
                do ie=1,n_energy
                   do ix=1,n_xi
                      gyrop_J0(is,ir,it,ie,ix) = 0.0
                      do jx=1,n_jxi
                         arg = k_perp(it,ir) * rho * vth(is) &
                              * mass(is) / (z(is) * Bmag(it)) &
                              * sqrt(2.0* energy(ie)) &
                              * sqrt(1-jxi(jx)**2)
                         call legendre(indx_xi(ix),jxi(jx),val)
                         gyrop_J0(is,ir,it,ie,ix) &
                              = gyrop_J0(is,ir,it,ie,ix) &
                              + 0.5 * w_jxi(jx) * BESJ0(arg) * val
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       deallocate(jxi)
       deallocate(w_jxi)

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
       deallocate(gyrop_J0)
       deallocate(gyrox_J0)
       initialized = .false.
    endif
    
  end subroutine gyro_alloc

end module gkcoll_gyro
