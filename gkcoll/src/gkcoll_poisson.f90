module gkcoll_poisson
  
  implicit none
  
  public :: POISSONp_do, POISSONx_do
  
contains
  
  subroutine POISSONp_do(ir,it)
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    implicit none
    integer, intent (in) :: ir, it
    integer :: is, ie, ix
    real :: sum_den
    
    phi(ir,it) = (0.0,0.0)
    do is=1,n_species
       do ie=1,n_energy
          do ix=1,n_xi
             phi(ir,it) = phi(ir,it) &
                  + gyrop_J0(is,ir,it,ie,ix) &
                  * z(is)*dens(is) * w_e(ie) &
                  * cap_h_p(is,ir,it,ie,ix)
          enddo
       enddo
    enddo

    sum_den = 0.0
    do is=1,n_species
       do ie=1,n_energy
          do ix=1,n_xi
             sum_den = sum_den &
                  + 0.5 * w_xi(ix) &
                  * z(is)**2/temp(is) *dens(is) * w_e(ie)
          enddo
       enddo
    enddo

    if(adiabatic_ele_model == 1) then
       sum_den = sum_den + dens_ele / temp_ele
    endif

    phi(ir,it) = phi(ir,it) / (-k_perp(it,ir)**2 * lambda_debye**2 &
         * dens_ele / temp_ele + sum_den) 
    
  end subroutine POISSONp_do
  
  subroutine POISSONx_do(ir,it)
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    implicit none
    integer, intent (in) :: ir, it
    integer :: is, ie, ix
    real :: sum_den
    
    phi(ir,it) = (0.0,0.0)
    do is=1,n_species
       do ie=1,n_energy
          do ix=1,n_xi
             phi(ir,it) = phi(ir,it) &
                  + 0.5 * w_xi(ix) &
                  * gyrox_J0(is,ir,it,ie,ix) &
                  * z(is)*dens(is) * w_e(ie) &
                  * h_x(is,ir,it,ie,ix)
          enddo
       enddo
    enddo
    
    sum_den = 0.0
    do is=1,n_species
       do ie=1,n_energy
          do ix=1,n_xi
             sum_den = sum_den &
                  + 0.5 * w_xi(ix) &
                  * (1.0 - gyrox_J0(is,ir,it,ie,ix)**2) &
                  * z(is)**2/temp(is) *dens(is) * w_e(ie)
          enddo
       enddo
    enddo

    if(adiabatic_ele_model == 1) then
       sum_den = sum_den + dens_ele / temp_ele
    endif

    phi(ir,it) = phi(ir,it) / (-k_perp(it,ir)**2 * lambda_debye**2 &
         * dens_ele / temp_ele + sum_den) 

  end subroutine POISSONx_do

end module gkcoll_poisson
