module gkcoll_neo

  implicit none

  public :: NEO_alloc, NEO_do

  logical, private :: initialized = .false.
  integer, parameter, private :: io_neo = 1
  character(len=80), private :: runfile_neo   = 'out.gkcoll.neo'
  complex, dimension(:), allocatable, private :: pflux, eflux, uparB
  
contains

  subroutine NEO_alloc(flag)
    use gkcoll_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
  
    if(neoclassical_model /= 1) return

    if(flag == 1) then
       if(initialized) return

       allocate(pflux(n_species))
       allocate(eflux(n_species))
       allocate(uparB(n_species))

       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io_neo,file=trim(path)//runfile_neo,status='replace')
          close(io_neo)
        endif

       initialized = .true.

    else
       if(.NOT. initialized) return

       deallocate(pflux)
       deallocate(eflux)
       deallocate(uparB)

       initialized = .false.
    endif

  end subroutine NEO_alloc

  subroutine NEO_do
    use gkcoll_globals
    use gkcoll_equilibrium
    implicit none
    integer :: is,it,ie,ix
    integer :: ir=1

    if(neoclassical_model /= 1) return

    do is=1,n_species
       pflux(is) = (0.0,0.0)
       eflux(is) = (0.0,0.0)
       uparB(is) = (0.0,0.0)
       do it=1,n_theta
          do ie=1,n_energy
             do ix=1,n_xi
         
                if(indx_xi(ix) == 0) then
                   pflux(is) = pflux(is) + cap_h_p(is,ir,it,ie,ix) &
                        * omega_rdrift(it,is) * energy(ie) &
                        * dens(is) * w_e(ie) * w_theta(it) * (4.0/3.0)
                   eflux(is) = eflux(is) + cap_h_p(is,ir,it,ie,ix) &
                        * omega_rdrift(it,is) * energy(ie) &
                        * energy(ie) * temp(is) &
                        * dens(is) * w_e(ie) * w_theta(it) * (4.0/3.0)
                else if(indx_xi(ix) == 2) then
                   pflux(is) = pflux(is) + cap_h_p(is,ir,it,ie,ix) &
                        * omega_rdrift(it,is) * energy(ie) &
                        * dens(is) * w_e(ie) * w_theta(it) * (2.0/15.0)
                   eflux(is) = eflux(is) + cap_h_p(is,ir,it,ie,ix) &
                        * omega_rdrift(it,is) * energy(ie) &
                        * energy(ie) * temp(is) &
                        * dens(is) * w_e(ie) * w_theta(it) * (2.0/15.0)
                else if(indx_xi(ix) == 1) then
                   uparB(is) = uparB(is) + cap_h_p(is,ir,it,ie,ix) &
                        * sqrt(2.0*energy(ie)) * vth(is) * Bmag(it) &
                        * w_e(ie) * w_theta(it) * (1.0/3.0)
                endif

             enddo
          enddo
       enddo
    enddo

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neo,file=trim(path)//runfile_neo,status='old',&
            position='append')
       do is=1,n_species
          write(io_neo,20) real(pflux(is)), real(eflux(is)), real(uparB(is))
       enddo
       close(io_neo)
       print *, real(pflux(1)), real(eflux(1)), real(uparB(1))
    endif
    
20  format(3(es11.4,1x))

  end subroutine NEO_do

end module gkcoll_neo
