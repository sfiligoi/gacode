module cgyro_allocate_profile

  implicit none

  public :: PROFILE_SIM_alloc, PROFILE_EXP_alloc

  logical, private :: initialized_sim = .false.
  logical, private :: initialized_exp = .false.

contains

  !  Allocate/Deallocate simulation-grid profile functions.
  subroutine PROFILE_SIM_alloc(flag)

    use cgyro_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if(flag == 1) then
       if(initialized_sim) return

       geo_numeq_flag = 0
       geo_ny = 0
       allocate(geo_yin(8,0:geo_ny))
       geo_yin(:,:) = 0.0

       initialized_sim = .true.

    else

       if(.NOT. initialized_sim) return

       deallocate(geo_yin)

       initialized_sim = .false.

    endif

  end subroutine PROFILE_SIM_alloc

  !  Allocate/Deallocate experimental-grid profile functions 
  !  (all _exp and _p) variables.
  subroutine PROFILE_EXP_alloc(flag)

    use cgyro_profile_exp
    use cgyro_globals, only: geo_ny
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if(flag == 1) then
       if(initialized_exp) return

       allocate(rhoN_torflux_exp(n_grid_exp))
       allocate(psiN_polflux_exp(n_grid_exp))
       allocate(rmin_exp(n_grid_exp))
       allocate(rmaj_exp(n_grid_exp))
       allocate(q_exp(n_grid_exp))
       allocate(kappa_exp(n_grid_exp))
       allocate(delta_exp(n_grid_exp))
       allocate(zeta_exp(n_grid_exp))
       allocate(zmag_exp(n_grid_exp))
       allocate(te_ade_exp(n_grid_exp))
       allocate(ne_ade_exp(n_grid_exp))
       allocate(tem_exp(n_species_exp,n_grid_exp))
       allocate(den_exp(n_species_exp,n_grid_exp))
       allocate(r_p(n_grid_exp))
       allocate(rmaj_p(n_grid_exp))
       allocate(zmag_p(n_grid_exp))
       allocate(b_unit_p(n_grid_exp))
       allocate(shift_p(n_grid_exp))
       allocate(s_kappa_p(n_grid_exp))
       allocate(s_delta_p(n_grid_exp))
       allocate(s_zeta_p(n_grid_exp))
       allocate(s_zmag_p(n_grid_exp))
       allocate(shat_p(n_grid_exp))
       allocate(dlnndr_p(n_species_exp,n_grid_exp))
       allocate(dlntdr_p(n_species_exp,n_grid_exp))     
       
       allocate(geo_yin_exp(8,0:geo_ny,n_grid_exp))
       geo_yin_exp(:,:,:)=0.0

       initialized_exp = .true.

    else
       if(.NOT. initialized_exp) return

       deallocate(rhoN_torflux_exp)
       deallocate(psiN_polflux_exp)
       deallocate(rmin_exp)
       deallocate(rmaj_exp)
       deallocate(q_exp)
       deallocate(kappa_exp)
       deallocate(delta_exp)
       deallocate(zeta_exp)
       deallocate(zmag_exp)
       deallocate(te_ade_exp)
       deallocate(ne_ade_exp)
       deallocate(tem_exp)
       deallocate(den_exp)
       deallocate(r_p)
       deallocate(rmaj_p)
       deallocate(zmag_p)
       deallocate(b_unit_p)
       deallocate(shift_p)
       deallocate(s_kappa_p)
       deallocate(s_delta_p)
       deallocate(s_zeta_p)
       deallocate(s_zmag_p)
       deallocate(shat_p)
       deallocate(dlnndr_p)
       deallocate(dlntdr_p)   
 
       deallocate(geo_yin_exp)

       initialized_exp = .false.

    endif

  end subroutine PROFILE_EXP_alloc


end module cgyro_allocate_profile
