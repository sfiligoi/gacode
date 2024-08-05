!\***************************************************************
! pt_mmm_mod.f90: pt_solver implementation of mmm
! series turbulent model.
!/***************************************************************
MODULE tgyro_mmm_mod
  implicit none

  ! public available varaibles
  PUBLIC
  real   :: mmm_chi_i       !Ion thermal diffusivity
  real   :: mmm_chi_e       !Electron thermal diffusivity
  real   :: mmm_chi_ne      !Electron particle diffusivity
  real   :: mmm_chi_nx      !Impurity ion diffusivity
  real   :: mmm_chi_phi     !Toroidal momentum diffusivity
  real   :: mmm_chi_theta   !Poloidal momentum diffusivity

  real   :: mmm_vheat_i     !Electron thermal
  real   :: mmm_vheat_e     !Ion thermal
  real   :: mmm_vgx_ne      !Particle
  real   :: mmm_vgx_nx      !Impurity particle
  real   :: mmm_vgx_phi     !Toroidal momentum
  !
  ! function accessibility
  PUBLIC :: tgyro_mmm_map
  PUBLIC :: tgyro_mmm_run

CONTAINS
  !\------------------------------------------
  ! load data for mmm mode
  !/
  SUBROUTINE tgyro_mmm_map(ierr)
    implicit none
    ! output
    integer, intent(out) :: ierr

    print *, "TGYRO is compiled without the MMM model"
    ierr = 1
    return
  END SUBROUTINE tgyro_mmm_map

  !\
  ! run mmm model
  !/
  SUBROUTINE tgyro_mmm_run(ierr)
    implicit none
    ! output
    integer, intent(out) :: ierr

    print *, "TGYRO is compiled without the MMM model"
    ierr = 1
    mmm_chi_i     = 0d0      ! Ion thermal diffusivity
    mmm_chi_e     = 0d0      ! Electron thermal diffusivity
    mmm_chi_ne    = 0d0      ! Electron particle diffusivity
    mmm_chi_nx    = 0d0      ! Impurity ion diffusivity
    mmm_chi_phi   = 0d0      ! Toroidal momentum diffusivity
    mmm_chi_theta = 0d0      ! Poloidal momentum diffusivity

    mmm_vheat_i = 0d0        ! Electron thermal
    mmm_vheat_e = 0d0        ! Ion thermal
    mmm_vgx_ne  = 0d0        ! Particle
    mmm_vgx_nx  = 0d0        ! Impurity particle
    mmm_vgx_phi = 0d0        ! Toroidal momentum
    return
  END SUBROUTINE tgyro_mmm_run

END MODULE tgyro_mmm_mod
