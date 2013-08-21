subroutine set_nbdims

  use nbi_types
  USE transp
  USE nbi_dimensions
  implicit NONE


  ! initialize dimensions in NBI calculation

  type (nbitype_dims) :: zdims        ! NBI code dimensions

  integer ierr

  !--------------------------------------------------------------

  call nbi_init_dims(zdims)           ! get defaults,source in nbi_init.f90

  ! set most dimensions (leave a few defaulted)...
  !undecalred ones are mibs,mmini,mifpsh   HSJ

  zdims%mib = beam_data%nbeam         ! max # of beams... (min of 1,0 not allowed)
  zdims%mj  = beam_data%nzone_nb + 3  ! # of radial zones + at least 3 ?
  zdims%mig = beam_data%ngmax         ! max # of thermal species not the same
  zdims%mibs = beam_data%nfast        ! max # of fast species (incl. minorities)
  zdims%mmini= nmini                  ! max # of minority species
  zdims%miefi= max(1,nerngfi)         ! # of energy ranges for some outputs

  zdims%mimxba = mimxba               ! fbm pitch angle dimension
  zdims%mimxbe = mimxbe               ! fbm energy dimension

! the following are not used ??
!  zdims%mifpen = mifpen               ! fpp energy grid (plasma frame)
!  zdims%mifpxi = mifpxi               ! fpp pitch angle grid (plasma frame)
!  nzone_fp = 50
!  zdims%mifpsh = nzone_fp             ! # of fpp radial zones

!  zdims%miboxm = miboxm               ! max # of boxes, beam-in-box calculation

  ! apply...

  call nbi_set_dims(zdims,ierr)       ! set array sizes in NBI code,
                                      ! source in nbi_set.f90
  if(ierr.ne.0) call errmsg_exit(' ?set_nbdims: nbi_set_dims call failed.')

  ! and allocate the arrays...
  ! most quantities cleared to 0,0.0,.FALSE.

  call nbi_alloc_init                 !source in nbi_alloc.f90

  ! note there are a few arrays whose allocation is deferred; sizes will
  ! be set from input values... (cf call to nbi_alloc2)
!  print *,'sub_nbdims,zdims ='
!  print *,'zdims%mib'    , zdims%mib
!  print *,'zdims%mj'     , zdims%mj
!  print *,'zdims%mig'    , zdims%mig
!  print *,'zdims%mibs'   , zdims%mibs
!  print *,'zdims%mimpt'  , zdims%mimpt
!  print *,'zdims%mmini'  , zdims%mmini  
!  print *,'zdims%mtheta' , zdims%mtheta 
!  print *,'zdims%miefi'  , zdims%miefi
!  print *,'zdims%mimxba' , zdims%mimxba
!  print *,'zdims%mimxbe' , zdims%mimxbe
!  print *,'zdims%mifpen' , zdims%mifpen
!  print *,'zdims%mifpxi' , zdims%mifpxi
!  print *,'zdims%mifpsh' , zdims%mifpsh
!  print *,'zdims%miboxm' , zdims%miboxm
end subroutine set_nbdims
