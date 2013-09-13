!odule eqi_geq_mod
 
  use ezspline_obj
  use ezspline
  use geqdsk_mod
 
  implicit NONE
  save
 
  type(ezspline2_r8) :: pspl   ! psi(R,Z) spline
 
  type(geqdsk) :: geq          ! GEQDISK data
 
  real*8 r0,z0                 ! mag. axis
 
  real*8 rpp,zpp               ! a point on large R side just beyond bdy
 
  real*8 :: iter_adjust = 0.9975_geq_r8
 
!nd module eqi_geq_mod
