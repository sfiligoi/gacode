subroutine p_cub_spline_deriv(x,y,n,yp)

  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: x,y
  double precision, intent(inout), dimension(n) :: yp
  
!f2py intent(in) n
!f2py intent(in) x
!f2py intent(in) y
!f2py intent(out) yp

  call cub_spline_deriv(x,y,n,yp)
  
end subroutine p_cub_spline_deriv

subroutine p_bound_deriv(df,f,r,n)

  integer, intent(in) :: n

  double precision, intent(inout), dimension(n) :: df
  double precision, intent(in), dimension(n) :: f
  double precision, intent(in),dimension(n) :: r

!f2py intent(in) n
!f2py intent(in) f
!f2py intent(in) r
!f2py intent(out) df

  call bound_deriv(df,f,r,n)
  
end subroutine p_bound_deriv
