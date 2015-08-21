subroutine tgyro_residual(f,g,res,n,method)

  use tgyro_globals, only : n_r, tgyro_resnorm, &
       loc_ti_feedback_flag,&
       loc_te_feedback_flag,&
       loc_ne_feedback_flag,&
       loc_er_feedback_flag,&
       loc_he_feedback_flag, q_gb,gamma_gb,pi_gb

  use tgyro_iteration_variables

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: method
  real, intent(in), dimension(n) :: f
  real, intent(in), dimension(n) :: g
  real, dimension(n) :: fw
  real, dimension(n) :: gw
  real, intent(inout), dimension(n) :: res

  p = 0
  do i=2,n_r
     if (loc_ti_feedback_flag == 1) then
        p = p+1
        weight(p) = q_gb(i)/q_gb(2)
     endif
     if (loc_te_feedback_flag == 1) then
        p = p+1
        weight(p) = q_gb(i)/q_gb(2)
     endif
     if (loc_ne_feedback_flag == 1) then
        p = p+1
        weight(p) = gamma_gb(i)/gamma_gb(2)
     endif
     if (loc_er_feedback_flag == 1) then
        p = p+1
        weight(p) = pi_gb(i)/pi_gb(2)
     endif
     if (loc_he_feedback_flag == 1) then
        p = p+1
        weight(p) = gamma_gb(i)/gamma_gb(2)
     endif
  enddo

  if (tgyro_resnorm == 1) then
     fw = f
     gw = g
  else
     fw = f*weight
     gw = g*weight
  endif

  select case (method)

  case (1) 

     ! ORIGINAL METHOD:
     res = (fw-gw)**2/(fw**2+gw**2)

  case (2)

     ! SIMPLE NORM:
     res = abs(fw-gw)

  case (3)

     ! SQUARE RESIDUAL
     res = 0.5*(fw-gw)**2

  case (4)

     ! BALANCED
     res = (fw-gw)**2/MAX((fw**2+gw**2),1.0)

  case (5)

     ! WEIGHTED
     res(:) = 0.0
     do p=1,n
        if (quant(p) == 'ne') then 
           res(p) = 0.005*(fw(p)-gw(p))**2
        else
           res(p) = 0.5*(fw(p)-gw(p))**2
        endif
     enddo

  case default

     print *,'Error in tgyro_residual'
     stop

  end select

end subroutine tgyro_residual
