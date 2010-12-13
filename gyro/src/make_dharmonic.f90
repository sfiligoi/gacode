!-----------------------------------------------------------
! make_dharmonic.f90
!
! PURPOSE:
!  Make pseudospectral (harmonic) derivative coefficients.
!-----------------------------------------------------------

subroutine make_dharmonic

  use gyro_globals
  use math_constants

  implicit none

  integer :: p
  complex :: z0
 

  w_d1(:) = (0.0,0.0)

  do m=-m_dx,m_dx-1
     do p=-m_dx,m_dx-1
        z0 = exp(i_c*p*m*(pi_2/n_x))
        w_d1(m) = w_d1(m)+p*z0
     enddo
  enddo

  w_d1(:) = -(1.0/n_x)*(pi_2*i_c/x_length)*w_d1(:)

  w_d2(:) = (0.0,0.0)

  do m=-m_dx,m_dx-1
     do p=-m_dx,m_dx-1
        z0 = exp(i_c*p*m*(pi_2/n_x))
        w_d2(m) = w_d2(m)+p*p*z0
     enddo
  enddo

  ! Note the sign change
  w_d2(:) = (1.0/n_x)*(pi_2*i_c/x_length)**2*w_d2(:)

end subroutine make_dharmonic

