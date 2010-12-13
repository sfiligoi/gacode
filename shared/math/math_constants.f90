module math_constants

  implicit none

  real, parameter ::  &
     root_1_by_pi = 5.641895835477563e-1, &
        two_by_pi = 6.366197723675813e-1, &
          pi_by_4 = 7.853981633974483e-1, &
     root_2_by_pi = 7.978845608028654e-1, &
          pi_by_2 = 1.570796326794897e0,  &
        pi_3_by_4 = 2.356194490192345e0,  &
               pi = 3.141592653589793e0,  &
             pi_2 = 6.283185307179586e0,  &
        twothirds = 0.666666666666667e0

  complex, parameter :: i_c = (0.0,1.0)

end module math_constants
