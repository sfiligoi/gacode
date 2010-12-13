void c_func(float ri, float theta, float phi, float *s, float *x)
{
  float bigR;

  //---------------------------------------------------
  // Shaped equilibrium:

  bigR = R0_halo+ri*cos(theta+asin(delta_halo)*sin(theta));

  x[0] = bigR*sin(phi);
  x[1] = bigR*cos(phi);
  x[2] = kappa_halo*ri*sin(theta);

  //---------------------------------------------------

  *s = 1.0;
}
