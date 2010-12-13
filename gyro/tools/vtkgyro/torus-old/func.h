void func(float ri, float theta, float phi, float *s, float *x)
{
  float bigR;

  complex<float> f_c;
  complex<float> phase;
  complex<float> cs;
  complex<float> IC(0.0,1.0);
  float dx;
  float sum;

  int i_loc;
  int i1,i2;
  int iL,iR;
  int wrap;

  //---------------------------------------------------
  // Shaped equilibrium:

  bigR = R0[i_r]+ri*cos(theta+asin(delta[i_r])*sin(theta));

  x[0] = bigR*sin(phi);
  x[1] = bigR*cos(phi);
  x[2] = kappa[i_r]*ri*sin(theta);

  //---------------------------------------------------


  i1 = (int) (n_theta_plot*(theta+pi)/(2*pi));
  i2 = i1+1;  

  wrap = 0;
  if (i2 >= n_theta_plot) {
    i1 = n_theta_plot-1;
    i2 = 0;
    wrap = 1;
  }

  sum = 0.0;

  for (int i_n=i_n_0; i_n<n_n; i_n++) {

    i_loc = i_r*(n_theta_plot)+ 
      i_field*(n_r*n_theta_plot) +
      i_n*(n_field*n_r*n_theta_plot) + 
      i_time*(n_n*n_field*n_r*n_theta_plot);

    phase = 1.0;

    if (wrap == 1) {
      phase = exp(-(float)(2.0*pi*n[i_n]*q[i_r])*IC);
    }

    dx = (theta-theta_plot[i1])/d_theta_plot;

    if (n_theta_plot > 1) {

    iL = i1+i_loc;
    iR = i2+i_loc;
 
    f_c = (float)((1.0-dx))*(f_real[iL]+IC*f_imag[iL]) 
      +dx*(f_real[iR]+IC*f_imag[iR])*phase;

    } else {

    f_c = (f_real[i_loc]+IC*f_imag[i_loc])*pow(phase,-dx); 
     
   }

   cs = exp(-(float)(n[i_n])*IC*(phi-q[i_r]*theta));

   if (i_n > 0) {
      sum = sum+real(cs*f_c);
    } else {
      sum = sum+0.5*real(cs*f_c);
    }

  }

  *s = 0.5+function_magnification*sum;
}
