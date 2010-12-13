void func_sq2(float rp, float alpha, float *s, float *x)
{
  complex<float> f_c;
  complex<float> f_c1;
  complex<float> f_c2;
  complex<float> cs;
  complex<float> IC(0.0,1.0);
  float sum;
  float dx;

  int i_loc;
  int i_loc1;
  int i_loc2;
  int i1;
  int i2;

  //---------------------------------------------------
  x[0] = 0.0;
  x[1] = 10.0*(rp-r[0])/(r[n_r-1]-r[0])-1.01;
  x[2] = 10.0*alpha/((n_alpha-1)*d_alpha)-0.5;
  //---------------------------------------------------

  sum = 0.0;

  for (int i_n=i_n_0; i_n<n_n; i_n++) {

    i1 = (int) ( n_r*(rp-r[0])/(r[n_r-1]-r[0]) );
    i2 = i1+1;

    if (i2 >= n_r) i2=0;

    i_loc1 = i1*(n_theta_plot)+ 
      i_field*(n_r*n_theta_plot) +
      i_n*(n_field*n_r*n_theta_plot);

    i_loc2 = i2*(n_theta_plot)+ 
      i_field*(n_r*n_theta_plot) +
      i_n*(n_field*n_r*n_theta_plot);

    f_c1 = (f_real[i_loc1]+IC*f_imag[i_loc1]); 
    f_c2 = (f_real[i_loc2]+IC*f_imag[i_loc2]); 
    
    dx = (rp-r[i1])/(r[2]-r[1]);    

    f_c = (float)((1.0-dx))*f_c1+dx*f_c2;

    cs = exp(-(float)(n[i_n])*IC*alpha);

    if (i_n > 0) {
      sum = sum+real(cs*f_c);
    } else {
      sum = sum+0.5*real(cs*f_c);
    }

  }

  *s = 0.5+2*function_magnification*sum;
}
