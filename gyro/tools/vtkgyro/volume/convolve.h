/*==========================================================
  module: convolve.h
 
  Compute the convolution f = (f1,f2).

  Here, f1 -> f1_real + i f1_imag  
        f2 -> f2_real + i f2_imag

  The result is f_real + i f_imag -> f. 
  ==========================================================*/

void convolve( void ) {    

  int i_loc;
  int i_m;

  complex<float> IC(0.0,1.0);
  complex<float> *f_n;
  complex<float> *f1_n;
  complex<float> *f2_n;
  complex<float> p1;
  complex<float> p2;


  f_n  = new complex<float>[n_n];
  f1_n = new complex<float>[n_n];
  f2_n = new complex<float>[n_n];

  for (i_time=0; i_time<n_time; i_time++) {
    for (i_r=0; i_r<n_r; i_r++) {
      for (i_theta=0; i_theta<n_theta_plot; i_theta++) {

	for (i_n=0; i_n<n_n; i_n++) {
 
	  i_loc = i_theta+i_r*(n_theta_plot)+ 
	    i_field*(n_r*n_theta_plot) +
	    i_n*(n_field*n_r*n_theta_plot) + 
	    i_time*(n_n*n_field*n_r*n_theta_plot);

          // (i k_rho) * conjg(phi)

	  f1_n[i_n] = IC*(float)(n[i_n])*kr[i_r]*(f1_real[i_loc]-IC*f1_imag[i_loc]); 

          // moment_n or moment_e

	  f2_n[i_n] = f2_real[i_loc]+IC*f2_imag[i_loc]; 

	}

	for (i_n=0; i_n<n_n; i_n++) {
	  f_n[i_n] = 0.0;
	  for (i_m=-n_n+1; i_m<n_n; i_m++) {

            if (abs(i_n-i_m) < n_n) {
 
	      if (i_n >= i_m) {
		p1 = f1_n[i_n-i_m];
	      } else {
		p1 = conj(f1_n[i_m-i_n]);
	      }

	      if (i_m >= 0) {
		p2 = f2_n[i_m];
	      } else {
		p2 = conj(f2_n[i_m]);
	      }
 
	      f_n[i_n] = f_n[i_n]+p1*p2; 

	    }
	  }
	}

	for (i_n=0; i_n<n_n; i_n++) {
 
	  i_loc = i_theta+i_r*(n_theta_plot)+ 
	    i_field*(n_r*n_theta_plot) +
	    i_n*(n_field*n_r*n_theta_plot) + 
	    i_time*(n_n*n_field*n_r*n_theta_plot);

	  f_real[i_loc] = real(f_n[i_n]); 
	  f_imag[i_loc] = imag(f_n[i_n]); 

	}

      }
    }
  }

  delete [] f_n;
  delete [] f1_n;
  delete [] f2_n;

}
