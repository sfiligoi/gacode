/*==========================================================
  module: read_all
 
  Manage all calls to get required input data.
  ==========================================================*/

#include "convolve.h"

void read_all(string path) 
{
  int i;
 
  //---------------------------------------------------
  // Read profile_vugyro.out from GYRO simulation
  //
  read_vugyro(path+"/profile_vugyro.out");

  cout << "----------------------------------" << endl;
  cout << "n_r          = " << n_r << endl;
  cout << "n_theta_plot = " << n_theta_plot << endl;
  cout << "n_field      = " << n_field << endl;
  cout << "r0           = " << r[n_r/2+1] << endl;
  cout << "----------------------------------" << endl;

  if (boundary_method == 1) {
    cout << "boundary_method: PERIODIC" << endl;
    i_buffer = 0;
  } else {
    cout << "boundary_method: NONPERIODIC" << endl;
    cout << "   buffer width: " << i_buffer << endl;
  }
  cout << "----------------------------------" << endl;

  //----------------------------------------------------

  //r_halo = grid[4*5*n_r];
  //R0_halo = grid[1+4*5*n_r]*r_halo; 

  //q_halo = grid[2+4*5*n_r];
  //delta_halo = grid[3+4*5*n_r];
  //kappa_halo = grid[4+4*5*n_r];

  //------------------------------------------------------------
  // n and theta_plot
  //
  theta_plot = new float[n_theta_plot];
 
  d_theta_plot = 2*pi/n_theta_plot;
  if (n_theta_plot == 1) {
    theta_plot[0] = 0.0;
  } else {
    for (i=0; i<n_theta_plot; i++) theta_plot[i] = -pi+i*d_theta_plot;
  }    

  n = new int[n_n];
  for (i=0; i<n_n; i++) n[i] = n0+i*dn;

  cout << "n = " ;
  for (i=0; i<n_n; i++) cout << n[i] << " ";
  cout << endl;
  //------------------------------------------------------------

  read_time(path+"/t.out",&n_time);

  cout << "----------------------------------" << endl;
  cout << "r     = " << 
    r[0] << " " << r[n_r-1] << " " << r_halo << endl;
  cout << "q     = " << 
    q[0] << " " << q[n_r-1] << " " << q_halo << endl;
  cout << "R0    = " << 
    R0[0] << " " << R0[n_r-1] << " " << R0_halo << endl;
  cout << "delta = " << 
    delta[0] << " " << delta[n_r-1] << " " << delta_halo << endl;
  cout << "kappa = " << 
    kappa[0] << " " << kappa[n_r-1] << " " << kappa_halo << endl;
  cout << "----------------------------------" << endl;
  cout << "Time slices in file: " << n_time << endl;
  cout << "Time slices to plot: " << t2-t1 << endl;
  cout << "----------------------------------" << endl;

  if (t2 < n_time) n_time = t2;

  //-----------------------------------------------------------
  // fields
  //
  n_elements = n_theta_plot*n_r*n_field*n_n*n_time;

  // In all nonzero cases, read the potential

  switch (plot_index) {

  case 0 : 

    i_n_0 = 1;

    cout << "** Setting field to zero **" << endl;

    for (i=0;i<n_elements;i++) { 
      f_real[i] = 0.0;
      f_imag[i] = 0.0;
    }
    break; 
    
  case 1 : 

    i_n_0 = 0;

    cout << "** Reading u.out **" << endl;

    f_real = new float[n_elements];
    f_imag = new float[n_elements];

    read_complex(path+"/u.out",n_elements,f_real,f_imag);
    break;

  case 2 :

    i_n_0 = 0;

    f_real  = new float[n_elements];
    f_imag  = new float[n_elements];
    f1_real = new float[n_elements];
    f1_imag = new float[n_elements];
    f2_real = new float[n_elements];
    f2_imag = new float[n_elements];

    cout << "Reading u.out: " << endl;
    read_complex(path+"/u.out",n_elements,f1_real,f1_imag);
    cout << "Reading moment_n.out: " << endl;
    read_complex(path+"/moment_n.out",n_elements,f2_real,f2_imag);

    // Convolve f1 and f2:

    convolve();

    delete [] f1_real;
    delete [] f1_imag;
    delete [] f2_real;
    delete [] f2_imag;

    break;

  case 3 :

    i_n_0 = 0;

    f_real  = new float[n_elements];
    f_imag  = new float[n_elements];
    f1_real = new float[n_elements];
    f1_imag = new float[n_elements];
    f2_real = new float[n_elements];
    f2_imag = new float[n_elements];

    cout << "Reading u.out: " << endl;
    read_complex(path+"/u.out",n_elements,f1_real,f1_imag);
    cout << "Reading moment_e.out: " << endl;
    read_complex(path+"/moment_e.out",n_elements,f2_real,f2_imag);

    // Convolve f1 and f2:

    convolve();

    delete [] f1_real;
    delete [] f1_imag;
    delete [] f2_real;
    delete [] f2_imag;

    break;

  default : 

    cout << "Bad plot_index." << endl; 
    exit(1);
    break;

  }

  //--------------------------------------------------------------

}
