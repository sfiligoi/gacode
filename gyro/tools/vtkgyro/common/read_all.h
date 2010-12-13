/*==========================================================
  module: read_all
 
  Manage all calls to get required input data.
  ==========================================================*/

void read_all( void ) 
{
  int i;
 
  //---------------------------------------------------
  // Read profile_vugyro.out from GYRO simulation
  //
  read_vugyro(path+simdir+"/profile_vugyro.out");

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

  R0_halo = R0[n_r-1];
  kappa_halo = kappa[n_r-1];
  delta_halo = delta[n_r-1];

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

  read_time(path+simdir+"/t.out",&n_time);

  cout << "----------------------------------" << endl;
  cout << "r     = " << 
    r[0] << " " << r[n_r-1] << " " << endl;
  cout << "q     = " << 
    q[0] << " " << q[n_r-1] << " " << endl;
  cout << "R0    = " << 
    R0[0] << " " << R0[n_r-1] << " " << endl;
  cout << "delta = " << 
    delta[0] << " " << delta[n_r-1] << endl;
  cout << "kappa = " << 
    kappa[0] << " " << kappa[n_r-1] << endl;
  cout << "----------------------------------" << endl;
  cout << "Time slices in file: " << n_time << endl;
  cout << "Time slices to plot: " << t2-t1 << endl;
  cout << "----------------------------------" << endl;

  if (t2 < n_time) n_time = t2;

}
