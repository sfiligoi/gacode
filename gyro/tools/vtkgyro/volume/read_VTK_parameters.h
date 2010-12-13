/*==========================================================
  module: read_VTK_parameters.h
 
  Read INPUT in current directory.
  These are, mostly, parameters used by VTK.
  ==========================================================*/

int read_VTK_parameters( void ) {

  int DIM_LABEL=40;
  char label[DIM_LABEL];


  ifstream in("INPUT");
  if (!in) {
    cout << "Cannot open INPUT" << endl;
    return 1;
  }

  in >> vis_method >> label;
  in >> plot_index >> label;
  in >> n_theta >> label;
  in >> n_phi >> label;
  in >> phi_max >> label;
  in >> t1 >> label;
  in >> t2 >> label;
  in >> n_colors >> label;
  in >> color_exp >> label;
  in >> t_stop >> label;
  in >> n_move >> label;
  in >> view_angle >> label;
  in >> dview_angle >> label;
  in >> x_pos >> label;
  in >> y_pos >> label;
  in >> z_pos >> label;
  in >> x_foc >> label;
  in >> y_foc >> label;
  in >> z_foc >> label;
  in >> xsize >> label;
  in >> ysize >> label;
  in >> image_magnification >> label;
  in >> x_light_pos1 >> label;
  in >> y_light_pos1 >> label;
  in >> z_light_pos1 >> label;
  in >> x_light_foc1 >> label;
  in >> y_light_foc1 >> label;
  in >> z_light_foc1 >> label;
  in >> x_light_pos2 >> label;
  in >> y_light_pos2 >> label;
  in >> z_light_pos2 >> label;
  in >> x_light_foc2 >> label;
  in >> y_light_foc2 >> label;
  in >> z_light_foc2 >> label;
  in >> intensity >> label;
  in >> function_magnification >> label;

  in.close();

  cout << t1 << " " << t2 << endl;

}
