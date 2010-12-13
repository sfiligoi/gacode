/*==========================================================
  module: globals.h
 
  Globals variables for VTKanim cap/torus rendering.
  ==========================================================*/

// variables contained in INPUT

int plot_index;

int t1;
int t2;
int n_colors;
float color_exp;
int t_stop;
int n_move;
float x_pos;
float y_pos;
float z_pos;
float x_foc;
float y_foc;
float z_foc;
int xsize;
int ysize;
int image_magnification;
float intensity;
float function_magnification;

// data dimensions (mainly queried from simulation 
// directory).

int n_rp;
int n_alpha;

int n_theta_plot;
int n_spec;
int n_kinetic;
int n_r;
int n_field;
int n_n;
int n_time;
int boundary_method;

int n_elements;

int i_part;
int i_alpha;
int i_r; 
int i_n;
int i_field;
int i_time;
int i_buffer;
int i_n_0;

float *f_real;
float *f_imag;
float *f1_real;
float *f1_imag;
float *f2_real;
float *f2_imag;
float *t;
float *theta_plot;
float ymin;
float ymax;
float r_center;
float ly;

int n0;
int dn;
int *n;

float *kr;

// Shape/profile parameters:

float *grid;

float *R0;
float *delta;
float *kappa;
float *r;
float *r_s;
float *q;

float R0_halo;
float delta_halo;
float kappa_halo;

// Color indices

float r_i;
float g_i;
float b_i;

// Angular working variables

float alpha;
float rp;
float d_alpha;
float d_r;
float d_theta_plot;

// Camera working variables

float x_cam;
float y_cam;
float eta;

// Motion variables

int i_move;
int i_frame;

// Constants

float pi = 3.1415926535;

string path;
string simdir;

