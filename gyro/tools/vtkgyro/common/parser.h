/*==========================================================
  module: read_input
 
  read all of input.dat
  ==========================================================*/

int read_vugyro(string file)
{
  int i;
  float *x;
  float *xs;

  ifstream in(file.c_str());

  if (!in) {
    cout << "Cannot open " << file << endl;
    return 1;
  }

  in >> n_r;
  in >> i;
  in >> i;
  in >> i;
  in >> i;
  in >> n_theta_plot;
  in >> n0;
  in >> n_n;
  in >> dn;
  in >> i_buffer;
  in >> i;
  in >> i;
  in >> n_field;
  in >> i;
  in >> n_kinetic;
  in >> n_spec;
  in >> i;
  in >> i;
  in >> i;
  in >> boundary_method;  

  x     = new float[n_r];
  xs    = new float[n_r*n_spec];
  r     = new float[n_r];
  r_s   = new float[n_r];
  q     = new float[n_r];
  delta = new float[n_r];
  kappa = new float[n_r];
  R0    = new float[n_r];
  kr    = new float[n_r];

  for (i=0; i<n_r; i++) in >> r[i]; 
  for (i=0; i<n_r; i++) in >> q[i];
  for (i=0; i<n_r; i++) in >> r_s[i];
  for (i=0; i<n_r; i++) in >> x[i];
  for (i=0; i<n_r*n_spec; i++) in >> xs[i];
  for (i=0; i<n_r*n_spec; i++) in >> xs[i];
  for (i=0; i<n_r*n_spec; i++) in >> xs[i];
  for (i=0; i<n_r*n_spec; i++) in >> xs[i];
  for (i=0; i<n_r; i++) in >> x[i];
  for (i=0; i<n_r; i++) {
    in >> R0[i];
    R0[i] = R0[i]*r_s[i]/2;
  }
  for (i=0; i<n_r; i++) in >> delta[i];
  for (i=0; i<n_r; i++) in >> x[i];
  for (i=0; i<n_r; i++) in >> kappa[i];

  in.close();

  return 0;
}

/*==========================================================
  module: read_progress
 
  Progress meter.
  ==========================================================*/

int read_progress(int i0,int n) 
{

  int i;

  i = i0+1;

  if (i == 0*n/20)  cout <<   "[                    ]" << flush;
  if (i == 1*n/20)  cout << "\r[>                   ]" << flush;
  if (i == 2*n/20)  cout << "\r[->                  ]" << flush;
  if (i == 3*n/20)  cout << "\r[-->                 ]" << flush;
  if (i == 4*n/20)  cout << "\r[--->                ]" << flush;
  if (i == 5*n/20)  cout << "\r[---->               ]" << flush;
  if (i == 6*n/20)  cout << "\r[----->              ]" << flush;
  if (i == 7*n/20)  cout << "\r[------>             ]" << flush;
  if (i == 8*n/20)  cout << "\r[------->            ]" << flush;
  if (i == 9*n/20)  cout << "\r[-------->           ]" << flush;
  if (i == 10*n/20) cout << "\r[--------->          ]" << flush;
  if (i == 11*n/20) cout << "\r[---------->         ]" << flush;
  if (i == 12*n/20) cout << "\r[----------->        ]" << flush;
  if (i == 13*n/20) cout << "\r[------------>       ]" << flush;
  if (i == 14*n/20) cout << "\r[------------->      ]" << flush;
  if (i == 15*n/20) cout << "\r[-------------->     ]" << flush;
  if (i == 16*n/20) cout << "\r[--------------->    ]" << flush;
  if (i == 17*n/20) cout << "\r[---------------->   ]" << flush;
  if (i == 18*n/20) cout << "\r[----------------->  ]" << flush;
  if (i == 19*n/20) cout << "\r[------------------> ]" << flush;
  if (i == n)       cout << "\r[------------------->]" << endl;

  return 0;

}

/*==========================================================
  module: read_time
 
  read t.out
  ==========================================================*/

int read_time(string file, int *n_elements)
{
  int dummy;
  int i;

  float fdummy;  


  ifstream in(file.c_str());

  if (!in) {
    cout << "Cannot open " << file << endl;
    return 1;
  }

  i = 0;
  while (!in.eof()) {
    in >> dummy >> fdummy;
    if (in.eof()) break;
    i++;
  }

  *n_elements = i;

  in.close();
  return 0;
}

