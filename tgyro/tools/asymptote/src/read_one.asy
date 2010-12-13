file fin=line(input(root+"/control.out"));
int n_r = fin;
int n_iter = fin;

int i_min=0;

file fin=line(input(root+"/gyrobohm.out"));

real [] r;
real [] chi_gb;
real [] q_gb;
real [] gamma_gb;

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,4);
   x = transpose(x);

   r        = x[0];
   chi_gb   = x[1];
   q_gb     = x[2];
   gamma_gb = x[3];
}

//----------------------------------
// flux_i.out

real [] r;
real [] gi_neo;
real [] gi_tur;
real [] qi_neo;
real [] qi_tur;

file fin=line(input(root+"/flux_i.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,5);
   x = transpose(x);

   r      = x[0];
   gi_neo = x[1];
   gi_tur = x[2];
   qi_neo = x[3];
   qi_tur = x[4];
}

//----------------------------------
// chi_i.out

real [] r;
real [] di_neo;
real [] di_tur;
real [] chii_neo;
real [] chii_tur;

file fin=line(input(root+"/chi_i.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,5);
   x = transpose(x);

   r      = x[0];
   di_neo = x[1];
   di_tur = x[2];
   chii_neo = x[3];
   chii_tur = x[4];
}

//----------------------------------
// flux_e.out

real [] r;
real [] ge_neo;
real [] ge_tur;
real [] qe_neo;
real [] qe_tur;

file fin=line(input(root+"/flux_e.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,5);
   x = transpose(x);

   r      = x[0];
   ge_neo = x[1];
   ge_tur = x[2];
   qe_neo = x[3];
   qe_tur = x[4];
}

//----------------------------------
// chi_e.out

real [] r;
real [] de_neo;
real [] de_tur;
real [] chie_neo;
real [] chie_tur;

file fin=line(input(root+"/chi_e.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,5);
   x = transpose(x);

   r      = x[0];
   de_neo = x[1];
   de_tur = x[2];
   chie_neo = x[3];
   chie_tur = x[4];
}

//----------------------------------
// profile.out

real [] r;
real [] ni;
real [] ne;
real [] ti;
real [] te;

file fin=line(input(root+"/profile.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,8);
   x = transpose(x);

   r  = x[0];
   ni = x[1];
   ne = x[2];
   ti = x[3];
   te = x[4];
}

//----------------------------------
// gradient.out

real [] r;
real [] alni;
real [] alne;
real [] alti;
real [] alte;

file fin=line(input(root+"/gradient.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,7);
   x = transpose(x);

   r  = x[0];
   alni = x[1];
   alne = x[2];
   alti = x[3];
   alte = x[4];
}

//----------------------------------
// nu_rho.out

real [] r;
real [] nui;
real [] nue;

file fin=line(input(root+"/nu_rho.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,1);
   real[][] x=dimension(fin,n_r,6);
   x = transpose(x);

   r  = x[0];
   nui = x[1];
   nue = x[2];
}

//----------------------------------
// geometry.out

real [] r;
real [] q;
real [] s;
real [] kappa;
real [] s_kappa;
real [] delta;
real [] s_delta;
real [] shift;
real [] rmaj;
real [] b_unit;

file fin=line(input(root+"/geometry.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,10);
   x = transpose(x);

   r       = x[0];
   q       = x[1];
   s       = x[2];
   kappa   = x[3];
   s_kappa = x[4];
   delta   = x[5];
   s_delta = x[6];
   shift   = x[7];
   rmaj    = x[8];
   b_unit  = x[9];
}

