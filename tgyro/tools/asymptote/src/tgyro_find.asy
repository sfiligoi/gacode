//=================================================================
// tgyro_find.asy
//
// PURPOSE: Read all profile data from vugyro at a given iteration. 
//
// INPUT:  string root   (TGYRO simdir containing data) 
//         int    i_grab (iteration to grab; -1 if minimum)
//=================================================================

//------------------------------
// Plot size defaults

defaultpen(1.5);
size(300,250,IgnoreAspect);

real xa=0.0;
real xb=1.0;
//-------------------------------

//======================================================
// First need to read dimensions from control.out
//
file fin=input("PATH");
string root=fin;
int i_grab=fin;

file fin=input(root+"/control.out").line();
int n_r = fin;
int n_evolve = fin;
int n_iter = fin;

if (i_grab >= 0) n_iter = i_grab;
//======================================================

//======================================================
// Now determine iteration which gives minimum residual
//
file fin=input(root+"/residual.out").line();

real z;
int i_min;
real z_min;

real[] r_res;
real[] res_ti;
real[] relax_ti;
real[] res_te;
real[] relax_te;
real[] res_ne;
real[] relax_ne;

z_min = 1000.0;

for (int i=0; i<=n_iter; ++i) {
   string[] t=fin.dimension(1);
   real[][] x=fin.dimension(n_r-1,2*n_evolve+1);
   x = transpose(x);

   r_res    = x[0];
   res_ti   = x[1];
   relax_ti = x[2];
   res_te   = x[3];
   relax_te = x[4];

   z = sum(res_ti+res_te)/(2*n_r-2);

   // Add density evolution residual if present
   if (n_evolve > 2) {
      res_ne   = x[5];
      relax_ne = x[6];
      z = z+sum(res_ne)/(2*n_r-2);
   }

   if (i > 0 && z < z_min) {
      z_min = z;
      i_min = i;
   }
 
   if (i == i_grab) {
      z_min = z;
      i_min = i;
   }

}
   
//======================================================

//======================================================
// gyrobohm.out

file fin=input(root+"/gyrobohm.out").line();

real [] r;
real [] chi_gb;
real [] q_gb;
real [] gamma_gb;
real [] mpi_gb;

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,6);
   x = transpose(x);

   r        = x[0];
   chi_gb   = x[1];
   q_gb     = x[2];
   gamma_gb = x[3];
   mpi_gb   = x[4];
}
//======================================================

//======================================================
// flux_i.out

real [] r;
real [] gi_neo;
real [] gi_tur;
real [] qi_neo;
real [] qi_tur;

file fin=input(root+"/flux_i.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,5);
   x = transpose(x);

   r      = x[0];
   gi_neo = x[1];
   gi_tur = x[2];
   qi_neo = x[3];
   qi_tur = x[4];
}
//======================================================

//======================================================
// chi_i.out

real [] r;
real [] di_neo;
real [] di_tur;
real [] chii_neo;
real [] chii_tur;

file fin=input(root+"/chi_i.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,5);
   x = transpose(x);

   r      = x[0];
   di_neo = x[1];
   di_tur = x[2];
   chii_neo = x[3];
   chii_tur = x[4];
}
//======================================================

//======================================================
// flux_e.out

real [] r;
real [] ge_neo;
real [] ge_tur;
real [] qe_neo;
real [] qe_tur;

file fin=input(root+"/flux_e.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,5);
   x = transpose(x);

   r      = x[0];
   ge_neo = x[1];
   ge_tur = x[2];
   qe_neo = x[3];
   qe_tur = x[4];
}
//======================================================

//======================================================
// flux_target.out

real [] r;
real [] qi;
real [] qi_t;
real [] qe;
real [] qe_t;
real [] ge;
real [] ge_t;

file fin=input(root+"/flux_target.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,7);
   x = transpose(x);

   r    = x[0];
   qi   = x[1];
   qi_t = x[2];
   qe   = x[3];
   qe_t = x[4];
   ge   = x[5];
   ge_t = x[6];
}
//======================================================

//======================================================
// chi_e.out

real [] r;
real [] de_neo;
real [] de_tur;
real [] chie_neo;
real [] chie_tur;

file fin=input(root+"/chi_e.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,5);
   x = transpose(x);

   r      = x[0];
   de_neo = x[1];
   de_tur = x[2];
   chie_neo = x[3];
   chie_tur = x[4];
}
//======================================================

//======================================================
// profile.out

real [] r;
real [] ni;
real [] ne;
real [] ti;
real [] te;

file fin=input(root+"/profile.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,8);
   x = transpose(x);

   r  = x[0];
   ni = x[1];
   ne = x[2];
   ti = x[3];
   te = x[4];
}
//======================================================

//======================================================
// gradient.out

real [] r;
real [] alni;
real [] alne;
real [] alti;
real [] alte;

file fin=input(root+"/gradient.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,8);
   x = transpose(x);

   r  = x[0];
   alni = x[1];
   alne = x[2];
   alti = x[3];
   alte = x[4];
}
//======================================================

//======================================================
// nu_rho.out

real [] r;
real [] nui;
real [] nue;

file fin=input(root+"/nu_rho.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(1);
   real[][] x=fin.dimension(n_r,6);
   x = transpose(x);

   r   = x[0];
   nui = x[1];
   nue = x[2];
}
//======================================================

//======================================================
// geometry.out

real [] r;
real [] rho;
real [] q;
real [] s;
real [] kappa;
real [] s_kappa;
real [] delta;
real [] s_delta;
real [] shift;
real [] rmaj;
real [] b_unit;

file fin=input(root+"/geometry.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,11);
   x = transpose(x);

   r       = x[0];
   rho     = x[1];
   q       = x[2];
   s       = x[3];
   kappa   = x[4];
   s_kappa = x[5];
   delta   = x[6];
   s_delta = x[7];
   shift   = x[8];
   rmaj    = x[9];
   b_unit  = x[10];
}
//======================================================

//======================================================
// power.out

real [] r;
real [] pwr_alpha;
real [] pwr_brem;
real [] pwr_exch;
real [] pwr_iaux;
real [] pwr_eaux;
real [] pwr_i;
real [] pwr_e;

file fin=input(root+"/power.out").line();

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,8);
   x = transpose(x);

   r         = x[0];
   pwr_alpha = x[1];
   pwr_brem  = x[2];
   pwr_exch  = x[3];
   pwr_iaux  = x[4];
   pwr_eaux  = x[5];
   pwr_i     = x[6];
   pwr_e     = x[7];
}
//======================================================

