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
real [] gi2_neo;
real [] gi2_tur;
real [] qi2_neo;
real [] qi2_tur;

file fin=line(input(root+"/flux_i2.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,5);
   x = transpose(x);

   r      = x[0];
   gi2_neo = x[1];
   gi2_tur = x[2];
   qi2_neo = x[3];
   qi2_tur = x[4];
}

//----------------------------------
// profile_2.out

real [] r;
real [] ni2;
real [] alni2;

file fin=line(input(root+"/profile_2.out"));

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,3);
   x = transpose(x);

   r   = x[0];
   ni2 = x[1];
   alni2 = x[2];
}

