// Smooth temperature profiles

import graph;

include tgyro_find ;

scale(Linear,Linear);

// Use approximately 256 points
int n=(int) (256/n_r);
real [] rf = new real[(n_r-1)*n+1];
real [] zif = new real[(n_r-1)*n+1];
real [] zef = new real[(n_r-1)*n+1];
real [] tif = new real[(n_r-1)*n+1];
real [] tef = new real[(n_r-1)*n+1];

real ymax;
real ymax_e;
real ymax_i;

int j=0;
for (int i=0; i<n_r-1; ++i) {
   for (int m=0; m<n; ++m) {
      z = m/(1.0*n);
      rf[j]  = r[i]*(1-z)+r[i+1]*z;
      zif[j] = alti[i]*(1-z)+alti[i+1]*z;
      zef[j] = alte[i]*(1-z)+alte[i+1]*z;
      j = j+1;
   } 
}
rf[j]  = r[n_r-1];
zif[j] = alti[n_r-1];
zef[j] = alte[n_r-1];

tif[j] = ti[n_r-1];
tef[j] = te[n_r-1];

for (int i=j; i>0; --i) {
   tif[i-1] = tif[i]*exp(0.5*(rf[i]-rf[i-1])*(zif[i]+zif[i-1]));
   tef[i-1] = tef[i]*exp(0.5*(rf[i]-rf[i-1])*(zef[i]+zef[i-1]));
}

draw(graph(rf,tif),black,"$T_i$");
draw(graph(rf,tef),red,"$T_e$");

for (int i=0; i<n_r; ++i) {
   dot((r[i],te[i]),black+4);
}

for (int i=0; i<n_r; ++i) {
   dot((r[i],ti[i]),black+4);
}

ymax_e = max(tef);
ymax_i = max(tif);
ymax = max(ymax_e,ymax_i);
ymax = 1.10*ymax;

xlimits(xa,xb);
ylimits(0.0,ymax);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$T$ [keV]",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);

