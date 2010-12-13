import graph;

include tgyro_find ;

scale(Linear,Linear);

file fin=input(root+"/profile.out").line();

real [] r;
real [] ti;
real [] te;

real ymax;

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,9);
   x = transpose(x);

   r  = x[0];
   ti = x[3];
   te = x[4];

   if (i == 0) {
      draw(graph(r,ti),black,"$T_i^{"+((string) i)+"}$");
      if (max(ti) > ymax) ymax = max(ti);
   }

   if (i == i_min) {
      draw(graph(r,ti),magenta,"$T_i^{"+((string) i)+"}$");
      if (max(ti) > ymax) ymax = max(ti);
   }
}

xlimits(xa,xb);
ylimits(0.0,ymax);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$T_i$ [keV]",LeftRight,RightTicks);

attach(legend(3,invisible),(point(S).x,truepoint(S).y),10S,UnFill);
