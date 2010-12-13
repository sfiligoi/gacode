import graph;

include tgyro_find ;

scale(Linear,Linear);

file fin=input(root+"/profile.out").line();

real [] r;
real [] ni;
real [] ne;

real ymax;

for (int i=0; i<=i_min; ++i) {
   string[] t=fin.dimension(2);
   real[][] x=fin.dimension(n_r,9);
   x = transpose(x);

   r  = x[0];
   ni = x[1];
   ne = x[2];

   if (i == 0) {
      draw(graph(r,ne),black,"$n_e^{"+((string) i)+"}$");
      if (max(ne) > ymax) ymax = max(ne);
   }

   if (i == i_min) {
      draw(graph(r,ne),magenta,"$n_e^{"+((string) i)+"}$");
      if (max(ne) > ymax) ymax = max(ne);
   }
}

xlimits(xa,xb);
ylimits(0.0,ymax);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$n_e$ $[1/cm^3]$",LeftRight,RightTicks);

attach(legend(3,invisible),(point(S).x,truepoint(S).y),10S,UnFill);
