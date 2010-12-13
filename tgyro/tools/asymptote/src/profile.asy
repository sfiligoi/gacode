import find;
import graph;

// BEGIN OPTIONS
real x_min=0.0;
real x_max=1.0;
real y_min=0.0;
real y_max=50.0;
// END OPTIONS

write("Generating profile.eps.");

defaultpen(1.5);

size(350,250,IgnoreAspect);
scale(Linear,Linear);

file fin=line(input(root+"/control.out"));
int n_r = fin;
int n_iter = fin;

file fin=line(input(root+"/profile.out"));

real [] r;
real [] ti;
real [] te;

for (int i=0; i<=i_min; ++i) {
   string[] t=dimension(fin,2);
   real[][] x=dimension(fin,n_r,8);
   x = transpose(x);

   r  = x[0];
   ti = x[3];
   te = x[4];
}

draw(graph(r,ti),black,"$T_i$");
draw(graph(r,te),red,"$T_e$");

xlimits(0.0,1.0);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$T$ [keV]",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);

