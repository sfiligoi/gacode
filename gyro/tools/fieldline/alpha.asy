import graph;

defaultpen(0.5);
//pen dashed=linetype("8 4");
//pen dotted=linetype("4 4");

size(600,150,IgnoreAspect);

scale(Linear,Linear);

real pi = 3.1415926535;

file fin=line(input("alpha.out"));
int n=fin;
real[][] a=dimension(fin,n,3);
a = transpose(a);

real[] alpha     = a[0];
real[] xi_dot    = a[1];
real[] alpha_dot = a[2];

draw(graph(alpha,xi_dot),black,"$d\xi/dt$");
draw(graph(alpha,alpha_dot),blue,"$d{\hat\alpha}/dt$");

xlimits(min(alpha),max(alpha),Crop);

xaxis("${\hat\alpha}/(2\pi)$",BottomTop,LeftTicks);
yaxis("",LeftRight,RightTicks);

attach(legend(3,invisible),(point(S).x,truepoint(S).y),10S,UnFill);

