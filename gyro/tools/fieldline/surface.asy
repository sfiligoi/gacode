import graph;

defaultpen(1.0);
pen dashed=linetype("8 4");
pen dotted=linetype("4 4");

size(400,400,IgnoreAspect);

scale(Linear,Linear);

write('Reading data in surface.out');

file fin=line(input("surface.out"));
int n_ic = fin;
real[][] a=dimension(fin,0,0);

a = transpose(a);

real[] xi    = a[0];
real[] alpha = a[1];

int n_point=xi.length;
int n0=(int) (n_point/n_ic);

for (int i=0; i<n_point; ++i) {

  int ic =(int) (i/n0);

  real x = ic*(1.0/n_ic);
  pen p  = black;

  dot((xi[i],alpha[i]),p+0.4);

}

xlimits(-0.5,0.5,Crop);
ylimits(-0.5,0.5,Crop);

xaxis("$\xi/(2\pi)$",BottomTop,LeftTicks);
yaxis("${\hat\alpha}/(2\pi)$",LeftRight,RightTicks);

//attach(legend(1,linelength=30,invisible),(50,105));

