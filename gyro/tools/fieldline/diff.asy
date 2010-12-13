import graph;

defaultpen(1.8);
pen dashed=linetype("8 4");
pen dotted=linetype("4 4");

size(300,300,IgnoreAspect);

scale(Linear,Linear);

file fin=line(input("surface.out"));
int n_ic = fin;
real[][] a=dimension(fin,0,0);

a = transpose(a);

real[] xi    = a[0];
real[] alpha = a[1];
real[] diff  = a[2];

int n_point=xi.length;
int n0=(int) (n_point/n_ic);

real[] t = new real[n0];
real[] d = new real[n0];

int ix=0;

for (int i=0; i<n_ic; ++i) {
  
   for (int j=0; j<n0; ++j) {

    t[j] = j;
    d[j] = diff[ix];

    ix = ix+1;

   }

  real x = i*(1.0/n_ic);
  pen p = rgb(0.5*x^2,0.5-0.5*x,0.5*x-0.5);
 
  draw(graph(t,d),p+0.3);

}

//xlimits(-0.5,0.5,Crop);
//ylimits(-0.5,0.5,Crop);

xaxis("$n_{\rm turn}$",BottomTop,LeftTicks);
yaxis("$\langle (\xi-\xi_0)^2 \rangle/n_{\rm turn}$",LeftRight,RightTicks);

//attach(legend(1,linelength=30,invisible),(50,105));

