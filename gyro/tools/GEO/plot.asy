import graph;
defaultpen(1.0);

//----------------------------------------
file fin=line(input("geov.out"));

real[][] a=dimension(fin,0,0);
a=transpose(a);

real[] t=a[0], u1=a[1], u2=a[2], u3=a[3], u4=a[4];
//----------------------------------------

//======================================================
size(300,200,IgnoreAspect);
draw(graph(t/pi,u1),"gsin");
draw(graph(t/pi,u2),blue,"usin");
draw(graph(t/pi,sin(t)),dashed,"$\sin\theta$");

real ymin=-1.1;
real ymax=1.1;

pair p = (1,(ymin+ymax)/2);
attach(legend(1,invisible),(1.2,0));

xlimits(-1,1);
ylimits(ymin,ymax);

xaxis("$\theta/\pi$",BottomTop,LeftTicks);
yaxis("${\rm sin}(\theta)$",LeftRight,RightTicks);
shipout("sin.eps"); 
erase();
//======================================================
//======================================================
size(300,200,IgnoreAspect);
draw(graph(t/pi,u3),"gcos");
draw(graph(t/pi,u4),blue,"ucos");
draw(graph(t/pi,cos(t)),dashed,"$\cos\theta$");

real ymin=-1.1;
real ymax=1.1;

pair p = (1,(ymin+ymax)/2);
attach(legend(1,invisible),(1.2,0));

xlimits(-1,1);
ylimits(ymin,ymax);

xaxis("$\theta/\pi$",BottomTop,LeftTicks);
yaxis("${\rm cos}(\theta)$",LeftRight,RightTicks);
shipout("cos.eps"); 
erase();
//======================================================
