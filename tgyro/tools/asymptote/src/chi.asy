// Neoclassical fraction of total diffusivity

import graph;

include tgyro_find ;

scale(Linear,Log);

real [] chii_f;
real [] chie_f;

chii_neo[0]=1.0;
chie_neo[0]=1.0;

chii_f = chii_neo/(chii_neo+chii_tur);
chie_f = chie_neo/(chie_neo+chie_tur);

draw(graph(r[1:n_r-1],chii_f[1:n_r-1]),black,
 "$\chi^{\rm neo}_i/\chi^{\rm tot}_i$");
draw(graph(r[1:n_r-1],chie_f[1:n_r-1]),red,
 "$\chi^{\rm neo}_e/\chi^{\rm tot}_e$");

xlimits(xa,xb);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);
