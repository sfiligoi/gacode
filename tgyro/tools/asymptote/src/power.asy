// Transport powers

import graph;

include tgyro_find ;

scale(Linear,Linear);

draw(graph(r,pwr_alpha),blue,"$P_\alpha$");
draw(graph(r,pwr_brem),purple,"$P_{\rm brem}$");
draw(graph(r,pwr_exch),orange,"$P_{\rm exch}$");
draw(graph(r,pwr_iaux),black,"$P_{\rm i,aux}$");
draw(graph(r,pwr_eaux),red,"$P_{\rm e,aux}$");
draw(graph(r,pwr_i),black+dashed,"$P_{\rm i}$");
draw(graph(r,pwr_e),red+dashed,"$P_{\rm e}$");

xlimits(xa,xb);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$P$ [MW]",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);

