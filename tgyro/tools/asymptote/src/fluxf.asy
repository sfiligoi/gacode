import graph;

include tgyro_find ;

scale(Linear,Linear);

draw(graph(r,ge*gamma_gb),red,"$F_e$");
draw(graph(r,ge_t*gamma_gb),red+dotted,"$F^T_e$");

xlimits(xa,xb);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$F_e$ (MW/keV/m$^2$)",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);
