// Simulation flux versus target flux (physical units)

import graph;

include tgyro_find ;

scale(Linear,Linear);

draw(graph(r,qi*q_gb),black,"$Q_i$");
draw(graph(r,qi_t*q_gb),black+dotted,"$Q^T_i$");

draw(graph(r,qe*q_gb),red,"$Q_e$");
draw(graph(r,qe_t*q_gb),red+dotted,"$Q^T_e$");

xlimits(xa,xb);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$Q$ (MW/m$^2$)",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);
