// Simulation flux versus target flux (GB units)

import graph;

include tgyro_find ;

scale(Linear,Linear);

draw(graph(r,qi),black,"$Q_i$");
draw(graph(r,qi_t),black+dotted,"$Q^T_i$");

draw(graph(r,qe),red,"$Q_e$");
draw(graph(r,qe_t),red+dotted,"$Q^T_e$");

xlimits(xa,xb);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$Q/Q_{GB}$",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);
