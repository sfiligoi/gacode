import graph;

include tgyro_find ;

write("residual = ",z_min);

scale(Linear,Log);

draw(graph(r_res,res_ti+1e-20),black,"$R_i$");
draw(graph(r_res,res_te+1e-20),red,"$R_e$");
if (n_evolve > 2) {
   draw(graph(r_res,res_ne+1e-20),orange,"$R_{ne}$");
}

xlimits(xa,xb);

xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$R$",LeftRight,RightTicks);

attach(legend(2,invisible),(point(S).x,truepoint(S).y),10S,UnFill);
