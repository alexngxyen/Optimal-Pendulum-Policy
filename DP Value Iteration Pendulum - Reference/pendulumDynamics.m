function ddth = pendulumDynamics(th,dth,u,p)
ddth = (p.g/p.l)*sin(th) + (-p.c*dth + u)/(p.m*p.l*p.l);
end