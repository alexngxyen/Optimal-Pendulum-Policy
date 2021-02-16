function dy = pendyn(y,m,L,g,b,u)
% Equation of motion for nonlinear simple pendulum dynamics which includes 
% damping

Sy = sin(y(1)); % unlinearized variable

dy(1,1) = y(2);
dy(2,1) = -b/m/L*y(2) - g/L*Sy + 1/m/L^2*u;
end