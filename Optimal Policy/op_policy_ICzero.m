% Optimal Policy for Simple Pendulum with Ball Mass - ICs = [0; 0]
% By: Alex Nguyen

clc; clear; close all;

% PENDULUM PARAMETERS
m = 0.5; % mass [kg]
b = 0.1; % damping coefficient [N/(m/s)] 
L = 1; % length [m]
g = 9.81; % gravity [m/s^2]

% SETTING UP DISCRETIZED MESH
u = -4:4; % action space
xset = linspace(-pi,pi,65); % \theta (x pos) 
yset = linspace(-3,3,65); % \dot{\theta} (y pos)
[x,y] = meshgrid(xset,yset);
Xcol = x(:); % \theta (position) values [rad]
Ycol = y(:); % \dot{\theta} (velocity) values [rad/s]
dt = 0.1; % time step [s]
N = [1:length(Xcol)]'; % # of discretized states

% COST FUNCTION
xgoal = 0; % want \theta to reach zero!! 
ygoal = 0; % want \dot{\theta} to reach zero!!
Vbig = 10000; % large initial penalty at each state
gamma = .95; % discount factor - try to vary this parameter
V = 0*Xcol + Vbig; % cost function value
fi = find((Xcol == xgoal) & (Ycol == ygoal));
V(fi) = 0; % GOAL has zero value/penalty
Uopt = 0*Xcol; % Optimal action, for each state

% BUILDING TRANSITION MATRIX
xnew = zeros(length(N),length(u)); ynew = xnew; % preallocation
Tmat = zeros(length(u),length(N)); 
for ni = 1:length(Xcol)
    for ui = 1:length(u)
        % \theta values
        xnew(ni,ui) = Xcol(ni) + dt*Ycol(ni); % new x-position
        
        if xnew(ni,ui) > xset(end) % "wrap" states to keep in xset range
            xnew(ni,ui) = xnew(ni,ui) - diff([xset(1) xset(end)]);
        elseif xnew(ni,ui) < xset(1)
            xnew(ni,ui) = xnew(ni,ui) + diff([xset(1) xset(end)]);
        end
        
        dx_temp = abs(Xcol - xnew(ni,ui));
        [~,Ix] = min(dx_temp);
        
        % \dot{\theta} values
        ynew(ni,ui) = Ycol(ni) + dt*(-b/m/L*Ycol(ni) - g/L*sin(Xcol(ni)) + ...
            1/m/L^2*u(ui)); % new y-position
        
        dy_temp = abs(Ycol - ynew(ni,ui));
        [~,Iy] = min(dy_temp);

        % Build "Transition Matrix"
        fi = find(((Xcol == Xcol(Ix)) .* (Ycol == Ycol(Iy))) == 1);
        Tmat(ui,ni) = fi;
        
        % ignore barycentric interpolation for now, using nearest neighbor
        % method
    end
end

% Then, keep iterating to find the "best possible option" to perform:
for n=1:100
    Vnext = 0*V;
    for xi = 1:length(Xcol)
        [vi,id] = min(V(Tmat(:,xi)));
        Vnext(xi) = gamma*vi + dt*(Xcol(xi) ~= xgoal || Ycol(xi) ~= ygoal);
        Uopt(xi) = id;
    end
    if abs(sum(V-Vnext)) < 1e-10
        fprintf('n = %d. No change in value!\n',n)
        break
    end
    V = Vnext;
end

% % Cost Functions
% Qf = 100; Q = 10; R = 0; % scalar weighting factors
% gn = @(th,dth) Qf*(abs(th) - pi).^2 + dth.^2; 
% gk = @(th,dth,u) Q*(abs(th) - pi).^2 + dth.^2 + R*u.^2;
% 
% % Find Optimal Policy
% a = 80; % optimal polciy iteration
% for n = a:-1:1
%     Vnext = 0*V;
%     for xi = 1:length(Xcol)
%         for ui = 1:length(ui)
%             if n == a
%                 Vnext(xi) = gn(xnew(xi,ui),ynew(xi,ui));
%                 Vn = Vnext;
%             else
%                 V = 0*V;
%                 V(xi) = gk(xnew(xi,ui),ynew(xi,ui),ui) + Vnext(xi);
%                 if V(xi) < Vnext(xi)
%                     Vnext(xi) = V(xi);
%                     Uopt(xi) = ui;
%                 end
%             end
%         end
%     end
% end  

fprintf('The value, V, gives time-to-go, until x = xgoal. It has converged.\n')    

% Visualization
Vfill = 0*x;
Ufill = Vfill;
Vfill(:) = V(:);
Ufill(:) = Uopt(:);

figure(1); clf
surf(x,y,Vfill)
view(0,90)
caxis([min(V) max(V)])
title('Cost-to-Go Function')
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')

figure(2); clf
surf(x,y,Ufill)
view(0,90)
caxis([0 11])
title('Optimal Policy')
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')

% figure(1); clf; colormap('jet')
% scatter3(Xcol, Ycol, Uopt,'filled');
% xlabel('$\theta$ [rad]','interpreter','latex')
% ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')
% zlabel('Torque (N*m)')
% title('Optimal Policy')
% view(2)
% axis equal
% 
% figure(2); clf; colormap('jet')
% scatter3(Xcol, Ycol, V.^(1/4),'filled');
% xlabel('$\theta$ [rad]','interpreter','latex')
% ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')
% zlabel('Value')
% title('Value (Cost-to-Go) Function')
% view(3)