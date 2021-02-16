% Optimal Policy for Simple Pendulum with Ball Mass
% By: Alex Nguyen

clc; clear; close all;

% PENDULUM PARAMETERS
m = 0.5; % mass [kg]
b = 0.1; % damping coefficient [N/(m/s)] 
L = 1; % length [m]
g = 9.81; % gravity [m/s^2]

% SETTING UP DISCRETIZED MESH
u = -5:5; % action space
xset = linspace(-pi,pi,45); % \theta (x pos) 
yset = linspace(-4,4,45); % \dot{\theta} (y pos)
[x,y] = meshgrid(xset,yset);
Xcol = x(:); % \theta (position) values [rad]
Ycol = y(:); % \dot{\theta} (velocity) values [rad/s]
dt = 0.1; % time step [s]
N = [1:length(Xcol)]'; % # of discretized states

% COST FUNCTION
xgoal = pi; % want \theta to reach \pi!!
ygoal = 0; % want \dot{\theta} to have zero velocity!!
Vbig = 500; % large initial penalty at each state
gamma = .95; % discount factor - try to vary this parameter
V = 0*Xcol + 0*Ycol + Vbig; % cost function value
fi = find((Xcol == xgoal) & (Ycol == ygoal)); 
V(fi) = 0; % GOAL has zero value/penalty
Uopt = 0*Xcol; % Optimal action, for each state

% BUILDING TRANSITION MATRIX
xnew = zeros(length(N),length(u)); ynew = xnew; % preallocation
Tmat = zeros(length(u),length(N)); 
for i = 1:length(N)
    for k = 1:length(u)
        % \theta values
        xnew(i,k) = Xcol(i) + dt*Ycol(i); % new x-position
        
        if xnew(i,k) > xset(end) % "wrap" states to keep in xset range
            xnew(i,k) = xnew(i,k) - diff([xset(1) xset(end)]);
        elseif xnew(i,k) < xset(1)
            xnew(i,k) = xnew(i,k) + diff([xset(1) xset(end)]);
        end
         
        [~,Ix] = min(abs(Xcol - xnew(i,k)));
        
        % \dot{\theta} values
        ynew(i,k) = Ycol(i) + dt*(-b/m/L*Ycol(i) - g/L*sin(Xcol(i)) + ...
            1/m/L^2*u(k)); % new y-position
        
        [~,Iy] = min(abs(Ycol - ynew(i,k)));
        
%         fi = find(((Xcol == Xcol(Ix)) .* (Ycol == Ycol(Iy))) == 1);
        fi = find((Xcol == Xcol(Ix)) & (Ycol == Ycol(Iy)));
        Tmat(k,i) = fi;
        
        % ignore barycentric interpolation for now, using nearest neighbor
        % method
    end
end

Tmat; % "Transition Matrix"

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

fprintf('The value, V, gives time-to-go, until x = xgoal. It has converged.\n')    

% Visualization

Vfill = 0*x;
Ufill = Vfill;
Vfill(:) = V(:);
Ufill(:) = Uopt(:);

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

%% Extra Code - Zoom Call Example
% xnew = zeros(length(N),length(u)); ynew = xnew; Tmat = xnew; % preallocation
% for i = 1:length(N)
%     for k = 1:length(u)
%         % \THETA VALUES
%         xnew(i,k) = Xcol(i) + dt*Ycol(i); % new x-position
%         [~,Ix] = min(abs(Xcol - xnew(i,k)));
%         
%         % \DOT{\THETA} VALUES
%         ynew(i,k) = Ycol(i) + dt*u(k); % new y-position
%         [~,Iy] = min(abs(Ycol - ynew(i,k)));
%         
%         fi = find((Xcol == Xcol(Ix)) .* (Ycol == Ycol(Iy)));
%         Tmat(i,k) = fi;
%         
%         % ignore barycentric interpolation for now, using nearest neighbor
%         % method
%     end
% end
