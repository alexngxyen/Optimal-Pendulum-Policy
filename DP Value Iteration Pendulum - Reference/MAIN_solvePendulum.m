% ECE 594D - Value Iteration for Simple Pendulum 
% Partially By: Alex Nguyen 
%
% This script finds the optimal policy for an inverted pendulum controller,
% by modeling it as a MDP and the solving with value iteration. The goal is
% to maintain upright balance, while minimizing actuator torque and
% preventing falls (nearing the edge of the boundry).
%
%
% Note: This code was modified from a GitHub page. I deleted unnecessary lines 
% and restructured the script to be more readable. Also, I edited the 
% functions for the main script as well. Finally, I changed some of the 
% parameters for the optimal policy portion of the script then added an 
% animation at the end for observing the pendulum's dynamics as well as
% seeing its phase portrait and displacement over time.
%
% GitHub = https://github.com/MatthewPeterKelly/MDP_Pendulum#readm

clc; clear; close all;

%% Optimal Policy Set-Up

% Pendulum Parameters
pendulum.dyn.m = 1; % mass [kg]
pendulum.dyn.g = 9.81; % gravity [m/s^2]
pendulum.dyn.c = 0.1; % damping coefficient [N/m/s]
pendulum.dyn.l = 0.5; % rod length [m]
pendulum.dyn.timeConstant = sqrt(pendulum.dyn.l/pendulum.dyn.g);
pendulum.dyn.timeStep = 0.05*pendulum.dyn.timeConstant;

% Markov Decision Process (MDP) Grid
pendulum.grid.th.bnd = [-pi/4,pi/4];  % angle [rad]
pendulum.grid.th.nGrid = 45;
pendulum.grid.w.bnd = [-1,1];  % velocity [rad/s]
pendulum.grid.w.nGrid = 45;
gravityTorque = pendulum.dyn.m*pendulum.dyn.g*pendulum.dyn.l;
pendulum.grid.u.bnd = gravityTorque*[-1,1];  % torque [N*m]
pendulum.grid.u.nGrid = 25;

% Value Iteration Parameters MDP Solver
pendulum.opt.discount = 0.9;  % decay after one time constant
pendulum.opt.convergeTol = 1e-8; % convergence tolerance
pendulum.opt.maxIter = 1e4;
pendulum.opt.dispMod = 25;  

% Cost Function Parameters
pendulum.opt.cost.angle = 100;  % Quadratic
pendulum.opt.cost.rate = 5; % Quadratic
pendulum.opt.cost.torque = 1; % Quadratic
pendulum.opt.cost.boundry.cost = 10;   % Penalty at the boundry
pendulum.opt.cost.boundry.width = 0.1;  % Boundry width (fraction)


%% Build MDP System Model and Solve                   

[V, P, S, A] = MDP_Pend(pendulum);

%% Plot Results

soln.th = S(1,:);
soln.w = S(2,:);
soln.u = A(P);
soln.markerSize = 50*ones(size(soln.u));

figure(1); clf; colormap('jet')
scatter3(soln.th, soln.w, soln.u, soln.markerSize, soln.u, 'filled');
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')
zlabel('Torque (N*m)')
title('Optimal Policy')
view(2)
axis equal

figure(2); clf; colormap('jet')
scatter3(soln.th, soln.w, V.^(1/4), soln.markerSize, V.^(1/4), 'filled');
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')
zlabel('Value')
title('Value (Cost-to-Go) Function')
view(3)

%% Run Test Simulations
figure(3); clf; hold on;
test.nSim = 100;
test.duration = 2;
test.nTime = ceil(test.duration/pendulum.dyn.timeStep);
test.t = linspace(0,test.duration,test.nTime);
test.th = pendulum.grid.th.bnd(1) + diff(pendulum.grid.th.bnd)*rand(1,test.nSim); % intial conditions
test.w = pendulum.grid.w.bnd(1) + diff(pendulum.grid.w.bnd)*rand(1,test.nSim);
test.sys = @(t,z) (pendulumSystem(z,A,P,pendulum.grid,pendulum.dyn));
test.z = rk4(test.sys,test.t,[test.th;test.w]);

for i=1:test.nSim
    thTest = reshape(test.z(1,i,:),test.nTime,1);
    wTest = reshape(test.z(2,i,:),test.nTime,1);
    u = pendulumController([thTest,wTest]',A,P,pendulum.grid);
    plot3(thTest, wTest, u, 'k-');
    plot3(thTest(1), wTest(1), u(1),'r.', 'MarkerSize',15);
    plot3(thTest(end), wTest(end), u(end), 'b.', 'MarkerSize',15);
end
xlabel('$\theta$ [rad]','interpreter','latex')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')
zlabel('Torque (N*m)')
title(sprintf('Simulation (Trials = %d)',test.nSim))
view(2)

%% Animation - Simple Pendulum
x = 20; % pendulum result from inital condition index 
fps = 3;           % Frames per second
interval = [0, test.duration]; % Time span
ivp = [test.th(x); test.w(x); pendulum.dyn.g; pendulum.dyn.m;  ...
    pendulum.dyn.l ; pendulum.dyn.c];     % Initial value's for the problem

% Position, Angular velocity, and Controller
t = test.t; % time
phi = reshape(test.z(1,x,:),test.nTime,1);
dtphi = reshape(test.z(2,x,:),test.nTime,1);
u = pendulumController([thTest,wTest]',A,P,pendulum.grid);
L = ivp(5); 

% To set the Range os Phase plane, time vs. depl plots
minu = 1.1*min(phi) ; maxu = 1.1*max(phi);
minv = 1.1*min(dtphi) ; maxv = 1.1*max(dtphi);

fh = figure ;
set(fh,'name','The Simple Pendulum','numbertitle','off','color', 'w','menubar','none') ;
stop = uicontrol('style','toggle','string','stop','background','w');

% Plot for Pendulum
subplot(121);
h = plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',1.5,'Color','b');
title('Simple Pendulum Animation','Color','b');
range = 1.1*L;
axis([-range range -range range]);
axis square;
set(gca,'XTickLabelMode', 'manual', 'XTickLabel', [],'YTickLabelMode', .....
    'manual', 'YTickLabel', [],'nextplot','replacechildren');

% Plot for Phase plane
subplot(222) ;
h1 = plot(ivp(1),ivp(2),'LineWidth',1,'Color','m') ;
axis([minu maxu minv maxv]) ;
xlabel('\phi') ;ylabel('$\dot{\phi}$','interpreter','latex') ;
set(get(gca,'YLabel'),'Rotation',0.0)
set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot','Color','m')

% Plot for time Vs. displacement 
subplot(224) ;
h2 = plot(t(1),ivp(1),'LineWidth',1,'Color','r') ;
axis([0 test.duration minu maxu]) ;
xlabel('t') ;ylabel('\phi') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Displacement Plot','Color','r');

% Animation starts
for i=1:length(phi)-1
    
 % Animation Plot
    if (ishandle(h)==1)
        Xcoord = [0, -L*sin(phi(i))];
        Ycoord = [0, L*cos(phi(i))];
        set(h,'XData',Xcoord,'YData',Ycoord);
        if get(stop,'value')==0
            drawnow;
        elseif get(stop,'value')==1
            break;
        end
        
        % Phase Plane Plot
        if (ishandle(h1)==1)
            PP(i,:) = [phi(i) dtphi(i)];
            set(h1,'XData',PP(:,1),'YData',PP(:,2));
            drawnow;    
            vasu = length(PP(:,1)) ;
      
        % Time Vs. displacement Plot  
            if (ishandle(h2)==1)
                DEPL(i,:) = [t(i) phi(i)] ;
                set(h2,'Xdata',DEPL(:,1),'YData',DEPL(:,2)) ;
                drawnow ;    
            end    
        end  
    end
    F(i) = getframe(fh) ;          
end

% Close the Figure window
set(stop,'style','pushbutton','string','close','callback','close(gcf)');
