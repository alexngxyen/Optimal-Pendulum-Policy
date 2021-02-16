% ECE 594D Project - Animating Simple Pendulum Dynamics
% By: Alex Nguyen

clc; clear; close all;

% Pendulum System Parameters
m = 1; b = 0.1; L = .5; g = 9.81;

% Initial Conditions
xset = [-pi/2 pi/2]; yset = [-1 1]; % limit initial conditions
x0 = [xset(1) + diff(xset)*rand; yset(1) + diff(yset)*rand]; % random initial conditions [rad; rad/s]

% solve ODE 
dt = sqrt(g/L)*0.005;
tspan = 0:dt:20;
[t,x] = ode45(@(t,x) pendyn(x,m,L,g,b,0),tspan,x0);

%% Animation
% Origin
O = [0 0];
axis(gca,'equal'); %aspect ratio of the plot
axis([-0.7 0.7 -0.7 0.2]); %limit of the plot

% Loop for animation
for i = 1:length(t)
    % Mass point
    P = L*[sin(x(i,1)) -cos(x(i,1))];
    
    % Circle in origin
    O_circ = viscircles(O,0.01);
    
    % Pendulum
    pend = line([O(1) P(1)],[O(2) P(2)]);
    
    % Ball
    ball = viscircles(P,0.05);
    
    % Time interval to update the plot
    pause(0.001);
    
    % Delete previous objects if it is not the final loop
    if i < length(t)
        delete(pend);
        delete(ball);
        delete(O_circ);
    end
end
