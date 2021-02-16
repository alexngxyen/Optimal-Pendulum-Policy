% ECE 594D Project - Simple Pendulum with Unlimited Voltage Applying Torque
% By: Alex Nguyen

clc; clear; close all;

% System Parameters
m = 0.5; % mass [kg]
b = 0.1; % damping coefficient [N/(m/s)] 
L = 1; % length [m]
g = 9.81; % gravity [m/s^2]

% % Solve ODE 
% tspan = 0:0.1:10; % time span [s]
% y0 = [pi/2; 0]; % initial conditions [rad;rad/s]
% 
% [t,y] = ode45(@(t,y) vargin(y,m,L,g,b,0),tspan,y0);

% Pendulum SS Matrices - Linearized about \pi
% EOM: \ddot{\theta} = (-b/m/L)\dot{\theta} - (g/L)sin(\theta) + (1/m/L^2)u
theta = pi; % angle [rad] 
A = [0 1; -g/L*cos(theta) -b/m/L];
B = [0; 1/m/L/L];

A_poles = eig(A);

% Check Controllability of A Matrix
n = rank(ctrb(A,B)); % see if full rank

% Pole Placement - Choose Eigenvalues
p = [-4 + 7*1i; -4 - 7*1i]; % imaginary parts (oscillations) - vary imaginary part
% p = [-5; -7]; % no imaginary parts (no overshoot)
K = place(A, B, p); % gain matrix

% Solve ODE - Pole Placement
tspan = 0:0.005:2; % time span [s]
y0 = [0; 0]; % initial conditions [rad;rad/s]
[t,y] = ode45(@(t,y) pendyn(y,m,L,g,b,-K*(y - [pi; 0])),tspan,y0);
u = -K*(y' - [pi; 0]); % torque input

% Print Results
fprintf('Assume: No limit on torque (i.e. unlimited voltage supply) \n')
fprintf('Torque Input Range: [%4.4f %4.4f]\n',min(u),max(u))

%% Animation
% Origin
O = [0 0];
axis(gca,'equal'); %aspect ratio of the plot
axis([-1.25 1.25 -1.25 1.25]); %limit of the plot

% Loop for animation

for i = 1:length(t)
    % Mass point
    P = L*[sin(y(i,1)) -cos(y(i,1))];
    
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