% ECE 594D Project - Torque-Limited Simple Pendulum Switching Controllers
% By: Alex Nguyen

clc; clear; close all;

%% Simple Pendulum Parameters
% System Parameters
m = 0.5; % mass [kg]
b = 0.1; % damping coefficient [N/(m/s)] 
L = 1; % length [m]
g = 9.81; % gravity [m/s^2]
tau = sqrt(g/L); % time constant [s]
dt = 0.001*tau; % time step [s]
T = m*g*L; %torque required [N*m]
ulim = sat(1)*T*[-1 1]; % saturated torque input range around [-4,4]

%% Energy-Shaping Control - Swing Up Portion
% E_tilde = 0.5*m*L^2*y(2)^2 - m*g*L*(1+cos(y(1)));
% u = -k*y(2)*E_tilde + b*y(2);
% u = -k*y(2)*(0.5*m*L^2*y(2)^2 - m*g*L*(1+cos(y(1)))) + b*y(2);

% Solve ODE 
tspan = 0:dt:20; % time span [s]
xset = [-pi/2 pi/2]; yset = [-0.5 0.5]; % limit initial conditions
% y1o = [xset(1) + diff(xset)*rand; yset(1) + diff(yset)*rand]; % random IC with a velocity value [rad; rad/s] 
y1o = [xset(1) + diff(xset)*rand; 0]; % random IC without a velocity value [rad; rad/s] 
k = 0.1; % saturates the torque input

% Check if Initial Conditions for Energy Control
if y1o(:) == 0 
    b = 0.5; % input to pendulum mass at IC = [0; 0]
    [t1,y1] = ode45(@(t,y) pendyn(y,m,L,g,b,-k*(y(2))*(0.5*m*(L*y(2))^2 - ...
        (m*g*L)*(1+cos(y(1)))) + b*(y(2))) + b,tspan,y1o);
else
    [t1,y1] = ode45(@(t,y) pendyn(y,m,L,g,b,-k*y(2)*(0.5*m*(L*y(2))^2 - ...
        (m*g*L)*(1+cos(y(1)))) + b*y(2)),tspan,y1o);
end

% Check Saturation - Need to be less than T (4.9 N*m)
u1 = -k*y1(:,2).*(0.5*m*(L*y1(:,2)).^2 - (m*g*L)*(1+cos(y1(:,1)))) + b*y1(:,2);
m1 = max(u1); m2 = min(u1); % actual min & max torque output

% Print Results
fprintf(['Simple Pendulum Input Torque Range: \n\tUnsaturated: [%4.2f %4.2f];' ...
    ' Saturated: [%4.2f %4.2f]\n\n'],-T,T,ulim(1),ulim(2))
fprintf('Energy Shape Control Torque Values: \n\t Max Torque = %4.2f; Min Torque = %4.2f\n\n',m1,m2)

for i = 1:length(y1(:,1))
    if y1(i,1) > pi/2 && y1(i,2) > 3.5
        z = i; %index of \theta greater than pi/2 & \dot{\theta} greater than 3.75
        break
    end
end

%% Two Methods: Pole Placement and LQR
% Controller Choice Options - Pole Placement (1) or LQR (0)

% c_stable = 0; % Hard Code

c_stable = input('Choose Either Pole Placement (1) or LQR (0): '); % Interactive


%-------------------------------------------------------------------------%
%                            POLE PLACEMENT                               %
%-------------------------------------------------------------------------%

if c_stable == 1    
% Pendulum SS Matrices - Linearized about \pi
% EOM: \ddot{\theta} = (-b/m/L)\dot{\theta} - (g/L)sin(\theta) + (1/m/L^2)u
theta = pi; % angle [rad] 
A = [0 1; -g/L*cos(theta) -b/m/L];
B = [0; 1/m/L^2];

% Check Controllability of A Matrix
A_poles = eig(A); % open-loop eigenvalues
n = rank(ctrb(A,B)); % see if full rank - yes

% Pole Placement - Choose Eigenvalues
% p = [-2 + 1i; -2 - 1i]; % imaginary parts (oscillations) 
p = [-3; -4]; % no imaginary parts (no overshoot)
K = place(A, B, p); % gain matrix

% Solve ODE - Pole Placement
tspan = 0:dt:20; % time span [s]
y2o = [y1(z,1); y1(z,2)]; % initial conditions [rad;rad/s]
[t,y2] = ode45(@(t,y) pendyn(y,m,L,g,b,-K*(y - [pi; 0])),tspan,y2o);

% Find Value Repetition 
len = L*[sin(y2(:,1)) -cos(y2(:,1))];
for i = 1:length(t)-1
if abs(len(i,1) - len(i+1,1)) < 1e-6 || abs(len(i,2) - len(i+1,2)) < 1e-6
    a = i;
    break
end
end

% Pole Placement Controller Values
u2 = -K*(y2(1:a,:)' - [pi; 0]); % torque input
fprintf('\nPole Placement Controller Torque Range: [%4.2f %4.2f] \n\n',min(u2),max(u2)) 


%-------------------------------------------------------------------------%
%                                 LQR                                     %
%-------------------------------------------------------------------------%

elseif c_stable == 0
% Pendulum SS Matrices - Linearized about \pi
% EOM: \ddot{\theta} = (-b/m/L)\dot{\theta} - (g/L)sin(\theta) + (1/m/L^2)u
theta = pi; % angle [rad] 
A = [0 1; -g/L*cos(theta) -b/m/L];
B = [0; 1/m/L^2];

% Check Controllability of A Matrix
A_poles = eig(A); % open-loop eigenvalues
n = rank(ctrb(A,B)); % see if full rank - yes

% LQR
Q = [.1 0; 0 .01]; % state cost value
R = 1000; % input cost value - high bc we want a low control input
K = lqr(A, B, Q, R); % gain matrix

% Solve ODE - LQR
tspan = 0:dt:20; % time span [s]
y2o = [y1(z,1); y1(z,2)]; % initial conditions [rad;rad/s]
[t,y2] = ode45(@(t,y) pendyn(y,m,L,g,b,-K*(y - [pi; 0])),tspan,y2o);

% Find Value Repetition 
len = L*[sin(y2(:,1)) -cos(y2(:,1))];
for i = 1:length(t)-1
if abs(len(i,1) - len(i+1,1)) < 1e-6 || abs(len(i,2) - len(i+1,2)) < 1e-6
    a = i;
    break
end
end

% LQR Controller Values
u2 = -K*(y2(1:a,:)' - [pi; 0]); % torque input
fprintf('\nLQR Controller Torque Range: [%4.2f %4.2f]\n\n',min(u2),max(u2)) 

end

%% Simple Pendulum Animation

% Updated State Vectors for the Simple Pendulum Solution
y(:,1) = [y1(1:z,1)' y2(1:a,1)']'; % position \theta [rad]
y(:,2) = [y1(1:z,2)' y2(1:a,2)']'; % velocity \dot{\theta} [rad/s]
u = [u1(1:z)' u2]'; % input  
t = [t1(1:z)' (t(1:a)+t1(z))']'; % time vector [s]

% Origin
O = [0 0];
axis(gca,'equal'); %aspect ratio of the plot
axis([-1.25 1.25 -1.25 1.25]); %limit of the plot
xlabel('x [m]'); ylabel('y [m]');

% Loop for animation
for i = 1:length(u)
    % Mass point
    P = L*[sin(y(i,1)) -cos(y(i,1))];
    
    % Circle in origin
    O_circ = viscircles(O,0.01);
    
    % Pendulum
    pend = line([O(1) P(1)],[O(2) P(2)]);
    
    % Ball
    ball = viscircles(P,0.075);
    
    % Time interval to update the plot
    title(sprintf('Simple Pendulum at t = %4.2f',t(i)))
    pause(0.001);
    
    % Delete previous objects if it is not the final loop
    if k < length(u)
        delete(pend);
        delete(ball);
        delete(O_circ);
    end
end


%% Plot State Evolution Over Time

% Create Plots
close all; figure;
subplot(3,1,1)
plot(t,y(:,1))
xlabel('t [s]')
ylabel('$\theta$ [rad]','interpreter','latex')
axis([min(t) max(t) min(y(:,1)) max(y(:,1))])
subplot(3,1,2)
plot(t,y(:,2))
xlabel('t [s]')
ylabel('$\dot{\theta} [\frac{rad}{s}]$','interpreter','latex')
axis([min(t) max(t) min(y(:,2)) max(y(:,2))])
subplot(3,1,3)
plot(t,u)
xlabel('t [s]')
ylabel('u')
axis([min(t) max(t) min(u) max(u)])
sgtitle(sprintf('Simple Pendulum State Evolution Plot'))
