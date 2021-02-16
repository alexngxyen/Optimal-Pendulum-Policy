% ECE 594D Project - Double Pendulum Animation 
% By: Alex Nguyen

clc; clear; close all;

%% Deriving Postition, Velocity, and Acceleration
% Assume: rods are massless!
syms theta_1(t) theta_2(t) L_1 L_2 m_1 m_2 g

% Displacements
x_1 = L_1*sin(theta_1);
y_1 = -L_1*cos(theta_1);
x_2 = x_1 + L_2*sin(theta_2);
y_2 = y_1 - L_2*cos(theta_2);

% Velocities and Accelerations
vx_1 = diff(x_1); ax_1 = diff(vx_1);
vy_1 = diff(y_1); ay_1 = diff(vy_1);
vx_2 = diff(x_2); ax_2 = diff(vx_2);
vy_2 = diff(y_2); ay_2 = diff(vy_2);

%% Define Equations of Motion

% Rod Tension
syms T_1 T_2

% Mass 1 Forces
eqx_1 = m_1*ax_1(t) == -T_1*sin(theta_1(t)) + T_2*sin(theta_2(t));
eqy_1 = m_1*ay_1(t) == T_1*cos(theta_1(t)) - T_2*cos(theta_2(t)) - m_1*g;

% Mass 2 Forces
eqx_2 = m_2*ax_2(t) == -T_2*sin(theta_2(t));
eqy_2 = m_2*ay_2(t) == T_2*cos(theta_2(t)) - m_2*g;

% Evalutate Forces
Tension = solve([eqx_1 eqy_1],[T_1 T_2]); % solve for tension forces
eqRed_1 = subs(eqx_2,[T_1 T_2],[Tension.T_1 Tension.T_2]);
eqRed_2 = subs(eqy_2,[T_1 T_2],[Tension.T_1 Tension.T_2]);

%% Solve System of Equations

% Define Pendulum Parameters
L_1 = 1; L_2 = 1.5; % length [m]
m_1 = 1; m_2 = 2; % mass [kg]
g = 9.81; % gravity [m/s^2]

% Reduced Equations
eqn_1 = subs(eqRed_1);
eqn_2 = subs(eqRed_2);

% Convert to First Order ODEs from Second Order then to MATLAB Function
[V,S] = odeToVectorField(eqn_1,eqn_2); % First Order ODE
pend = matlabFunction(V,'vars',{'t','Y'}); % MATLAB function

% Solve for State Outputs
IC = [pi/4; 0; pi/3; 0];
tspan = [0 10]; % time span [s]
dpen = ode45(pend,tspan,IC);

%% Plot Solutions

% Defining Variables from Structure
t = dpen.x'; % time [s]
x = dpen.y'; % states [rad] or [rad/s]

% Each State on Subplot
list = {'$\theta_2$'; '$\frac{d\theta_2}{dt}$'; '$\theta_1$'; '$\frac{d\theta_1}{dt}$'};
figure(1); 
for i = 1:4
    subplot(4,1,i)
    plot(t,x(:,i))
    xlabel('t [s]')
    ylabel(list{i},'interpreter','latex')
end
sgtitle('Double Pendulum State Output')

% All States on One Figure
figure(2); 
plot(t,x)
legend('\theta_2','d\theta_2/dt','\theta_1','d\theta_1/dt')
title('Double Pendulum State Output')
xlabel('t [s]')
ylabel('Output [rad] or [rad/s]')

%% Create Animations

% Double Pendulum Oscillation
x_1 = @(t) L_1*sin(deval(dpen,t,3));
y_1 = @(t) -L_1*cos(deval(dpen,t,3));
x_2 = @(t) L_1*sin(deval(dpen,t,3))+L_2*sin(deval(dpen,t,1));
y_2 = @(t) -L_1*cos(deval(dpen,t,3))-L_2*cos(deval(dpen,t,1));

% Stop Motion Animation Object
figure(3);
fanimator(@(t) plot(x_1(t),y_1(t),'ro','MarkerSize',m_1*10,'MarkerFaceColor','r'));
axis equal; hold on;
fanimator(@(t) plot([0 x_1(t)],[0 y_1(t)],'r-'));
fanimator(@(t) plot(x_2(t),y_2(t),'go','MarkerSize',m_2*10,'MarkerFaceColor','g'));
fanimator(@(t) plot([x_1(t) x_2(t)],[y_1(t) y_2(t)],'g-'));
fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)));
hold off;


playAnimation
