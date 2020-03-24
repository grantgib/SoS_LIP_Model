%% Dynamics ODE Test
clear; clc; close all;

%% Define LIPM parameters
args = struct;
args.g = 9.81;
args.r_foot = 0.05;
args.z_bar = 1;
args.r_step = 0.7;

%% Simulate Dynamics
tspan = [0 10];
x_init = [0.2; 0];
[t,x] = ode45(@(t,x) lip_dynamics(t,x,args), tspan, x_init);

%% Plots and Animations
figure
subplot(1,2,1)
sgtitle('dynamics plot') % puts title over both
plot(t, x(:,1))
title('CoM position')
subplot(1,2,2)
plot(t,x(:,2))
title('CoM velocity')

%% ODE FUNCTION
function dx = lip_dynamics(t,x,args)
    % Extract args
    g = args.g;
    z_bar = args.z_bar;
    r_foot = args.r_foot;

    %solving for the input
    u1 = 0;

    % State space
    q = x(1);   % x center of mass position
    dq = x(2);  % x center of mass velocity
    ddq = g/z_bar*(q+r_foot*u1);

    x = [q; dq];
    dx = [dq; ddq];

end


