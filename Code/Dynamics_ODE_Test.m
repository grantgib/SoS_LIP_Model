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

%% Plots
Plot_States(t,x);

%% Animations
Animate_LIPM(t,x,args);

%% Functions
% ODE FUNCTION
function dx = lip_dynamics(t,x,args)
% Extract args
g = args.g;
z_bar = args.z_bar;
r_foot = args.r_foot;

% Compute Input
u1 = 0;

% State space
q = x(1);   % x center of mass position
dq = x(2);  % x center of mass velocity
ddq = g/z_bar*(q+r_foot*u1);

x = [q; dq];
dx = [dq; ddq];
end

% Plot position and velocity trajectory
function [] = Plot_States(t,x)
figure
subplot(1,2,1)
sgtitle('dynamics plot') % puts title over both
plot(t, x(:,1))
title('CoM position')
subplot(1,2,2)
plot(t,x(:,2))
title('CoM velocity')
end

% Animate LIPM
function [] = Animate_LIPM(t,x,args)
% Extract arguments
z_bar = args.z_bar;

% Initialize animation
f = figure;
ax = axes(f);
hold(ax,'on');
pos_stanceFoot = [0; 0];
scatter(ax,pos_stanceFoot(1), pos_stanceFoot(2),'or');
h1 = line([pos_stanceFoot(1) x(1,1)], [pos_stanceFoot(2) z_bar]);
h2 = scatter(ax,x(1,1), z_bar,'or');

% Animate for all t values
for i = 1:length(t)
    pos_COM = [x(i,1); z_bar];
    set(h1,'XData',[pos_stanceFoot(1) pos_COM(1)],'YData',[pos_stanceFoot(2) pos_COM(2)]);
    set(h2,'XData',pos_COM(1),'YData',pos_COM(2));
    drawnow;
    pause(0.1);
    % set(ax,'YLim',[0,1.6]);
    % set(ax,'XLim',[-1,6]);
end

end





