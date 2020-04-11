%% Plot Results
clear; clc; close all;

% Load Results
syms x1 x2 t real
sol_0step = load('0step_outer');
sol_1step = load('1step_outer');
V_0step = sol_0step.V_opt;
V_1step = sol_1step.V_opt;

V_1max = double(subs(V_1step,[x1,x2,t],[0,0,0]))

%% Plots
figure;
% fs0 = fsurf(V_0step + V_1max,'EdgeColor','none');
fs0 = fsurf(V_0step,'EdgeColor','none');
fs0.FaceColor = 'k';
fs0.FaceAlpha = 0.5;
hold on;
fs1 = fsurf(subs(V_1step,t,0), 'EdgeColor','none');
fs1.FaceColor = [0.9290, 0.6940, 0.1250];
fs1.FaceAlpha = 0.5;

sz = 25;
xlim([-0.5 0.5]); xlabel("$x_{cm}$",'interpreter','latex','FontSize',sz);
ylim([-1 1]); ylabel("$\dot{x}_{cm}$",'interpreter','latex','FontSize',sz); 
zlim([0 inf]); zlabel("V^*",'FontSize',sz/2); 



