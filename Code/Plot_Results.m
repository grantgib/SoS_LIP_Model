%% Plot Results
clear; clc; close all;

% Load Results
syms x1 x2 t real
sol_0step = load('0step_outer');
sol_1step = load('1step_outer');
V_0step = sol_0step.V_opt;
V_1step = sol_1step.V_opt;

% Plots
figure;
fs0 = fsurf(V_0step,'EdgeColor','none');
fs0.FaceColor = 'b';
hold on;
fs1 = fsurf(subs(V_1step,t,0), 'EdgeColor','none');
fs1.FaceColor = 'y';
fs1.FaceAlpha = 0.5;
xlim([-0.5 0.5]); xlabel("$x_{cm}$",'interpreter','latex');
ylim([-1 1]); ylabel("$\dot{x}_{cm}$",'interpreter','latex');
zlim([0 inf]); zlabel("V^*");




