%% LIPM Outer Approximation - N-step viable capture region
clear; clc;
restoredefaultpath;
addpath(genpath('../Toolboxes/SOSTOOLS.303'));
addpath(genpath('../Toolboxes/sedumi'));

%% =============================================
clear; clc;
% Initialize symbolics and other variables
syms x1 x2 p t s h 'real'
syms sigma_R sigma_T sigma_V sigma_s q_h sigma_Rp sigma_Tp sigma_Rn sigma_Tn 'real'
x = [x1; x2];
R = 2;  % 
T = 1; % step time
disp('Symbolics');

%% =============================================
% Dynamics in control affine form --> xdot = f(x) + g(x)*u. (f is the drift vector)
% xdot = [x2; (grav/z_bar)*(x1 + r_foot*u1)];
grav = 9.81;       % gravity
z_bar = 1;      % CoM height
r_foot = 0.05;  % stance foot max width
r_step = 0.7;
r = [x1-r_step*s;x2];
f_x = [x2; (grav/z_bar)*x1];
g_x = [0; (grav/z_bar)*r_foot];
disp('Dynamics');

%% =============================================
% First, initialize the sum of squares program
prog = sosprogram([x;t;s;h]);

% Declare decision variables
vars =  [sigma_R; sigma_T;
         sigma_V; sigma_s; q_h;
         sigma_Rp; sigma_Tp; 
         sigma_Rn; sigma_Tn];
prog = sosdecvar(prog, vars);
poly_deg = 2;
if poly_deg == 4
    [prog,W] = sossosvar(prog,[1; x1; x2; x1^2; x1*x2; x2^2],'wscoeff');  
    [prog,V] = sospolyvar(prog, monomials([1; x1; x2],poly_deg));
    [prog,p] = sospolyvar(prog, monomials([1; x1; x2],poly_deg));
elseif poly_deg == 2
    [prog,V_N] = sospolyvar(prog, monomials([1; x1; x2; t],poly_deg));
    [prog,W] = sossosvar(prog,[1; x1; x2; x1^2; x1*x2; x2^2],'wscoeff'); 
%     [prog,V_N1] = sospolyvar(prog, monomials([1; x1; x2; t],poly_deg));
    [prog,p] = sospolyvar(prog, monomials([1; x1; x2],poly_deg));
end
disp('SOS program and vars');

%% =============================================
jacVx = jacobian(V_N,x);
jacVt = jacobian(V_N,t);
V_0step = load('0step_outer','V_opt');
V_1step = load('1step_outer','V_opt');
V_N1 = V_0step.V_opt;

% Constraints
constr1 = -jacVx*f_x - jacVt-1*p - sigma_R*(R^2 - x'*x) - sigma_T*(T*t-t^2);     % need to define variables p, sigma_R
prog = sosineq(prog,constr1);

Vtemp = subs(V_N1,[t;x1;x2],[0;r(1);r(2)]);
constr2 = V_N - sigma_V*Vtemp - sigma_s*(1-s^2) - q_h*h;
% constr2 = -sigma_V.*Vtemp;
prog = sosineq(prog,constr2);

constr3 = p - jacVx*g_x - sigma_Rp*(R^2-x'*x) - sigma_Tp*(T*t-t^2); % assume m = 1
prog = sosineq(prog,(constr3));

constr4 = p + jacVx*g_x - sigma_Rn*(R^2-x'*x) - sigma_Tn*(T*t-t^2);
prog = sosineq(prog,(constr4));

constr5 = W - subs(V_N,t,0) - 1;
prog = sosineq(prog,constr5);
disp("Constraints");

% setup the objective function
W0 = int(int(W,x1,-R,R),x2,-R,R);  % W = c1 integral(x) + c2 integral(x^2)
prog = sossetobj(prog,W0);
disp('Objective');

%% And call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
disp('Solver');

%% Extract solution values
V_opt = sosgetsol(prog,V_N);

disp("*************** V_OPT RESULTS *****************");
scale = 0;
if scale
   V_opt = subs(V_opt,[x1,x2],[x1/z_bar, x2/sqrt(grav*z_bar)]); 
end
disp("Vopt = "); disp(V_opt)
[coeff_V,mono_V] = coeffs(V_opt, [x1 x2], 'All');
coeff_V = subs(coeff_V,t,0); 
coeff_V = double(coeff_V);
disp("Coefficients of V = "); disp(coeff_V);
disp("Monomials of V = "); disp(mono_V);

disp("*************** PLOT V RESULTS *****************");
figure;
fsurf(subs(V_opt,t,0));
% xlim([-0.5 0.5]); xlabel("$x_{cm}$",'interpreter','latex');
% ylim([-1 1]); ylabel("$\dot{x}_{cm}$",'interpreter','latex');
zlim([0 inf]); zlabel("V^*");
title("V_{opt} Outer Approximation (R = " + R + ", T = " + T + ")");

% save('1step_outer','V_opt')




























