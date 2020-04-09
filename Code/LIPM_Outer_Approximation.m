%% LIPM Outer Approximation - 0-step viable capture region
clear; clc;
restoredefaultpath;
addpath(genpath('../Toolboxes/SOSTOOLS.303'));
addpath(genpath('../Toolboxes/sedumi'));

%% =============================================
clear; clc;
% Initialize symbolics and other variables
syms x1 x2 p sigma_R sigma_p sigma_n 'real'
x = [x1; x2];
R = 2;
disp('Symbolics');

%% =============================================
% Dynamics in control affine form --> xdot = f(x) + g(x)*u. (f is the drift vector)
% xdot = [x2; (grav/z_bar)*(x1 + r_foot*u1)];
grav = 9.81;       % gravity
r_foot = 0.05;  % stance foot max width
z_bar = 1;      % CoM height
f_x = [x2; (grav/z_bar)*x1];
g_x = [0; (grav/z_bar)*r_foot];
disp('Dynamics');

%% =============================================
% First, initialize the sum of squares program
prog = sosprogram(x);

% Declare decision variables
vars =  [sigma_R; sigma_p; sigma_n];
prog = sosdecvar(prog, vars);

poly_deg = 2;
if poly_deg == 4
    [prog,W] = sossosvar(prog,[1; x1; x2; x1^2; x1*x2; x2^2],'wscoeff');  % constr 5
    [prog,V] = sospolyvar(prog, monomials([1; x1; x2],poly_deg));
    [prog,p] = sospolyvar(prog, monomials([1; x1; x2],poly_deg));
elseif poly_deg == 2
    [prog,W] = sossosvar(prog,[1; x1; x2; x1^2; x1*x2; x2^2],'wscoeff');  % constr 5
    [prog,V] = sospolyvar(prog, monomials([1; x1; x2],poly_deg));
    [prog,p] = sospolyvar(prog, monomials([1; x1; x2],poly_deg));
end
disp('SOS program and vars');

%% =============================================
% Constraints
jacV = jacobian(V,x);
constr1 = -jacV*f_x - 1*p - sigma_R*(R^2 - x'*x);     % need to define variables p, sigma_R
prog = sosineq(prog,constr1);

constr2 = subs(V - 0.01,[x1;x2],[0;0]); %
prog = sosineq(prog,constr2);

constr3 = p - jacV*g_x - sigma_p*(R^2-x'*x); % assume m = 1
prog = sosineq(prog,(constr3));

constr4 = p + jacV*g_x - sigma_n*(R^2-x'*x);
prog = sosineq(prog,(constr4));

% constr6 = W - V - 1
prog = sosineq(prog,W-V-1);

% setup the objective function
W0 = int(int(W,x1,-R,R),x2,-R,R);  % W = c1 integral(x) + c2 integral(x^2)
prog = sossetobj(prog,W0);
disp('Constraints and Objective');

%% And call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
disp('Solver');

%% Extract solution values
V_opt = sosgetsol(prog,V);
W_opt = sosgetsol(prog,W);
p_opt = sosgetsol(prog,p);
sigmaR_opt = sosgetsol(prog,sigma_R);
sigmaP_opt = sosgetsol(prog,sigma_p);
sigmaN_opt = sosgetsol(prog,sigma_n);

disp("*************** V_OPT RESULTS *****************");
disp("Vopt = "); disp(V_opt)
[coeff_V,mono_V] = coeffs(V_opt, [x1 x2], 'All');
coeff_V = double(coeff_V); 
disp("Coefficients of V = "); disp(coeff_V);
disp("Monomials of V = "); disp(mono_V);
disp("*************** PLOT RESULTS *****************");
figure;
fsurf(V_opt);
xlim([-0.5 0.5]);
ylim([-1 1]);
zlim([0 inf]);
title("V_{opt} Outer Approximation (R = " + R + ")");
% 
% disp("*************** V_OPT w/ SCALE RESULTS *****************");
% % scale
% a = 1e6;
% b = 1e6;
% V_opt = subs(V_opt,[x1,x2],[a*(x1/z_bar),b*(x2/sqrt(grav*z_bar))]);
% disp("Vopt = "); disp(V_opt)
% [coeff_V,mono_V] = coeffs(V_opt, [x1 x2], 'All');
% coeff_V = double(coeff_V); 
% disp("Coefficients of V = "); disp(coeff_V);
% disp("Monomials of V = "); disp(mono_V);
% disp("*************** PLOT RESULTS *****************");
% figure;
% fsurf(V_opt);
% xlim([-0.5 0.5]);
% ylim([-1 1]);
% zlim([0 inf]); 
% title("V_{opt} Outer Approximation (R = " + R + ")");


































