%% LIPM Outer Approximation - 0-step viable capture region
clear; clc;
restoredefaultpath;
addpath(genpath('../Toolboxes/SOSTOOLS.303'));
addpath(genpath('../Toolboxes/sedumi'));

%% =============================================
clear; clc;

% Initialize symbolics and other variables
syms x1 x2 p sigma_R sigma_p sigma_n 'real'
Q = sym('Q',[2,2]);
x = [x1; x2];
% V = x'*Q*x; %note since this is 2 states it's actually SOS

ep = 0.1;
R = 1;
disp('Symbolics');

%% =============================================
% Dynamics in control affine form --> xdot = f(x) + g(x)*u. (f is the drift vector)
grav = 9.81;       % gravity
r_foot = 0.05;  % stance foot max width
z_bar = 1;      % CoM height
% xdot = [x2; (grav/z_bar)*(x1 + r_foot*u1)];
f_original = [x2; (grav/z_bar)*x1];
g_original = [0; (grav/z_bar)];
% scale = (grav/z_bar);
% f_tilde = f_original./scale;
% g_tilde = g_original./scale;

disp('Dynamics');

%% =============================================
% First, initialize the sum of squares program
prog = sosprogram(x);

% Declare decision variables
vars =  [p; sigma_R; sigma_p; sigma_n];
syms a b;
% vars2 = [a;b];
prog = sosdecvar(prog, vars);

% The W(x): 
[prog,W] = sossosvar(prog,[x1; x2],'wscoeff');
[prog,V] = sossosvar(prog,[x1; x2],'wscoeff');

disp('SOS program and vars');

%% =============================================
% Constraints
jacV = jacobian(V,x);
constr1 = -jacV*f_original - 1*p - sigma_R*(R^2 - x'*x);     % need to define variables p, sigma_R
% constr1 = -jacV*f_tilde - 1*p - sigma_R*(R^2 - x'*x);     % need to define variables p, sigma_R
prog = sosineq(prog,constr1);

constr2 = V; %
prog = sosineq(prog,constr2);

constr3 = p - jacV*g_original - sigma_p*(R^2-x'*x);
% constr3 = p - jacV*g_tilde - sigma_p*(R^2-x'*x);
prog = sosineq(prog,(constr3));

constr4 = p + jacV*g_original - sigma_n*(R^2-x'*x);
% constr4 = p + jacV*g_tilde - sigma_n*(R^2-x'*x);
prog = sosineq(prog,(constr4));

% constr5 = W
prog = sosineq(prog,W);

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









































