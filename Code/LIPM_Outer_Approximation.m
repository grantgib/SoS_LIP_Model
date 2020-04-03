%% LIPM Outer Approximation - 0-step viable capture region
clear; clc;
restoredefaultpath;
addpath(genpath('../Toolboxes/SOSTOOLS.303'));
addpath(genpath('../Toolboxes/sedumi'));

%% =============================================
clear; clc;

% Initialize symbolics and other variables
syms x1 x2 p sigma_R sigma_p sigma_n t 'real'
Q = sym('Q',[2,2]);
x = [x1; x2];
% V = x'*Q*x; %note since this is 2 states it's actually SOS

R = 1;
disp('Symbolics');

%% =============================================
% Dynamics in control affine form --> xdot = f(x) + g(x)*u. (f is the drift vector)
grav = 9.81;       % gravity
r_foot = 0.05;  % stance foot max width
z_bar = 1;      % CoM height
% xdot = [x2; (grav/z_bar)*(x1 + r_foot*u1)];
f_original = [x2; (grav/z_bar)*x1];
g_original = [0; (grav/z_bar)*r_foot];
% scale = (grav/z_bar);
% f_tilde = f_original./scale;
% g_tilde = g_original./scale;

disp('Dynamics');

%% =============================================
% First, initialize the sum of squares program
prog = sosprogram(x);
% prog = sosprogram([x;t]);

% Declare decision variables
vars =  [sigma_R; sigma_p; sigma_n];
% vars2 = [a;b];
prog = sosdecvar(prog, vars);

% The W(x): 
[prog,W] = sossosvar(prog,[x1; x2],'wscoeff');  % constr 5
% [prog,V] = sospolyvar(prog, monomials([1; x1; x2;t],2));
[prog,V] = sospolyvar(prog, monomials([1; x1; x2],2));
[prog,p] = sospolyvar(prog, monomials([x1; x2],2:4));

disp('SOS program and vars');

%% =============================================
% Constraints
jacV = jacobian(V,x);
constr1 = -jacV*f_original - 1*p - sigma_R*(R^2 - x'*x);     % need to define variables p, sigma_R
% constr1 = -jacobian(V,t) - jacV*f_original - 1*p - sigma_R*(R^2 - x'*x);     % need to define variables p, sigma_R
prog = sosineq(prog,constr1);

constr2 = subs(V,[x1;x2],[0;0]); %
prog = sosineq(prog,constr2);

constr3 = p - jacV*g_original - sigma_p*(R^2-x'*x); % assume m = 1
prog = sosineq(prog,(constr3));

constr4 = p + jacV*g_original - sigma_n*(R^2-x'*x);
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









































