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

% R = 2;
% R = 1; % too small
% R = 0.1; % too large
% R = 0.5; % too large
R = 0.1;
disp('Symbolics');

%% =============================================
% Dynamics in control affine form --> xdot = f(x) + g(x)*u. (f is the drift vector)
grav = 9.81;       % gravity
r_foot = 0.05;  % stance foot max width
z_bar = 1;      % CoM height

%scaling the states
xdim = x.*[1/z_bar;1/(grav*z_bar)^0.5];
% x = x.*[1/z_bar;1/(grav*z_bar)^0.5];

% xdot = [x2; (grav/z_bar)*(x1 + r_foot*u1)];
f_original = [xdim(2); (grav/z_bar)*xdim(1)];
g_original = [0; (grav/z_bar)*r_foot];


disp('Dynamics');

%% =============================================
% First, initialize the sum of squares program
prog = sosprogram(x);

% Declare decision variables
vars =  [sigma_R; sigma_p; sigma_n];
prog = sosdecvar(prog, vars);

use_quartic = 0;
if use_quartic
    [prog,W] = sossosvar(prog,[1; x1; x2; x1^2; x1*x2; x2^2],'wscoeff');  % constr 5
    [prog,V] = sospolyvar(prog, monomials([1; x1; x2],4));
    [prog,p] = sospolyvar(prog, monomials([1; x1; x2],4));
else
    [prog,W] = sossosvar(prog,[1; x1; x2],'wscoeff');  % constr 5
    [prog,V] = sospolyvar(prog, monomials([1; x1; x2],2));
    [prog,p] = sospolyvar(prog, monomials([1; x1; x2],2));
end
    disp('SOS program and vars');

%% =============================================
% Constraints
jacV = jacobian(V,x);
constr1 = -jacV*f_original - 1*p - sigma_R*(R^2 - xdim'*xdim);     % need to define variables p, sigma_R
prog = sosineq(prog,constr1);

constr2 = subs(V,[x1;x2],[0;0]); %
prog = sosineq(prog,constr2);

constr3 = p - jacV*g_original - sigma_p*(R^2-xdim'*xdim); % assume m = 1
prog = sosineq(prog,(constr3));

constr4 = p + jacV*g_original - sigma_n*(R^2-xdim'*xdim);
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
zlim([0 inf]); title("V_{opt} Outer Approximation (R = " + R + ")");

% disp("*************** W_OPT RESULTS *****************");
% disp("Wopt = "); disp(W_opt)
% [coeff_W,mono_W] = coeffs(W_opt, [x1 x2], 'All');
% coeff_W = double(coeff_W); 
% disp("Coefficients of W = "); disp(coeff_W);
% disp("Monomials of W = "); disp(mono_W);
% disp("*************** W_OPT PLOT RESULTS *****************");
% figure;
% fsurf(W_opt);
% xlim([-0.5 0.5]);
% ylim([-1 1]);
% zlim([0 inf]);
% 

































