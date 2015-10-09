
clc; close all;

% simulation parameters

N      = 80;        % horizon length
h      = 0.2;       % sampling time [s]
BS     = 10;        % block size
rho    = 1;         % regularization parameter
TOL    = 1e-8;      % tolerance for termination condition
SQP_it = 1;         % number of SQP iterations for NLP subproblems

% input and state constraints
v_max  = 0.25;
u_max  = 10;
du_max = 15;

% save problem data in struct
PROBLEM.nx   =  8;
PROBLEM.nu   =  2;
PROBLEM.Q    =  diag([6e-1, 1.5e-1, 5e-2, 1e-3, 1e-4, 1e-4, 3e-3, 1e-1]);
PROBLEM.R    =  diag([1e-4, 1e-4]);
PROBLEM.QN   =  PROBLEM.Q;
PROBLEM.x0   =  [0.0 0 0.8 0 0 0 0 0]';
PROBLEM.xref =  [0.55 0 0.4 0 0 0 0 0]';
PROBLEM.uref =  zeros(PROBLEM.nu,1);
PROBLEM.xLow = -[inf; v_max; inf; v_max; u_max; u_max; inf; inf];
PROBLEM.xUpp =  [inf; v_max; inf; v_max; u_max; u_max; inf; inf];
PROBLEM.uLow = -[du_max; du_max];
PROBLEM.uUpp =  [du_max; du_max];

% run block ALADIN
[ X, U, L, i] = BLOCK_ALADIN_fun(N, TOL, rho, SQP_it, BS, PROBLEM);

figure;
plotResults(h, N, X, U, PROBLEM);

niters = i;
disp(' ')
disp('---------------------------------------------------------------------------------------');
disp(['--> TOTAL NUMBER OF ITERATIONS: ' num2str(niters)]);
disp(' ')
disp('---------------------------------------------------------------------------------------');
disp(' ')
