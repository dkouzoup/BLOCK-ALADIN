% =======
% LICENSE
% =======
%
%    This file is part of BLOCK ALADIN 
%    (https://github.com/dkouzoup/BLOCK-ALADIN).
%
%    BLOCK ALADIN: A Block Based Augmented Lagrangian Algorithm for Highly 
%    Parallelizable Optimal Control.
%
%    Copyright (C) 2015-2016 by Dimitris Kouzoupis and Rien Quirynen, 
%    Albert Ludwigs University of Freiburg and K.U.Leuven.
%    Developed under the supervision of Boris Houska and Moritz Diehl. 
%    All rights reserved.
%
%    BLOCK ALADIN is a free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    BLOCK ALADIN is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with BLOCK ALADIN; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
%    MA  02110-1301  USA
%
%    \authors Dimitris Kouzoupis, Rien Quirynen
%    \date 2015

clear all; close all; clc

EXPORT_INTEGRATOR = 1;
GENERATE_MEX      = 1;

% Simulation parameters
rho   = 1;      % Scaling factor
SQPit = 1;      % Number of SQP iterations per subproblem
TOL   = 1e-8;   % Tolerance for termination condition

Ts = 0.2;   % Sampling time [s]
N  = 80;    % Horizon length
NX = 8;     % Number of states
NU = 2;     % Number of inputs
BS = 10;    % Block size

if GENERATE_MEX
    % generate mex files for condensing and dual step
    generateMex(NX,NU,BS,N);
else
    % check that files exist for the selected block size
    checkMex(BS);
end

if EXPORT_INTEGRATOR
    % generate integrator (requires ACADO toolkit)
    generateSim()
end

% constraints
v_max  = 0.25;
u_max  = 10;
du_max = 15;

% save problem data in struct
PROBLEM.nx   =  NX;
PROBLEM.nu   =  NU;
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

% run block ALADIN or block SQP
[ X, U, L, timelog ] = BLOCK_ALADIN_fun(N, TOL, rho, SQPit, BS, PROBLEM);

% plot tajectories
plotResults(Ts, N, X, U, PROBLEM);

% print information
niters = length(timelog);
disp(' ')
disp('---------------------------------------------------------------------------------------');
disp(['--> NUMBER OF ITERATIONS: ' num2str(niters)]);
disp(' ')
disp('---------------------------------------------------------------------------------------');
disp(' ')
