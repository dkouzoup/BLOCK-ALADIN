function [ X, U, Lc, i ] = BLOCK_ALADIN_fun(N, TOL, rho, ALADIN_iter, BLOCKSIZE, PROBLEM)

% BLOCK_ALADIN_FUN A prototype implementation of the block ALADIN scheme
%                  for optimal control.
% 
% =======
% INPUTS
% =======
%
% h:             Sampling time [s]
% N:             Number of shooting intervals
% TOL:           Tolerance for termination condition (distance from optimal solution)
% rho:           Regularization term (tuning parameter)
% ALADIN_iter:   Number of SQP iterations in stage-NLPs
% BLOCKSIZE:     Number of combined subproblems
% PROBLEM:       Structure with problem data
%
% =======
% OUTPUTS
% =======
%
% X:             Optimal state trajectory
% U:             Optimal input trajectory
% Lc             Optimal multipliers for linear coupling constraints
%
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

MAXIT      = 40;    % maximum number of ALADIN iterations
active_TOL = 1e-8;  % tolerance for active inequality constraints
INFTY      = 1e12;  % infinity

% add qpDUNES to path
qpDUNES_PATH = 'qpdunes_dev/interfaces/matlab';
addpath(qpDUNES_PATH);

% setup option struct
qpOptions_NLP = qpDUNES_options(...
    'default',              ...
    'printLevel',  0,       ...
    'logLevel',    1,       ...
    'maxIter',     200,     ...
    'lsType',      4,       ... % accelerated gradient biscection LS
    'regType',     2,       ... % regularize only singular directions; 1 is normalized Levenberg Marquardt
    'stationarityTolerance', 1e-6 ...
    );

qpOptions_CQP = qpOptions_NLP;

% load problem data
nu   = PROBLEM.nu;
nx   = PROBLEM.nx;
Q    = PROBLEM.Q;
R    = PROBLEM.R;
QN   = PROBLEM.QN;
x0   = PROBLEM.x0;
xref = PROBLEM.xref;
uref = PROBLEM.uref;
xLow = PROBLEM.xLow;
xUpp = PROBLEM.xUpp;
uLow = PROBLEM.uLow;
uUpp = PROBLEM.uUpp;

nz = nx+nu;                  % number of stage variables
ns = (nx+nu)*BLOCKSIZE + nx; % number of block variables
nb = nx + BLOCKSIZE*nu;      % number of block condensed variables

% convert MATLAB infs to INFTYs
xLow(xLow < -INFTY) = -INFTY;
xUpp(xUpp >  INFTY) =  INFTY;
uLow(uLow < -INFTY) = -INFTY;
uUpp(uUpp >  INFTY) =  INFTY;

ziLow = [xLow; uLow];
ziUpp = [xUpp; uUpp];

% load optimal solution if available
try
    eval(['load sol_' num2str(N) '.mat X U'])
    X_sol = X; U_sol = U;
    if size(X_sol,1) ~= nx || size(U_sol,1) ~= nu
        error('Wrong number of states or inputs in loaded solution')
    end
    if size(X_sol,2) ~= N+1
        error('Loaded solution with wrong horizon lengh')
    end
catch
    warning('Optimal solution from SQP not available, continue?')
    keyboard
end

disp(['ALGORITHM:        ' 'BLOCK-ALADIN'])
disp(['HORIZON:          ' num2str(N)])
disp(['BLOCK SIZE:       ' num2str(BLOCKSIZE)])
disp(['NUMBER OF BLOCKS: ' num2str(N/BLOCKSIZE)])

% INITIALIZE SOLUTION:
X  = repmat(xref,1,N+1);
U  = zeros(nu, N);
L  = zeros(nx, N);             % multipliers of all equality constraints
Lc = zeros(nx, N/BLOCKSIZE-1); % multipliers of coupling constraints


Hi = blkdiag(Q, R);
H  = repmat( Hi, 1, N );
P  = Q;
ziRef = [xref; uref];
g = [repmat( -Hi*ziRef, N, 1 ); -Q*xref ];
C = ones(nx,N*(nx+nu));
cFull = ones(N*nx,1);

DX = 1;
DU = 1;
i  = 0;

zOpts   = []; % to store all solutions of stage NLPs
low_act = []; % flag for active lower bounds
upp_act = []; % flag for active upper bounds


while max(DX, DU) > TOL && i <= MAXIT
    
    disp(['!!! ALADIN iteration ' num2str(i+1) ' !!!']);
    
    % Initialize upper and lower bounds to infinity
    zLow = -repmat( INFTY, N*(nx+nu)+nx, 1 );
    zUpp =  repmat( INFTY, N*(nx+nu)+nx, 1 );
        
    % % % % % % % % % % % % % DECOUPLING STEP % % % % % % % % % % % % % % %
    
    Xl = X;
    Ul = U;
    
    for k = 1:N/BLOCKSIZE
        for NLP = 1:ALADIN_iter
            
            gi = [];
            
            % -- for xi --:
            if k == 1
                Hi = repmat(blkdiag(Q, R), 1, BLOCKSIZE);
                gi = repmat([gi; -blkdiag(Q, R)*ziRef], BLOCKSIZE, 1);
            else
                Hi = repmat(blkdiag(Q, R), 1, BLOCKSIZE);
                Hi(1:nx,1:nx) = Hi(1:nx,1:nx)+rho*eye(nx);
                
                gi = repmat([gi; -blkdiag(Q, R)*ziRef], BLOCKSIZE, 1);
                gi(1:nz) = gi(1:nz) - [rho*X(:,(k-1)*BLOCKSIZE+1)+Lc(:,k-1);zeros(nu,1)];
            end
            
            % -- for si+1 --:
            if k == N/BLOCKSIZE
                Pi = Q;
                gi = [gi; -Q*xref];
            else
                Pi = rho*eye(nx);
                gi = [gi; -rho*X(:,k*BLOCKSIZE+1) + Lc(:,k)];
            end
            
            % -- the nonlinear equality constraint --:
            Ci = zeros(nx,BLOCKSIZE*nz);
            ci = zeros(BLOCKSIZE*nx,1);
            for k_i = 1:BLOCKSIZE
                input.x = Xl(:,(k-1)*BLOCKSIZE+k_i);
                input.u = Ul(:,(k-1)*BLOCKSIZE+k_i);
                [states info] = integrate(input);
                
                Ci(:,(k_i-1)*nz+1:k_i*nz) = [states.sensX states.sensU];
                ci((k_i-1)*nx+1:k_i*nx)   = states.value - states.sensX*Xl(:,(k-1)*BLOCKSIZE+k_i) - states.sensU*Ul(:,(k-1)*BLOCKSIZE+k_i);
                
                k_index = (k-1)*BLOCKSIZE+k_i;
                C(:,(k_index-1)*nz+1:k_index*nz)   = [states.sensX states.sensU];
                cFull((k_index-1)*nx+1:k_index*nx) = states.value - states.sensX*Xl(:,k_index) - states.sensU*Ul(:,k_index);

            end
                        
            lbi = [repmat(ziLow,BLOCKSIZE,1);xLow];
            ubi = [repmat(ziUpp,BLOCKSIZE,1);xUpp];
            
            if k == 1
                % add an initial value constraint (update on block index 0)
                lbi(1:nx) = x0;
                ubi(1:nx) = x0;
            end
            
            qpDUNES( 'init', BLOCKSIZE, Hi, Pi, gi, Ci, ci, ...
                lbi, ubi, [], [], [], qpOptions_NLP);
            
            % SOLVE:
            tic
            [zOpt, stat, lambda, mu, objFctnVal] = qpDUNES( 'solve' );
            qpSolveTime = toc;
                                    
            zOpts(:,k) = zOpt;
            
            % TAKE A FULL STEP:
            for k_i = 1:BLOCKSIZE
                Xl(:,(k-1)*BLOCKSIZE+k_i) = zOpt((k_i-1)*nz+1:(k_i-1)*nz+nx);
                Ul(:,(k-1)*BLOCKSIZE+k_i) = zOpt((k_i-1)*nz+nx+1:k_i*nz);
            end
            Sl = zOpt(end-nx+1:end);
            if k == N/BLOCKSIZE
                Xl(:,N+1) = zOpt(BLOCKSIZE*nz+1:BLOCKSIZE*nz+nx);
            end
        end
        
        % check for active inequality constraints in block
        block_low = [repmat(ziLow,BLOCKSIZE,1);xLow];
        err_low = abs(zOpt - block_low);
        ind_low = find(err_low < active_TOL);
        low_act(:,k) = err_low < active_TOL;
        for ind = 1:length(ind_low)
            if ind_low(ind) <= BLOCKSIZE*nz || k == N/BLOCKSIZE
                index = (k-1)*BLOCKSIZE*nz + ind_low(ind);
                zLow(index) = block_low(ind_low(ind));
                zUpp(index) = block_low(ind_low(ind));
            end
        end
        
        block_upp = [repmat(ziUpp,BLOCKSIZE,1);xUpp];
        err_upp = abs(block_upp - zOpt);
        ind_upp = find(err_upp < active_TOL);
        upp_act(:,k) = err_upp < active_TOL;
        for ind = 1:length(ind_upp)
            if ind_upp(ind) <= BLOCKSIZE*nz || k == N/BLOCKSIZE
                index = (k-1)*BLOCKSIZE*nz + ind_upp(ind);
                zLow(index) = block_upp(ind_upp(ind));
                zUpp(index) = block_upp(ind_upp(ind));
            end
        end
        
        % RE-LINEARIZE FOR THE CONSENSUS STEP:
        for k_i = 1:BLOCKSIZE
            input.x = Xl(:,(k-1)*BLOCKSIZE+k_i);
            input.u = Ul(:,(k-1)*BLOCKSIZE+k_i);
            [states info] = integrate(input);
            
            k_index = (k-1)*BLOCKSIZE+k_i;
            
            % calculate c_k - xl_{k+1} for the relative dynamics
            if k_i < BLOCKSIZE
                x_next = Xl(:,k_index+1);
            else
                x_next = Sl;
            end
            res((k_i-1)*nx+1:k_i*nx,k) = states.value - x_next;
            
            C(:,(k_index-1)*nz+1:k_index*nz)   = [states.sensX states.sensU];
            cFull((k_index-1)*nx+1:k_index*nx) = states.value - states.sensX*Xl(:,k_index) - states.sensU*Ul(:,k_index);
        end
    
    end
       
    norm_diff = [norm([Xl(:,1:N)-X(:,1:N);Ul-U(:,1:N)])];
    disp(['Norm decoupled diff :  ' num2str(norm_diff) ' ']);
    
    % % % % %  % % % % % % % % CONSENSUS STEP % % % % % % % % % % % % % % %
                 
    zLow(1:nx) = x0;
    zUpp(1:nx) = x0;
    qpDUNES( 'init', N, H, P, g, C, cFull, ...
        zLow, zUpp, [], [], [], qpOptions_CQP);
    
    % SOLVE:
    tic
    [zOpt, stat, lambda_abs, mu, objFctnVal] = qpDUNES( 'solve' );
    qpSolveTime = toc;
    
    % extract multipliers of coupling constraints
    L = reshape(lambda_abs,nx,N);
    Lc = L(:,BLOCKSIZE:BLOCKSIZE:N-BLOCKSIZE);
        
    % TAKE A FULL STEP:
    for k = 1:N
        temp = zOpt((k-1)*nz+1:(k-1)*nz+nx,1);
        X(:,k) = temp;
        
        temp = zOpt((k-1)*nz+nx+1:k*nz,1);
        U(:,k) = temp;
    end
    temp = zOpt(N*nz+1:end,1);
    X(:,N+1) = temp;
          
    norm_diff = [norm([Xl(:,1:N)-X(:,1:N);Ul-U(:,1:N)])];
    disp(['Norm consensus step:  ' num2str(norm_diff) ' ']);
    
    % COMPUTE MAX DISTANCE TO OPTIMAL SOLUTION:
    try
        DX = max(max(abs(X-X_sol)));
        DU = max(max(abs(U-U_sol)));
    catch
        warning('Termination condition by-passed (missing or inconsistent optimal solution)')
    end
    
    maxStep = max(DX,DU);
    
    disp(['Distance to solution:  ' num2str(maxStep) ' ']);
    disp(' ');
    
    i = i+1;
end


end

