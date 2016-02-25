function [ X, U, Lc, timelog ] = BLOCK_ALADIN_fun(N, TOL, rho, ALADIN_iter, BLOCKSIZE, PROBLEM)

% BLOCK_ALADIN_FUN A prototype implementation of the block ALADIN scheme
%                  for optimal control.
% 
% =======
% INPUTS
% =======
%
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
% timelog:       Structure with logging of timings
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
WARMSTART  = 1;     % to use qpOASES_sequence in each seperate block
active_TOL = 1e-8;  % tolerance for active inequality constraints
INFTY      = 1e12;  % infinity

% setup option struct for qpOASES in subproblems
qpOptions_NLP = qpOASES_options('MPC');

if WARMSTART == 1
    QPhandles = cell(N/BLOCKSIZE,1);
else
    QPhandles = [];
end

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

C = ones(nx,N*(nx+nu)); % to store concatenated sensitivities (qpDUNES-style)
cAbs = ones(N*nx,1);

DX = 1;
DU = 1;
i  = 0;
timelog = [];

zOpts   = []; % to store all solutions of stage NLPs
low_act = []; % flag for active lower bounds
upp_act = []; % flag for active upper bounds


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Build H_k, H_N, zref_k and zref_N (g_k, g_N are calculated online)
Htemp = [];
ztemp = [];
Hnlps = [];

for ii = 1:BLOCKSIZE
    Htemp = blkdiag(Htemp,Hi);
    ztemp = [ztemp; xref; uref];
end
Hstagek = blkdiag(Htemp,zeros(nx));

HstageN = blkdiag(Htemp,P);
zstagek = [ztemp; zeros(nx,1)];
zstageN = [ztemp; xref];

for ii = 1:N/BLOCKSIZE-1
    Hnlps = blkdiag(Hnlps, Hstagek);
end
Hnlps = blkdiag(Hnlps, HstageN);

% Build coupling matrix and vector for the sparse case

% first coupling matrix
E0 = zeros((N/BLOCKSIZE-1)*nx,ns);
E0(1:nx,end-nx+1:end) = -eye(nx);
E = E0;

% second coupling matrix (rest is derived by shifting this one)
Ek = zeros((N/BLOCKSIZE-1)*nx,ns);
Ek(1:nx,1:nx) = eye(nx);
if N/BLOCKSIZE >= 3
    Ek(nx+1:2*nx,end-nx+1:end) = -eye(nx);
    for ii = 2:N/BLOCKSIZE
        E  = [E Ek];
        Ek = circshift(Ek,nx,1);
        if ii == N/BLOCKSIZE-1
            Ek(Ek == -1) = 0; % as there is no -eye(nx) term in last Ek
        end
    end
else
    E = [E Ek];
end

e = zeros((N/BLOCKSIZE-1)*nx,1);

% Derivatives of box constraints for sparse case (order: x_up u_up x_low u_low)
D_states = [ eye(nx) zeros(nx,nu)];
D_inputs = [zeros(nu,nx)  eye(nu)];

D = [];
for ii = 1:BLOCKSIZE
    D = blkdiag(D,[D_states; D_inputs; -D_states; -D_inputs]);
end
D = blkdiag(D,[eye(nx);-eye(nx)]);

% store derivatives of inequalities in cells (as the dimensions are varying)
Dineq = {};
for ii = 1:N/BLOCKSIZE
    Dineq{ii} = D;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

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
                cAbs((k_index-1)*nx+1:k_index*nx) = states.value - states.sensX*Xl(:,k_index) - states.sensU*Ul(:,k_index);
                
                % log integration time
                timelog(i+1).integrationTime(NLP,k_index) = 1000*info.clockTime;
            end
                        
            % PERFORM AN SQP ITERATION:
            [zOpt, stat, lambda_nlp, mu_nlp, obj_nlp, itlog, QPhandles ] = ...
                iterate_NLP(Hi, Pi, gi, Ci, ci, ziLow, xLow, ziUpp, xUpp, k, BLOCKSIZE, nx, nu, x0, qpOptions_NLP, WARMSTART, QPhandles);
                        
            % log QP solution and block condensing times
            timelog(i+1).QPTime(NLP,k)         = 1000*itlog.QPTime;
            timelog(i+1).condensingTime(NLP,k) = 1000*itlog.condensingTime;
            
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
        
        % !! Check for active inequality constraints in this block: !!
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
            
            % log re-linearization times
            timelog(i+1).relinearize(k,k_i) = 1000*info.clockTime;
            
            k_index = (k-1)*BLOCKSIZE+k_i;
            
            % calculate c_k - xl_{k+1} for the relative dynamics
            if k_i < BLOCKSIZE
                x_next = Xl(:,k_index+1);
            else
                x_next = Sl;
            end
            cRel((k_i-1)*nx+1:k_i*nx,k) = states.value - x_next;
            
            C(:,(k_index-1)*nz+1:k_index*nz)   = [states.sensX states.sensU];
            cAbs((k_index-1)*nx+1:k_index*nx) = states.value - states.sensX*Xl(:,k_index) - states.sensU*Ul(:,k_index);
        end
        
    end
    
    % % % % % % % % % % % END DECOUPLING STEP % % % % % % % % % % % % % % %
   
    norm_diff = [norm([Xl(:,1:N)-X(:,1:N);Ul-U(:,1:N)])];
    disp(['Norm decoupled diff :  ' num2str(norm_diff) ' ']);
    
    % % % % %  % % % % % % % % CONSENSUS STEP % % % % % % % % % % % % % % %
    
    % get local NLP solutions and split them
    uindx = [repmat([zeros(nx,1); ones(nu,1)],BLOCKSIZE,1); zeros(nx,1)];
    uindx(1:nx) = 1;
    uOpts = zOpts(find(uindx),:);
    xOpts = zOpts(find(~uindx),:);
    
    uNLP  = uOpts(:);
    xNLP  = xOpts(:);
    zNLP  = zOpts(:);
    
    % get A,B matrices of the last linearization
    for ii = 1:N
        ABtmp(:,:,ii) = C(:,(ii-1)*(nx+nu)+1:ii*(nx+nu));
    end
    
    % use local NLP solutions to build gradient
    gnlps = [];
    for ii = 1:N/BLOCKSIZE
        gtemp = [];
        if ii < N/BLOCKSIZE
            gtemp = Hstagek*(zOpts(:,ii) - zstagek);
        else
            gtemp = HstageN*(zOpts(:,ii) - zstageN);
        end
        gnlps = [gnlps; gtemp];
    end
    
    % setup G term
    G = E*zNLP-e;
    
    upp_act(end-nx+1:end,1:end-1) = 0;
    low_act(end-nx+1:end,1:end-1) = 0;
    
    % get indicies of active constraints
    zCon = [];
    for ii = 1:N/BLOCKSIZE
        upp_tmp = upp_act(:,ii);
        low_tmp = low_act(:,ii);
        uppN    = upp_tmp(end-nx+1:end); upp_tmp(end-nx+1:end) = [];
        lowN    = low_tmp(end-nx+1:end); low_tmp(end-nx+1:end) = [];
        
        upp_tmp = reshape(upp_tmp,nx+nu,BLOCKSIZE);
        low_tmp = reshape(low_tmp,nx+nu,BLOCKSIZE);
        zC_tmp  = [upp_tmp; low_tmp];
        zC_tmp  = zC_tmp(:);
        
        zCon(:,ii) = [zC_tmp; uppN; lowN];
    end
    
    
    % block condense
    Hnlps_CON = [];
    gnlps_CON = [];
    
    for ii = 1:N/BLOCKSIZE
        
        ABtmpp = [];
        Htmp   =  Hnlps((ii-1)*ns+1:ii*ns,(ii-1)*ns+1:ii*ns);
        Htmpp  = [];
        
        for jj = 1:BLOCKSIZE
            ABtmpp     = [ABtmpp ABtmp(:,:,(ii-1)*BLOCKSIZE+jj)];
            stageindex = (jj-1)*(nx+nu)+1:jj*(nx+nu);
            Htmpp      = [Htmpp Htmp(stageindex,stageindex)];
        end
        
        HtmpN = Htmp(end-nx+1:end,end-nx+1:end);
        cinput.H = [Htmpp(:); HtmpN(:)];
        cinput.g = gnlps((ii-1)*ns+1:ii*ns,1);
        cinput.AB= ABtmpp;
        cinput.b = cRel(:,ii);
        cinput.lb= -inf(ns,1);
        cinput.ub= inf(ns,1);
        
        eval(['coutput = condense_N' num2str(BLOCKSIZE) '(cinput);']);
        
        timelog(i+1).blockCondenseTime(ii) = 1000*coutput.info.cpuTime;
        
        % X{ii}_opt = cdatas{ii}.A*U{ii}_opt + cdatas.C
        cdatas{ii}.A    = [eye(nx) zeros(nx,size(coutput.Ac,2)-nx); coutput.Ac];
        cdatas{ii}.C    = [zeros(nx,1); coutput.d-coutput.Ac*uOpts(:,ii)+xOpts(:,ii)];
        cdatas{ii}.Crel = [zeros(nx,1); coutput.d];
        
        Hnlps_CON  = blkdiag(Hnlps_CON,coutput.Hc);
        gnlps_CON  = [gnlps_CON; coutput.gc];
        
    end
    
    % Build coupling matrix and vector for the condensed case
    
    % initial block
    Ek = zeros((N/BLOCKSIZE-1)*nx,nb);
    Ek(1:nx,:) = -cdatas{1}.A(end-nx+1:end,:);
    E_CON = Ek;
    e_CON = zeros((N/BLOCKSIZE-1)*nx,1);
    e_CON(1:nx) = cdatas{1}.C(end-nx+1:end);
    
    for ii = 2:N/BLOCKSIZE
        % build and append intermediate blocks
        Ek = zeros((N/BLOCKSIZE-1)*nx,nb);
        Ek(1:nx,1:nx)   = eye(nx);
        if N/BLOCKSIZE >= 3
            Ek(nx+1:2*nx,:) = -cdatas{ii}.A(end-nx+1:end,:);
            Ek = circshift(Ek,(ii-2)*nx,1);
            if ii == N/BLOCKSIZE
                Ek(1:nx,:) = 0; % modify last block
            end
        end
        E_CON  = [E_CON Ek];
        if ii < N/BLOCKSIZE
            e_CON((ii-1)*nx+1:ii*nx)= cdatas{ii}.C(end-nx+1:end);
        end
    end
    
    % Derivatives of box constraints for the condensed case (order: x_up u_up x_low u_low)
    Dineq_CON = {};
    dineq_CON = {};
    for ii =1:N/BLOCKSIZE
        D_inputs = [zeros(nu,nx)  eye(nu) zeros(nu,(BLOCKSIZE-1)*nu)];
        D_temp   = [];
        d_temp   = [];
        for jj = 1:BLOCKSIZE+1
            D_states = cdatas{ii}.A((jj-1)*nx+1:jj*nx,:);
            d_states = cdatas{ii}.Crel((jj-1)*nx+1:jj*nx);
            if jj < BLOCKSIZE + 1
                D_temp = [D_temp; D_states; D_inputs; -D_states; -D_inputs];
                d_temp = [d_temp; d_states; zeros(nu,1); -d_states; zeros(nu,1)];
            else
                D_temp = [D_temp; D_states; -D_states];
                d_temp = [d_temp; d_states; -d_states];
            end
            D_inputs = circshift(D_inputs,nu,2);
        end
        Dineq_CON{ii} = D_temp;
        dineq_CON{ii} = d_temp;
    end
    
    G_CON = E_CON*uNLP-e_CON;
    
    % keep only derivatives of active constraints (condensed case)
    Diter_CON = Dineq_CON;
    diter_CON = dineq_CON;
    for ii = 1:N/BLOCKSIZE
        zTemp = zCon(:,ii);
        DTemp = Diter_CON{ii};
        dTemp = diter_CON{ii};
        rindx = zTemp == 0;
        DTemp(rindx,:) = [];
        dTemp(rindx)   = [];
        Diter_CON{ii}  = DTemp;
        diter_CON{ii}  = dTemp;
    end
    
    % add initial condition constraint
    Diter_CON{1} = [eye(nx) zeros(nx,nb-nx); Diter_CON{1}];
    diter_CON{1} = [zeros(nx,1); diter_CON{1}];
    
    % restructure data to pass in mex file
    inputd.C = [];
    for ii = 1:N/BLOCKSIZE-1
        inputd.C = [inputd.C E_CON((ii-1)*nx+1:ii*nx,(ii-1)*nb+1:ii*nb)];
    end
    inputd.c = -G_CON; % \tilde{c} = c - E*Z = - G_CON
    
    inputd.H = [];
    for ii = 1:N/BLOCKSIZE
        inputd.H = [inputd.H Hnlps_CON((ii-1)*nb+1:ii*nb,(ii-1)*nb+1:ii*nb)];
    end
    inputd.g = gnlps_CON;
    
    DT   = []; % all D's transposed [D_0' D_1' ...]
    d    = []; % all d's [d_0; d_1; ...]
    nact = []; % number of active constraints per block
    
    for ii = 1:(N/BLOCKSIZE)
        Dtmp = Diter_CON{ii};
        dtmp = diter_CON{ii};
        DT   = [DT Dtmp'];
        d    = [d; dtmp];
        nact = [nact; size(Dtmp,1)];
    end
    DT = [DT zeros((BLOCKSIZE*nu+nx), ((N/BLOCKSIZE)*(BLOCKSIZE*nu+nx)*(BLOCKSIZE*nu+nx)-numel(DT))/(BLOCKSIZE*nu+nx),1)];
    d  = [d; zeros((BLOCKSIZE*nu+nx)*(N/BLOCKSIZE)-sum(nact),1)];
    
    inputd.DT   = DT;
    inputd.d    = d;
    inputd.nact = nact;
    
    eval(['outputd = dualStep_N' num2str(BLOCKSIZE) '(inputd);'])
    
    % log timings
    timelog(i+1).formDualHessianTime       = 1000*outputd.info.formTime';
    timelog(i+1).factorizeSolveTime        = 1000*outputd.info.solveTime;
    timelog(i+1).formEliminationMatrixTime = 1000*outputd.info.PTime';
    timelog(i+1).nullspaceTime             = 1000*outputd.info.nullspaceTime';
    
    % update primal and dual variables
    dualVars_CON   = outputd.lambda;
    primalStep_CON = outputd.primalStep;
    
    % calculate new primal iterate and expand solution
    uVars = primalStep_CON + uNLP;
    uVars = reshape(uVars,nu*BLOCKSIZE+nx,N/BLOCKSIZE);
    utmps = uVars(nx+1:end,:);
    
    primalVars_CON = [];
    for ii = 1:N/BLOCKSIZE
        xtmp = cdatas{ii}.A*uVars(:,ii)+cdatas{ii}.C;
        xtmp = reshape(xtmp,nx,BLOCKSIZE+1);
        if ii == N/BLOCKSIZE
            xNtmp = xtmp(:,end);
        end
        xtmp = xtmp(:,1:end-1);
        utmp = reshape(utmps(:,ii),nu,BLOCKSIZE);
        ztmp = [xtmp; utmp];
        ztmp = ztmp(:);
        primalVars_CON = [primalVars_CON; ztmp];
    end
    primalVars_CON =[primalVars_CON; xNtmp];
    
        
    % % % % %  % % % % % % END CONSENSUS STEP % % % % % % % % % % % % % % %

    zOpt = primalVars_CON;
    Lc   = -reshape(dualVars_CON,nx,N/BLOCKSIZE-1);
    
    % TAKE A FULL STEP:
    for k = 1:N
        temp = zOpt((k-1)*nz+1:(k-1)*nz+nx,1);
        X(:,k) = temp;
        
        temp = zOpt((k-1)*nz+nx+1:k*nz,1);
        U(:,k) = temp;
    end
    temp = zOpt(N*nz+1:end,1);
    X(:,N+1) = temp;
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
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

% clean qpOASES memory
if WARMSTART == 1
    for ii = 1: N/BLOCKSIZE
        qpOASES_sequence( 'c',QPhandles{ii});
    end
end

end

