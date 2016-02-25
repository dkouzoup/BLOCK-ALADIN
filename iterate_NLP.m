function [zOpt, stat, lambda, mu, objFctnVal, itlog, QPhandles] = iterate_NLP(Hi, Pi, gi, Ci, ci, ziLow, xLow, ziUpp, xUpp, iB, BLOCKSIZE, nx, nu, x0, qpOptions, WARMSTART, QPhandles)

%ITERATE_NLP Perform on SQP iteration, i.e., solve one QP, for the NLP
%            subproblem using the selected solver
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

nz = nx+nu;

lbi = [repmat(ziLow,BLOCKSIZE,1);xLow];
ubi = [repmat(ziUpp,BLOCKSIZE,1);xUpp];

if iB == 1
    % add an initial value constraint (update on block index 0)
    lbi(1:nx) = x0;
    ubi(1:nx) = x0;
end


% build Hessian
Hinew = [];
for ii = 1:BLOCKSIZE
    Hinew = blkdiag(Hinew,Hi(:,(ii-1)*nz+1:ii*nz));
end
Hinew = blkdiag(Hinew,Pi);

% build equality constraints
Aeq = zeros(BLOCKSIZE*nx,BLOCKSIZE*(nx+nu)+nx);
for kk = 1:BLOCKSIZE
    AB = Ci(:,(kk-1)*nz+1:kk*nz);
    st.sensX = AB(:,1:nx);
    st.sensU = AB(:,nx+1:end);
    Aeq((kk-1)*nx+1:kk*nx,(kk-1)*(nx+nu)+1:kk*(nx+nu)+nx) = [st.sensX st.sensU -eye(nx)];
end
beq = -ci;

input.H  = [Hi(:);Pi(:)];
input.g  = gi;
input.AB = Ci;
input.b  = ci;
input.lb = lbi;
input.ub = ubi;  %#ok<STRNU>

% call mexed condensing routine
eval(['output = condense_N' num2str(BLOCKSIZE) '(input);']);

condTime = output.info.cpuTime;

lbU = [lbi(1:nx); repmat(ziLow(nx+1:end),BLOCKSIZE,1)];
ubU = [ubi(1:nx); repmat(ziUpp(nx+1:end),BLOCKSIZE,1)];

% solve with qpOASES
if WARMSTART == 0
    [usol, objFctnVal, stat, iter, lambda, out] = ...
        qpOASES(output.Hc,output.gc,output.Ac,lbU,ubU,output.lbA,output.ubA, qpOptions);
    
elseif WARMSTART == 1
    if isempty(QPhandles{iB})
        [QPhandles{iB}, usol, objFctnVal, stat, iter, lambda, out] = ...
            qpOASES_sequence('i', output.Hc,output.gc,output.Ac,lbU,ubU,output.lbA,output.ubA, qpOptions);
    else
        [usol, objFctnVal, stat, iter, lambda, out] = ...
            qpOASES_sequence('m', QPhandles{iB}, output.Hc,output.gc,output.Ac,lbU,ubU,output.lbA,output.ubA, qpOptions);
    end
end

qpSolveTime = out.cpuTime;

xsol  = [usol(1:nx); output.Ac*usol + output.d];
usol  = reshape(usol(nx+1:end),nu,BLOCKSIZE);
xsol  = reshape(xsol,nx,BLOCKSIZE+1);
xsolN = xsol(:,end);
zsol  = [xsol(:,1:end-1); usol];
zOpt  = [zsol(:); xsolN];
mu    = [];

itlog.condensingTime = condTime;
itlog.QPTime         = qpSolveTime;

end

