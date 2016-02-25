function generateMex(NX, NU, M, N)

% GENERATEMEX Function to generate mexed code for given problem dimensions
%
% =======
% INPUTS
% =======
%
% NX: Number of states
% NU: Number of inputs
% M:  Block size
% N:  Horizon length
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

% get seperator char from OS
home = pwd;
sep  = home(1);

% create names for generated files (with relative path)
storepath_c = [home sep 'condense_N'];
storepath_d = [home sep 'dualStep_N'];

% generate condensing
cd('condensing');

updateField('NX',NX);
updateField('NU',NU);

for ii = 1:length(M)
    updateField('N',M(ii));
    make_condensing([storepath_c num2str(M(ii))])
end

cd(home);

% generate dual step
cd('dual_step');

updateField('NX',NX);
updateField('NU',NU);

for ii = 1:length(M)
    updateField('NB',N/M(ii));
    updateField('BS',M(ii));
    make_dual_step([storepath_d num2str(M(ii))])
end

cd(home);

end

