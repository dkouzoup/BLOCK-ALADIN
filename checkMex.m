function checkMex( BS )

% CHECKMEX  Function to check if Mex files exist in directory
%
% =======
% INPUTS
% =======
%
% BS: Block size
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

cond_ok = 0;
dual_ok = 0;

listing = dir();

for ii = 1:length(listing)
    check_cond = strfind(listing(ii).name,'condense_N');
    check_dual = strfind(listing(ii).name,'dualStep_N');
    if check_cond == 1
        pos = strfind(listing(ii).name,'condense_N');
        pos = pos + numel('condense_N');
        siz = listing(ii).name(pos:pos+numel(num2str(BS))-1);
        if str2double(siz) == BS
            cond_ok = 1;
        end
    end
    if check_dual == 1
        pos = strfind(listing(ii).name,'dualStep_N');
        pos = pos + numel('dualStep_N');
        siz = listing(ii).name(pos:pos+numel(num2str(BS))-1);
        if str2double(siz) == BS           
            dual_ok = 1;
        end
    end
end
if dual_ok == 0 || cond_ok == 0
    error('Mexed files not found or compiled with wrong block size. Set GENERATE_MEX to 1.')
end

end

