
function updateField(arg, value)

% UPDATEFIELD Set one of the following fields N,NX,NU to given value and update NVC
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

f   = fileread('condensing.h');

if numel(arg) == 1
    arg = [arg ' '];
end

if ~ischar(value)
    value = num2str(value);
end

value = ['#define ' arg ' ' value '               '];

argpos = strfind(f,['#define ' arg]);

f(argpos:argpos+numel(value)-1) = value;

% update NVC

argpos = strfind(f,'#define N ');
N = str2double(f(argpos+10:argpos+20));
argpos = strfind(f,'#define NX ');
NX = str2double(f(argpos+10:argpos+20));
argpos = strfind(f,'#define NU ');
NU = str2double(f(argpos+10:argpos+20));

NVC = N*NU+NX;

argpos = strfind(f,'#define NVC ');
value  = ['#define NVC  ' num2str(NVC) '        '];
f(argpos:argpos+numel(value)-1) = value;

% Write new file and overwrite

fid = fopen('temp.h','wt');
fprintf(fid,'%s\n',f);
fclose(fid);

copyfile('temp.h', 'condensing.h');
delete('temp.h');

end




