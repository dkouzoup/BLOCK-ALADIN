function [ ] = plotResults( h, N, X, U, PROBLEM )

% =======
% INPUTS
% =======
%
% h: sampling time in [s]
% N: number of shooting intervals
% X: State trajectory
% U: Input trajectory
% PROBLEM: problem data (to extract constraint values and reference)
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

constr_values(1) = PROBLEM.xUpp(2);
constr_values(2) = PROBLEM.xUpp(5);
constr_values(3) = PROBLEM.uUpp(1);

xref = PROBLEM.xref;

subplot(2,4,1);
plot(h*[0:N], X(1,:)); hold on;
plot(h*[0 N], [xref(1) xref(1)], 'k--');
xlabel('t')
ylabel('xT')
title('Position trolley')

subplot(2,4,2);
plot(h*[0:N], X(2,:)); hold on;
plot(h*[0 N], [xref(2) xref(2)], 'k--');
plot(h*[0 N], -[constr_values(1) constr_values(1)], 'g--');
plot(h*[0 N], [constr_values(1) constr_values(1)], 'g--');
xlabel('t')
ylabel('vT')
title('Velocity trolley')

subplot(2,4,3);
plot(h*[0:N], X(5,:)); hold on;
plot(h*[0 N], -[constr_values(2) constr_values(2)], 'g--');
plot(h*[0 N], [constr_values(2) constr_values(2)], 'g--');
xlabel('t')
ylabel('uT')
title('Input trolley')

subplot(2,4,5);
plot(h*[0:N], X(3,:)); hold on;
plot(h*[0 N], [xref(3) xref(3)], 'k--');
xlabel('t')
ylabel('xL')
title('Position cable')

subplot(2,4,6);
plot(h*[0:N], X(4,:)); hold on;
plot(h*[0 N], [xref(4) xref(4)], 'k--');
plot(h*[0 N], -[constr_values(1) constr_values(1)], 'g--');
plot(h*[0 N], [constr_values(1) constr_values(1)], 'g--');
xlabel('t')
ylabel('vL')
title('Velocity cable')

subplot(2,4,7);
plot(h*[0:N], X(6,:)); hold on;
plot(h*[0 N], -[constr_values(2) constr_values(2)], 'g--');
plot(h*[0 N], [constr_values(2) constr_values(2)], 'g--');
xlabel('t')
ylabel('uL')
title('Input cable')

subplot(2,4,4);
stairs(h*[0:N-1], U(1,:), 'r'); hold on;
plot(h*[0 N-1], -[constr_values(3) constr_values(3)], 'g--');
plot(h*[0 N-1], [constr_values(3) constr_values(3)], 'g--');
xlabel('t')
ylabel('duT')
title('Input trolley controller')

subplot(2,4,8);
stairs(h*[0:N-1], U(2,:), 'r'); hold on;
plot(h*[0 N-1], -[constr_values(3) constr_values(3)], 'g--');
plot(h*[0 N-1], [constr_values(3) constr_values(3)], 'g--');
xlabel('t')
ylabel('duL')
title('Input cable controller')


end

