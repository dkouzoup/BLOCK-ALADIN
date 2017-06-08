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

function generateSim()

clc;
close all;

h = 0.2;
N = 5;

acadoSet('problemname', 'crane');

DifferentialState xT vT xL vL uT uL phi omega;
Control duT duL;

tau1 = 0.012790605943772;
a1   = 0.047418203070092;
tau2 = 0.024695192379264;
a2   = 0.034087337273386;
g = 9.81;
c = 0.0;
m = 1318.0;

aT = -1.0/tau1*vT + a1/tau1*uT;
aL = -1.0/tau2*vL + a2/tau2*uL;

%% Differential Equation
f = dot([phi; omega]) == [ omega; ...
    1.0/xL*(-g*sin(phi) - aT*cos(phi) - 2*vL*omega - c*omega/(m*xL)) ];

A1 = zeros(6,6);
B1 = zeros(6,2);
A1(1,2) = 1.0;
A1(2,2) = -1.0/tau1;
A1(2,5) = a1/tau1;
A1(3,4) = 1.0;
A1(4,4) = -1.0/tau2;
A1(4,6) = a2/tau2;

B1(5,1) = 1.0;
B1(6,2) = 1.0;

numSteps = 1;

%% SIMexport
sim = acado.SIMexport( h );
sim.setLinearInput( A1, B1 );
sim.setModel( f );

sim.set( 'INTEGRATOR_TYPE',           'INT_IRK_GL4'     );
sim.set( 'DYNAMIC_SENSITIVITY',       'FORWARD'         );
sim.set( 'NUM_INTEGRATOR_STEPS',      numSteps          );

sim.exportCode('sim_export')
cd sim_export
make_acado_integrator('../integrate');
cd ..

end

