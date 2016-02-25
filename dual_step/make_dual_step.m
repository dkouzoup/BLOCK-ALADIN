function make_dual_step( name )

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

	% Output file name, and also function name
	if (nargin == 1)
		fileOut = name;
	else
		fileOut = 'dual_step';
	end;
		
	% Root folder of code generation
	CGRoot = '.';	
		
	% Auto-generated files
	CGSources = [ ...
		'CGRoot/dual_step_mex.c ' ...
		'CGRoot/dual_step.c ' ...
		'CGRoot/timing.c ' ...
		];
		
	% Adding additional linker flags for Linux
	ldFlags = '-lm -lstdc++';
	if (isunix() && ~ismac())
		ldFlags = '-lrt';
	end;

	% Recipe for compilation
	CGRecipe = [ ...
		'mex -v CC="gcc" CXX="g++"' ...
		' -I. -I''CGRoot''' ...
        ' ldFlags' ...
        ' COPTIMFLAGS="-Ofast -ffast-math -finline-functions" ' ...
        ' CXXOPTIMFLAGS="-Ofast -ffast-math -finline-functions" ' ...
		' -D__MATLAB__ CGSources -output %s.%s' ...
		];

% Compilation
CGSources = regexprep(CGSources, 'CGRoot', CGRoot);

CGRecipe = regexprep(CGRecipe, 'CGRoot', CGRoot);
CGRecipe = regexprep(CGRecipe, 'CGSources', CGSources);
CGRecipe = regexprep(CGRecipe, 'ldFlags', ldFlags);

% disp( sprintf( CGRecipe, fileOut, mexext ) ); 
fprintf( 'compiling... ' );
% disp( sprintf(CGRecipe, fileOut, mexext) );
eval( sprintf(CGRecipe, fileOut, mexext) );
fprintf( ['done! --> ' fileOut '.' mexext '\n'] );
