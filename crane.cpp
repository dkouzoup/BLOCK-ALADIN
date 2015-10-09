/*
*    This file is part of ACADO Toolkit.
*
*    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
*    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
*    Developed within the Optimization in Engineering Center (OPTEC) under
*    supervision of Moritz Diehl. All rights reserved.
*
*    ACADO Toolkit is free software; you can redistribute it and/or
*    modify it under the terms of the GNU Lesser General Public
*    License as published by the Free Software Foundation; either
*    version 3 of the License, or (at your option) any later version.
*
*    ACADO Toolkit is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public
*    License along with ACADO Toolkit; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*
*/


/**
*    Author David Ariens, Rien Quirynen
*    Date 2009-2013
*    http://www.acadotoolkit.org/matlab 
*/

#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado/utils/matlab_acado_utils.hpp>

USING_NAMESPACE_ACADO

#include <mex.h>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
 { 
 
    MatlabConsoleStreamBuf mybuf;
    RedirectStream redirect(std::cout, mybuf);
    clearAllStaticCounters( ); 
 
    mexPrintf("\nACADO Toolkit for Matlab - Developed by David Ariens and Rien Quirynen, 2009-2013 \n"); 
    mexPrintf("Support available at http://www.acadotoolkit.org/matlab \n \n"); 

    if (nrhs != 0){ 
      mexErrMsgTxt("This problem expects 0 right hand side argument(s) since you have defined 0 MexInput(s)");
    } 
 
    TIME autotime;
    DifferentialState xT;
    DifferentialState vT;
    DifferentialState xL;
    DifferentialState vL;
    DifferentialState uT;
    DifferentialState uL;
    DifferentialState phi;
    DifferentialState omega;
    Control duT;
    Control duL;
    DMatrix acadodata_M1;
    acadodata_M1.read( "crane_data_acadodata_M1.txt" );
    DMatrix acadodata_M2;
    acadodata_M2.read( "crane_data_acadodata_M2.txt" );
    SIMexport ExportModule1( 1, 0.2 );
    ExportModule1.set( GENERATE_MATLAB_INTERFACE, 1 );
    uint options_flag;
    options_flag = ExportModule1.set( INTEGRATOR_TYPE, INT_IRK_GL4 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: INTEGRATOR_TYPE");
    options_flag = ExportModule1.set( DYNAMIC_SENSITIVITY, FORWARD );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: DYNAMIC_SENSITIVITY");
    options_flag = ExportModule1.set( NUM_INTEGRATOR_STEPS, 1 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: NUM_INTEGRATOR_STEPS");
    DifferentialEquation acadodata_f1;
    acadodata_f1 << dot(phi) == omega;
    acadodata_f1 << dot(omega) == (-((-7.818238E+01)*vT+3.707268E+00*uT)*cos(phi)+(-9.810000E+00)*sin(phi)-2.000000E+00*omega*vL)/xL;

    ExportModule1.setLinearInput( acadodata_M1, acadodata_M2 );
    ExportModule1.setModel( acadodata_f1 );

    uint export_flag = 0;
    ExportModule1.setTimingSteps( 0 );
    export_flag = ExportModule1.exportCode( "sim_export" );
    if(export_flag != 0) mexErrMsgTxt("ACADO export failed because of the above error(s)!");


    clearAllStaticCounters( ); 
 
} 

