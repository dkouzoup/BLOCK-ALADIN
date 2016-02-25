/*
 *  =======
 *  LICENSE
 *  =======
 
 *  This file is part of BLOCK ALADIN
 *  (https://github.com/dkouzoup/BLOCK-ALADIN).
 * 
 *  BLOCK ALADIN: A Block Based Augmented Lagrangian Algorithm for Highly
 *  Parallelizable Optimal Control.
 
 *  Copyright (C) 2015-2016 by Dimitris Kouzoupis and Rien Quirynen,
 *  Albert Ludwigs University of Freiburg and K.U.Leuven.
 *  Developed under the supervision of Boris Houska and Moritz Diehl.
 *  All rights reserved.
 * 
 *  BLOCK ALADIN is a free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 * 
 *  BLOCK ALADIN is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with BLOCK ALADIN; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA 
 *
 */
    
// authors: Dimitris Kouzoupis, Rien Quirynen
// date: 2015

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "mex.h"
#include "dual_step.h"

/** Internally used floating point type */
typedef double real_t;

#define FREE( mem ) { if( mem ) { mxFree( mem ); mem = NULL; } }

/* Define number of outputs */
#define NOO 5

data_struct data;

/** A bit more advanced printing function. */
void mexErrMsgTxtAdv(	char* string,
						...
						)
{
	static char buffer[ 128 ];
	
	va_list printArgs;
	va_start(printArgs, string);
	
	vsprintf(buffer, string, printArgs);
	va_end( printArgs );

	mexErrMsgTxt( buffer );
}

/** A simple helper function. */
void printMatrix(	const char* name,
					real_t* mat,
					unsigned nRows,
					unsigned nCols
					)
{
    unsigned r, c;
    mexPrintf("%s: \n", name);
    for (r = 0; r < nRows; ++r)
    {
        for(c = 0; c < nCols; ++c)
            mexPrintf("\t%f", mat[c * nRows + r]);
        mexPrintf("\n");
    }
}

/** A function for copying data from MATLAB to C array in COLUMN MAJOR. */
int getArray(	const unsigned mandatory,
				const mxArray* source,
				const int index,
				const char* name,
				real_t* destination,
				const unsigned nRows,
				const unsigned nCols
				)
{
	mxArray* mxPtr = mxGetField(source, index, name);
	unsigned i, j;
	double* dPtr;
	
	if (mxPtr == NULL)
	{
		if ( !mandatory )
			return -1;
		else
			mexErrMsgTxtAdv("Field %s not found.", name);
	}

    if ( !mxIsDouble( mxPtr ) )
		mexErrMsgTxtAdv("Field %s must be an array of doubles.", name);

    if (mxGetM( mxPtr ) != nRows || mxGetN( mxPtr ) != nCols )
		mexErrMsgTxtAdv("Field %s must be of size: %d x %d.", name, nRows, nCols);

	dPtr = mxGetPr( mxPtr );
	
	if (destination == NULL)
		destination = (real_t*)mxCalloc(nRows * nCols, sizeof( real_t ));

	if (nRows == 1 && nCols == 1)
		*destination = *dPtr;
	else
		for (i = 0; i < nRows; ++i)
			for (j = 0; j < nCols; ++j)
				destination[j * nRows + i] = (real_t)dPtr[j * nRows + i];
			
	return 0;
}

void setArray( 	mxArray* destination,
				const int index,
				const char* name,
				real_t* source,
				const unsigned nRows,
				const unsigned nCols
				)
{
	mxArray* mxPtr = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
	double* dPtr = mxGetPr( mxPtr );
	unsigned i, j;
	
	if (nRows == 1 && nCols == 1)
		*dPtr = *source;
	else
		for (i = 0; i < nRows; ++i)
			for(j = 0; j < nCols; ++j)
				dPtr[j * nRows + i] = (double)source[j * nRows + i];

	mxSetField(destination, index, name, mxPtr);
}

/** The MEX interface function. */
void mexFunction(	int nlhs,
					mxArray *plhs[],
					int nrhs,
					const mxArray *prhs[]
					)
{
    int ii;
    
	const mxArray* src = prhs[ 0 ];
	
	const char *infoNames[ 5 ] = {"status", "nullspaceTime", "PTime", "formTime", "solveTime"};
	mxArray* info;
	double status, solveTime;
    //double formTime[NB];
	mxArray* shPtr;
	timer tmr;
	
	const char *outNames[ NOO ];
    outNames[ 0 ] = "info";
	outNames[ 1 ] = "lambda";
    outNames[ 2 ] = "primalStep";
    outNames[ 3 ] = "Ld";
	outNames[ 4 ] = "Ll";
    
	if (nrhs != 1)
		mexErrMsgTxt(
			"This function requires exactly one input: a structure with parameters.");
			
	if (nlhs != 1)
		mexErrMsgTxt(
			"This function returns one output.");
			
	if( !mxIsStruct( src ) )
		mexErrMsgTxt("The function argument must be a structure.");

    memset(&data,0,sizeof(data_struct));

	/* Copy MATLAB arrays to C arrays. */
	getArray(1, src, 0, "C", data.C, NX, (NB-1)*(BS*NU+NX));
	getArray(1, src, 0, "H", data.H, BS*NU+NX, NB*(BS*NU+NX));
    getArray(1, src, 0, "g", data.g, NB*(BS*NU+NX), 1);
    getArray(1, src, 0, "nact", data.nact, NB, 1); 
    getArray(1, src, 0, "DT", data.DT, BS*NU+NX, NB*(BS*NU+NX));
    getArray(1, src, 0, "d", data.d, NB*(BS*NU+NX), 1); 
    getArray(1, src, 0, "c", data.c, (NB-1)*NX, 1);

    /* Fix sign convention */
    for ( ii = 0 ; ii < NB*(BS*NU+NX) ; ii++ )
        data.d[ii] = -data.d[ii];
    
    /* Copy DT matrices for QR factorization */
    for ( ii = 0 ; ii < NB*(BS*NU+NX)*(BS*NU+NX) ; ii++ )
        data.DT_tmp[ii] = data.DT[ii];

    /* Create an identity matrix */
    for ( ii = 0 ; ii < (BS*NU+NX) ; ii++ )
        data.eye[(BS*NU+NX)*ii+ii] = 1.0;
    
    /* GET NULLSPACE AND ROW SPACE OF EACH D_i */
    form_nullspace( );
    
    /* FORM ELIMINATION MATRICES */
    form_P( );
        
    /* FORM DUAL HESSIAN */
    form_dual_hessian( );
    /* ----------------------- */
    
    /* FORM RHS */
    form_rhs( );
    /* ----------------------- */


    /* CALL FACTORIZATION AND SOLVE ROUTINE */
    tic( &tmr );
    solve_dual_system( );
    solveTime = toc( &tmr );
    /* ----------------------- */   
	
    /* UPDATE PRIMAL VARIABLES */
    //tic( &tmr );
    update_primal( );
    //updateTime = toc( &tmr );
    /* ----------------------- */
    
	/* Prepare return argument */
	plhs[ 0 ] = mxCreateStructMatrix(1, 1, NOO, outNames);

    for ( ii = 0 ; ii < (NB-1)*NX ; ii++ )
        data.lambda[ii] = data.RHS[ii]; /* -data.m[ii]; */
    
    setArray(plhs[ 0 ], 0, "lambda", data.lambda, (NB-1)*NX, 1);
    setArray(plhs[ 0 ], 0, "primalStep", data.primalStep, NB*(BS*NU+NX),1);
    setArray(plhs[ 0 ], 0, "Ld", data.Md, NX, (NB-1)*NX);
	setArray(plhs[ 0 ], 0, "Ll", data.Ml, NX, (NB-2)*NX);

	/* Create the info structure. */
	info = mxCreateStructMatrix(1, 1, 5, infoNames);
		
	setArray(info, 0, "status", &status, 1, 1);
	setArray(info, 0, "solveTime", &solveTime, 1, 1);
    setArray(info, 0, "formTime", &data.formTimes, NB, 1);
	setArray(info, 0, "nullspaceTime", &data.nullspaceTimes, NB, 1);
	setArray(info, 0, "PTime", &data.PTimes, NB, 1);
	
	mxSetField(plhs[ 0 ], 0, "info", info);
}
