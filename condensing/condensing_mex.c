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
/*
input
{
	H
    g
    AB
    b
    lb
    ub
}

output
{
    info
	Hc
    gc
    Ac
    lbA
    ubA
    lbU
    ubU
    C
    c
}
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "mex.h"
#include "condensing.h"

#define FREE( mem ) { if( mem ) { mxFree( mem ); mem = NULL; } }

/* Define number of outputs */
#define NOO 10

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
	const mxArray* src = prhs[ 0 ];
	
	const char *infoNames[ 2 ] = {"status", "cpuTime"};
	mxArray* info;
	double status, cpuTime;
	mxArray* shPtr;
	timer tmr;
	uint i, j, k, k2;
	
	const char *outNames[ NOO ];
    outNames[ 0 ] = "info";
	outNames[ 1 ] = "Hc";
	outNames[ 2 ] = "gc";
	outNames[ 3 ] = "Ac";
	outNames[ 4 ] = "lbA";
	outNames[ 5 ] = "ubA";
	outNames[ 6 ] = "lbU";
	outNames[ 7 ] = "ubU";
	outNames[ 8 ] = "C";
	outNames[ 9 ] = "d";
	
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
	getArray(1, src, 0, "H", data.H, N*(NX+NU)*(NX+NU)+NX*NX, 1);
	getArray(1, src, 0, "g", data.g, N*(NX+NU)+NX, 1);
	getArray(1, src, 0, "AB", data.AB, NX, N*(NX+NU));
	getArray(1, src, 0, "b", data.b, N*NX, 1);
	getArray(1, src, 0, "lb", data.lb, N*(NX+NU)+NX, 1);
	getArray(1, src, 0, "ub", data.ub, N*(NX+NU)+NX, 1);
    
    /* copy H to Q, S and R */
    for( k = 0; k < N; k++ ) {
		for( j = 0; j < NX; j++ ) {
			for( i = 0; i < NX; i++ ) {
				data.Q[k*(NX*NX)+j*NX+i] = data.H[k*(NX+NU)*(NX+NU)+j*(NX+NU)+i];
			}
			for( i = 0; i < NU; i++ ) {
				data.S[k*(NX*NU)+j*NU+i] = data.H[k*(NX+NU)*(NX+NU)+j*(NX+NU)+NX+i];
			}
		}
		for( j = 0; j < NU; j++ ) {
			for( i = 0; i < NU; i++ ) {
				data.R[k*(NU*NU)+j*NU+i] = data.H[k*(NX+NU)*(NX+NU)+(NX+j)*(NX+NU)+NX+i];
			}
		}
	}
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			data.Q[N*(NX*NX)+j*NX+i] = data.H[N*(NX+NU)*(NX+NU)+j*NX+i];
		}
	}
	/* copy AB to A and B */
	for( k = 0; k < N; k++ ) {
		for( j = 0; j < NX; j++ ) {
			for( i = 0; i < NX; i++ ) {
				data.A[k*(NX*NX)+j*NX+i] = data.AB[k*NX*(NX+NU)+j*NX+i];
			}
		}
		for( j = 0; j < NU; j++ ) {
			for( i = 0; i < NX; i++ ) {
				data.B[k*(NX*NU)+j*NX+i] = data.AB[k*NX*(NX+NU)+(NX+j)*NX+i];
			}
		}
	}
// 	memset(&C, 0, sizeof( C ));
// 	memset(&d, 0, sizeof( d ));
// 	memset(&W1_x, 0, sizeof( W1_x ));
// 	memset(&W2_x, 0, sizeof( W2_x ));
// 	memset(&W1_u, 0, sizeof( W1_u ));
// 	memset(&W2_u, 0, sizeof( W2_u ));
// 	memset(&w1, 0, sizeof( w1 ));
// 	memset(&w2, 0, sizeof( w2 ));
    
	
    /* CALL CONDENSING ROUTINE */
	tic( &tmr );
    block_condensing( );
	cpuTime = toc( &tmr );
    /* ----------------------- */
    
	
	/* Make condensed Hessian a symmetric Matrix: */
	for( k = 0; k < N; k++ ) {
		for( j = 0; j < NX; j++ ) {
			for( i = 0; i < NU; i++ ) {
					data.Hc[(NX+k*NU+i)*NVC+j] = data.Hc[j*NVC+NX+k*NU+i];
			}
		}
	}
	for( k2 = 0; k2 < N-1; k2++ ) {
		for( k = k2+1; k < N; k++ ) {
			for( j = 0; j < NU; j++ ) {
				for( i = 0; i < NU; i++ ) {
					data.Hc[(NX+k*NU+i)*NVC+NX+k2*NU+j] = data.Hc[(NX+k2*NU+j)*NVC+NX+k*NU+i];
				}
			}
		}
	}
	
	/* Prepare return argument */
	plhs[ 0 ] = mxCreateStructMatrix(1, 1, NOO, outNames);
		
	setArray(plhs[ 0 ], 0, "Hc", data.Hc, NVC, NVC);
	setArray(plhs[ 0 ], 0, "gc", data.gc, NVC, 1);
	setArray(plhs[ 0 ], 0, "Ac", data.C, N*NX, NVC);
	setArray(plhs[ 0 ], 0, "lbA", data.lbA, N*NX, 1);
	setArray(plhs[ 0 ], 0, "ubA", data.ubA, N*NX, 1);
	setArray(plhs[ 0 ], 0, "lbU", data.lbU, NVC, 1);
	setArray(plhs[ 0 ], 0, "ubU", data.ubU, NVC, 1);
	setArray(plhs[ 0 ], 0, "C", data.C, N*NX, NVC);
	setArray(plhs[ 0 ], 0, "d", data.d, N*NX, 1);

	/* Create the info structure. */
	info = mxCreateStructMatrix(1, 1, 2, infoNames);
		
	setArray(info, 0, "status", &status, 1, 1);
	setArray(info, 0, "cpuTime", &cpuTime, 1, 1);
		
	mxSetField(plhs[ 0 ], 0, "info", info);
}
