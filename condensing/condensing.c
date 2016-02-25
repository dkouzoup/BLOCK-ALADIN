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

#include "mex.h"

#include "linear_algebra.h"

extern data_struct data;
        
extern void printMatrix( const char* name, real_t* mat, unsigned nRows, unsigned nCols );

/* This contains an implementation of our block condensing algorithm. */
void block_condensing( ) {
	uint i, j, r, c;
    /* mexPrintf("%d \n", NVC); */
    
    /* Copy bound values */
    for( i = 0; i < NX; i++ ) data.lbU[i] = data.lb[i];
    for( i = 0; i < NX; i++ ) data.ubU[i] = data.ub[i];
    for( i = 0; i < N; i++ ) {
		for( j = 0; j < NU; j++ ) {
			data.lbU[NX+i*NU+j] = data.lb[i*(NX+NU)+NX+j];
			data.ubU[NX+i*NU+j] = data.ub[i*(NX+NU)+NX+j];
		}
	}
    
    /* Create matrix G, NOTE: this is a sparse matrix which is currently stored as a dense one! */
    /* propagate x0: */
    for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			data.C[j*N*NX+i] = data.A[j*NX+i]; /* A_0 */
		}
	}
    for( i = 1; i < N; i++ ) {
		propagateCX(&data.C[i*NX], &data.A[i*NX*NX]); /* G_{i,0} <- A_{i-1}*G_{i-1,0}  */
	}
    
    /* propagate controls: */
    for( j = 0; j < N; j++ ) {
		for( c = 0; c < NU; c++ ) {
			for( r = 0; r < NX; r++ ) {
					data.C[(NX+j*NU+c)*N*NX+j*NX+r] = data.B[j*NX*NU+c*NX+r]; /* B_j */
			}
		}
		for( i = j+1; i < N; i++ ) {
			propagateCU(&data.C[(NX+j*NU)*N*NX+i*NX], &data.A[i*NX*NX]); /* G_{i,0} <- A_{i-1}*G_{i-1,0}  */
		}
		if( j > 0 ) {
			propagatec(&data.d[j*NX], &data.A[j*NX*NX], &data.b[j*NX]);
		}
		else {
			for( i = 0; i < NX; i++ ) data.d[i] = data.b[i];
		}
	}
	
	/* correct lbA and ubA with d values */
	for( i = 0; i < N; i++ ) {
		for( j = 0; j < NX; j++ ) {
			data.lbA[i*NX+j] = data.lb[(i+1)*(NX+NU)+j] - data.d[i*NX+j];
			data.ubA[i*NX+j] = data.ub[(i+1)*(NX+NU)+j] - data.d[i*NX+j];
		}
	}
	
	/* !! Hessian propagation !! */
	/* propagate x0: */
	computeWx(&data.Q[N*NX*NX], &data.C[(N-1)*NX], &data.A[0]); /* A is unused for this operation because the W's are zero */
	for( i = N-1; i > 0; i-- ) {
		computeH_offDX(&data.Hc[NX+i*NU], &data.S[i*NX*NU], &data.C[(i-1)*NX], &data.B[i*NX*NU]);
		computeWx(&data.Q[i*NX*NX], &data.C[(i-1)*NX], &data.A[i*NX*NX]);
	}
	computeH_DX( );
	
	/* propagate controls: */
	for( j = 0; j < N; j++ ) {
		for( i = 0; i < NX*NU; i++ ) data.W2_u[i] = 0.0;
		computeWu(&data.Q[N*NX*NX], &data.C[(NX+j*NU)*N*NX+(N-1)*NX], &data.A[0]); /* A is unused here because W is zero */
		for( i = N-1; i > j; i-- ) {
			computeH_offDU(&data.Hc[(NX+j*NU)*NVC+NX+i*NU], &data.S[i*NX*NU], &data.C[(NX+j*NU)*N*NX+(i-1)*NX], &data.B[i*NX*NU]);
			computeWu(&data.Q[i*NX*NX], &data.C[(NX+j*NU)*N*NX+(i-1)*NX], &data.A[i*NX*NX]);
		}
		computeH_DU(&data.Hc[(NX+j*NU)*NVC+NX+j*NU], &data.R[j*NU*NU], &data.B[j*NX*NU]);
	}
	
	/* !! gradient propagation !! */
	computeWg(&data.g[N*(NX+NU)], &data.Q[N*NX*NX], &data.d[(N-1)*NX], &data.A[0]); /* A is unused for this operation because the W's are zero */
	for( i = N-1; i > 0; i-- ) {
		computeG_off(&data.gc[NX+i*NU], &data.g[i*(NX+NU)+NX], &data.S[i*NX*NU], &data.d[(i-1)*NX], &data.B[i*NX*NU]);
		computeWg(&data.g[i*(NX+NU)], &data.Q[i*NX*NX], &data.d[(i-1)*NX], &data.A[i*NX*NX]);
	}
	computeG( );
	
}




