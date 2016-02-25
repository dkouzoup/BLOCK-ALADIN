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

#ifndef DIMENSIONS_H
#define DIMENSIONS_H

#include "timing.h"
#include <math.h>
#include <string.h>

/*
 * Common definitions
 */
/** Number of control/estimation intervals. */
#define N  10                                                              
/** Number of control variables. */
#define NU 2                                                                     
/** Number of differential variables. */
#define NX 8                                                               
/** !!!!! NVC = N*NU+NX !!!!! **/
#define NVC  28                                                            

typedef struct data_struct_ {
            
    /* INPUT DATA */
    real_t H[N*(NX+NU)*(NX+NU)+NX*NX];
    real_t g[N*(NX+NU)+NX];
    real_t AB[N*NX*(NX+NU)];
    real_t b[N*NX];
    real_t lb[N*(NX+NU)+NX];
    real_t ub[N*(NX+NU)+NX];

    /* OUTPUT DATA */
    real_t Hc[NVC*NVC];
    real_t gc[NVC];
    /*real_t Ac[N*NX*NVC];*/
    real_t lbA[N*NX];
    real_t ubA[N*NX];
    real_t lbU[NVC];
    real_t ubU[NVC];
    real_t C[N*NX*NVC];
    real_t d[N*NX];

    /* EXTRA DATA */
    real_t Q[(N+1)*NX*NX];
    real_t R[N*NU*NU];
    real_t S[N*NX*NU];
    real_t A[N*NX*NX];
    real_t B[N*NX*NU];

    real_t W1_x[NX*NX];
    real_t W2_x[NX*NX];
    real_t W1_u[NX*NU];
    real_t W2_u[NX*NU];

    real_t w1[NX];
    real_t w2[NX];
 
} data_struct;

void block_condensing( );

#endif /* DIMENSIONS_H */



























































































































































































































































































































































































































































