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

#ifndef LIN_ALG_H
#define LIN_ALG_H

#include "condensing.h"

extern data_struct data;

static void propagatec( real_t* d_, real_t* A_, real_t* b_ ) {
	int i, j;
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			d_[i] += A_[j*NX+i]*d_[j-NX];
		}
		d_[j] += b_[j];
	}
}

static void propagateCU( real_t* C_, real_t* A_ ) {
	int i, j, k;
	for( j = 0; j < NU; j++ ) {
		for( i = 0; i < NX; i++ ) {
			for( k = 0; k < NX; k++ ) {
				C_[j*N*NX+i] += A_[k*NX+i]*C_[j*N*NX-NX+k];
			}
		}
	}
}

static void propagateCX( real_t* C_, real_t* A_ ) {
	int i, j, k;
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			for( k = 0; k < NX; k++ ) {
				C_[j*N*NX+i] += A_[k*NX+i]*C_[j*N*NX-NX+k];
			}
		}
	}
}

static void computeWx( real_t* Q_, real_t* C_, real_t* A_ ) {
	int i, j , k;
	for( i = 0; i < NX*NX; i++ ) data.W1_x[i] = data.W2_x[i];
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			data.W2_x[j*NX+i] = 0.0;
			for( k = 0; k < NX; k++ ) {
				data.W2_x[j*NX+i] += Q_[k*NX+i]*C_[j*N*NX+k];
				data.W2_x[j*NX+i] += A_[i*NX+k]*data.W1_x[j*NX+k];
			}
		}
	}
}

static void computeWu( real_t* Q_, real_t* C_, real_t* A_ ) {
	int i, j , k;
	for( i = 0; i < NX*NU; i++ ) data.W1_u[i] = data.W2_u[i];
	for( j = 0; j < NU; j++ ) {
		for( i = 0; i < NX; i++ ) {
			data.W2_u[j*NX+i] = 0.0;
			for( k = 0; k < NX; k++ ) {
				data.W2_u[j*NX+i] += Q_[k*NX+i]*C_[j*N*NX+k];
				data.W2_u[j*NX+i] += A_[i*NX+k]*data.W1_u[j*NX+k];
			}
		}
	}
}

static void computeH_offDX( real_t* Hc_, real_t* S_, real_t* C_, real_t* B_ ) {
	int i, j , k;
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NU; i++ ) {
			Hc_[j*NVC+i] = 0.0;
			for( k = 0; k < NX; k++ ) {
				Hc_[j*NVC+i] += S_[k*NU+i]*C_[j*N*NX+k];
				Hc_[j*NVC+i] += B_[i*NX+k]*data.W2_x[j*NX+k];
			}
		}
	}
}

static void computeH_offDU( real_t* Hc_, real_t* S_, real_t* C_, real_t* B_ ) {
	int i, j , k;
	for( j = 0; j < NU; j++ ) {
		for( i = 0; i < NU; i++ ) {
			Hc_[j*NVC+i] = 0.0;
			for( k = 0; k < NX; k++ ) {
				Hc_[j*NVC+i] += S_[k*NU+i]*C_[j*N*NX+k];
				Hc_[j*NVC+i] += B_[i*NX+k]*data.W2_u[j*NX+k];
			}
		}
	}
}

static void computeH_DX( ) {
	int i, j , k;
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NU; i++ ) {
			data.Hc[j*NVC+NX+i] = data.S[j*NU+i];
			for( k = 0; k < NX; k++ ) {
				data.Hc[j*NVC+NX+i] += data.B[i*NX+k]*data.W2_x[j*NX+k];
			}
		}
	}
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			data.Hc[j*NVC+i] = data.Q[j*NX+i];
			for( k = 0; k < NX; k++ ) {
				data.Hc[j*NVC+i] += data.A[i*NX+k]*data.W2_x[j*NX+k];
			}
		}
	}
}

static void computeH_DU( real_t* Hc_, real_t* R_, real_t* B_ ) {
	int i, j , k;
	for( j = 0; j < NU; j++ ) {
		for( i = 0; i < NU; i++ ) {
			Hc_[j*NVC+i] = R_[j*NU+i];
			for( k = 0; k < NX; k++ ) {
				Hc_[j*NVC+i] += B_[i*NX+k]*data.W2_u[j*NX+k];
			}
		}
	}
}


static void computeWg( real_t* q_, real_t* Q_, real_t* d_, real_t* A_ ) {
	int i, j;
	for( i = 0; i < NX; i++ ) {
		data.w1[i] = data.w2[i];
		data.w2[i] = q_[i];
	}
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			data.w2[i] += Q_[j*NX+i]*d_[j];
			data.w2[i] += A_[i*NX+j]*data.w1[j];
		}
	}
}

static void computeG_off( real_t* gc_, real_t* r_, real_t* S_, real_t* d_, real_t* B_ ) {
	int i, j;
	for( i = 0; i < NU; i++ ) gc_[i] = r_[i];
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NU; i++ ) {
			gc_[i] += S_[j*NU+i]*d_[j];
			gc_[i] += B_[i*NX+j]*data.w2[j];
		}
	}
}

static void computeG( ) {
	int i, j;
	/* x0 */
	for( i = 0; i < NX; i++ ) data.gc[i] = data.g[i];
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NX; i++ ) {
			data.gc[i] += data.A[i*NX+j]*data.w2[j];
		}
	}
	/* first control */
	for( i = 0; i < NU; i++ ) data.gc[NX+i] = data.g[NX+i];
	for( j = 0; j < NX; j++ ) {
		for( i = 0; i < NU; i++ ) {
			data.gc[NX+i] += data.B[i*NX+j]*data.w2[j];
		}
	}
}


#endif /* LIN_ALG_H */
