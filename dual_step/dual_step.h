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

/** Block size */
#define BS 10                                                                        

/** Number of blocks (N/BLOCKSIZE) */
#define NB 8                                                                         

/** Number of differential variables */
#define NX 8                                                                       

/** Number of control inputs */
#define NU 2                                                               

typedef struct data_struct_ {
    
    /* INPUT DATA */    
    
    real_t C[(NB-1)*NX*(BS*NU+NX)];          /* Coupling matrices without the identity */
    real_t H[NB*(BS*NU+NX)*(BS*NU+NX)];      /* Hessian matrices */
     
    real_t g[NB*(BS*NU+NX)];                 /* Gradient vectors*/
    real_t c[(NB-1)*NX];                     /* Coupling vectors */
    
    real_t nact[NB];                         /* Number of active constraints per block */
    real_t d[NB*(BS*NU+NX)];                 /* Right hand sides d_i. Actual size is sum(nact) */
    real_t DT[NB*(BS*NU+NX)*(BS*NU+NX)];     /* Jacobians of active constraints D_i TRANSPOSED */
       
    /* INTERMEDIATE RESULTS */
    
    real_t DT_tmp[NB*(BS*NU+NX)*(BS*NU+NX)]; /* Copy of DT matrices to be overwritten with R factors in QR factorization */
    real_t RZ[NB*(BS*NU+NX)*(BS*NU+NX)];     /* Q_i matrices from QR factorization of DT_i 
                                              * Q_i = [R_i Z_i] with ncols(R_i) = nact(i) */
    real_t P_NEW[NB*(BS*NU+NX)*(BS*NU+NX)];  /* Elimination matrices */
    real_t Hred[NB*(BS*NU+NX)*(BS*NU+NX)];   /* Reduced Hessian matrices */
    real_t eye[(BS*NU+NX)*(BS*NU+NX)];       /* Identity matrix to use as nullspace when there are no active constraints */

    real_t Md[(NB-1)*NX*NX];                 /* Diagonal blocks of M */
    real_t Ml[(NB-2)*NX*NX];                 /* LOWER diagonal blocks of M */
    real_t CP[(NB-1)*NX*(BS*NU+NX)];         /* Only diagonal blocks of C*P, off-diagonal are trivial to compute*/
    
    real_t DR[NB*(BS*NU+NX)*(BS*NU+NX)];     /* Matrices D*R */
    real_t DRDR[NB*(BS*NU+NX)*(BS*NU+NX)];   /* Pos def matrices (D*R)'*(D*R) */
    real_t DRd[NB*(BS*NU+NX)];               /* (D*R)'*d  */
    real_t RDRd_NEW[NB*(BS*NU+NX)];          /* R*inv(DR)*d vectors */
    real_t HRDRdg[NB*(BS*NU+NX)];            /* Vector  H*RDRd+g */
    real_t PHRDRdg[NB*(BS*NU+NX)];           /* Vector P*(H*RDRd+g) */
    real_t RHS_temp[NB*(BS*NU+NX)];          /* Part of right hand side, before multiplying with C (RDRd - PHRDRdg) */
    real_t RHS[(NB-1)*NX];                   /* Right hand side */
    real_t primal_temp[NB*(BS*NU+NX)];       /* Intermediate result for primal step */
    
    /* OUTPUT DATA */
    
    real_t lambda[(NB-1)*NX];                /* New dual variables*/
    real_t primalStep[NB*(BS*NU+NX)];        /* Delta Z (primal variables) */
    real_t nullspaceTimes[NB];               /* Timings for nullspace and rowspace computation */
    real_t PTimes[NB];                       /* Timings for computation of elimination matrices */
    real_t formTimes[NB];                    /* Timings to form Hessian contributions */
    
    
} data_struct;

void form_dual_hessian( );

void form_rhs( );

void solve_dual_system( );

void update_primal( );

void calculateDiagonalBlock(real_t *Md, real_t *CP, const real_t *C, const real_t *P);

void matVecMult(real_t *res, const real_t *A, const real_t *x, const int m, const int n, const int trans);

void matVecMult_NV(real_t *res, const real_t *A, const real_t *x);

void matVecMult_NX(real_t *res, const real_t *A, const real_t *x);

void matTransVecMult_NX(real_t *res, const real_t *A, const real_t *x);

void matMatMult(real_t *C, const real_t *A, const real_t *B, const int m, const int n, const int p, const int transA, const int transB);

void CtimesVector(real_t *res, const real_t *x );

void CTtimesVector(real_t *res, const real_t *x );

void couplingTimesEliminationTruncated(real_t *res, const real_t *C, const real_t *P);

void blockCholeskyFactorization(real_t *Md, real_t *Ml);

void myblas_dtrsm(const int M, const int N, const real_t *A, real_t *B);

void myblas_dgemm(const int M, const int N, const int K, const real_t *A, 
        const real_t *B, real_t *C);

void blockBackwardSubstitution(real_t *Md, real_t *Ml, real_t *x);

void blockForwardSubstitution(real_t *Md, real_t *Ml, real_t *x);

void choleskyFactorization(const int *_n, real_t *a, const int *_lda, int *info);

void copyMatrix(real_t *B, const real_t *A, int nRows, int nCols);

void qr(real_t *Q, real_t *R, const int n);

real_t norm(const real_t *v, const int n);

void Qtimesw(real_t *res, const real_t *A, const real_t *x, int n);

void ZHZ(const real_t *Z, const real_t *H, real_t *res, int m, int n);

void matrixSubstitutionUpper(const int nz, const real_t *L, real_t *Bt);

void matrixSubstitutionLower(const int nz, const real_t *L, real_t *Bt);

void form_RDRd(const real_t *R, const real_t *D, real_t *d, const int dim, const int nblock);

void forwardSubstitution(const int dim, const real_t *L, real_t *v);

void backwardSubstitution(const int dim, const real_t *L, real_t *v);

#endif /* DIMENSIONS_H */


















