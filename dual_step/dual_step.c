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
#include "dual_step.h"
#include <math.h>

extern data_struct data;

extern void printMatrix( const char* name, real_t* mat, unsigned nRows, unsigned nCols );

/* Copy a matrix A (nRows x nCols) to B */
void copyMatrix(real_t *B, const real_t *A, int nRows, int nCols)
{
	int i;

	for (i=0;i<nRows*nCols;i++)
		B[i] = A[i];
}


/* Multiply a matrix A(m x n) with a matrix B(n x p) and save to matrix C(m x p) */
void matMatMult(real_t *C, const real_t *A, const real_t *B, const int m, const int n, const int p, const int transA, const int transB)
{
    int i,j,k;
    
    if(transA == 0 && transB == 0)
    {
        for(k=0;k<p;k++)
        {
            for(j=0;j<n;j++)
            {
                for(i=0;i<m;i++)
                {
                    C[k*m+i] += A[j*m+i]*B[k*n+j];
                }
            }
        }
    }
    else if(transA == 1 && transB == 0)
    {
        for(k=0;k<p;k++)
        {
            for(j=0;j<n;j++)
            {
                for(i=0;i<m;i++)
                {
                    C[k*m+i] += A[i*n+j]*B[k*n+j];
                }
            }
        }
    }
    else if(transA == 0 && transB == 1)
    {
        for(k=0;k<p;k++)
        {
            for(j=0;j<n;j++)
            {
                for(i=0;i<m;i++)
                {
                    C[k*m+i] += A[j*m+i]*B[j*p+k];
                }
            }
        }
    }
    else if(transA == 1 && transB == 1)
    {
        for(k=0;k<p;k++)
        {
            for(j=0;j<n;j++)
            {
                for(i=0;i<m;i++)
                {
                    C[k*m+i] += A[i*n+j]*B[j*p+k];
                }
            }
        }
    }
}


/* Multiply a matrix A (m x n) with a vector x (n x 1) 
 * - if trans == 1, A' is of dimension m x n 
 * - res does not need to be initialized */
void matVecMult(real_t *res, const real_t *A, const real_t *x, const int m, const int n, const int trans)
{
    int i,j;
    
    if ( trans == 0 )
    {
        /* first column seperately to initialize result*/
        for( i = 0 ; i < m ; i++ )
        {
            res[i] = A[i]*x[0];
        }
        /* add all other columns to first one*/
        for( j = 1 ; j < n ; j++ )
        {
            for( i = 0 ; i < m ; i++ )
            {
                res[i] += A[j*m+i]*x[j];
            }
        }
    }
    else if ( trans == 1 )
    {
        for( j = 0 ; j < m ; j++ )
        {
            res[j] = 0;

            for( i = 0 ; i < n ; i++ )
            {
                res[j] += A[j*n+i]*x[i];
            }
        }

    }
    
}


/* Square matrix A(NV x NV) times vector, with NV = BS*NU+NX */
void matVecMult_NV(real_t *res, const real_t *A, const real_t *x)
{
    int i,j;
    
    /* first column seperately to initialize result*/
    for(i=0;i<BS*NU+NX;i++)
    {
        res[i] = A[i]*x[0];
    }
    /* add all other columns to first one*/
    for(j=1;j<BS*NU+NX;j++)
    {
        for(i=0;i<BS*NU+NX;i++)
        {
            res[i] += A[j*(BS*NU+NX)+i]*x[j];
        }
    }
    
}


/* Matrix A(NX x NV) times vector */
void matVecMult_NX(real_t *res, const real_t *A, const real_t *x)
{
    int i,j;
    
    /* first column seperately to initialize result*/
    for(i=0;i<NX;i++)
    {
        res[i] = A[i]*x[0];
    }
    /* add all other columns to first one*/
    for(j=1;j<BS*NU+NX;j++)
    {
        for(i=0;i<NX;i++)
        {
            res[i] += A[j*NX+i]*x[j];
        }
    }
    
}


/* Matrix A(NX x NV)' times vector */
void matTransVecMult_NX(real_t *res, const real_t *A, const real_t *x)
{
    int i,j;
    
    for(j=0;j<BS*NU+NX;j++)
    {
        res[j] = 0;
        for(i=0;i<NX;i++)
        {
            res[j] += A[j*NX+i]*x[i];
        }
    }
    
}


/* Get the nullspace and row space of each D_i with QR factorization of D_i' */
void form_nullspace( ) 
{
    int i;
    int debug = 0;
    timer tmr;
    
    /* index of accumulated active constraints */
    int totact = 0;
    
    for(i=0;i<NB;i++)
    {        
        if (data.nact[i] > 0)
        {
            tic( &tmr );          
            /* D_i' = Q_i*R_i, Q_i = [rowspace_i nullspace_i] */
            qr(&data.RZ[i*(BS*NU+NX)*(BS*NU+NX)], &data.DT_tmp[totact*(BS*NU+NX)], (int)data.nact[i]);
            data.nullspaceTimes[i] = toc( &tmr );
        }      
        totact += (int)data.nact[i];
    } 
}


/* Helper function for form_nullspace:
 * QR factorization of matrix A (m x n, stored in R) with Householder transformations
 * Assumptions: 
 * - m (fixed dimension = BS*NU+NX) >= n (varying dimension)
 * - Q initialized with zeros
 * - R contains matrix to be factorized on input and factor R on output */
void qr(real_t *Q, real_t *R, const int n)
{
	int i,j,k,s;
    
	real_t nrm, u, tau;
	real_t w[BS*NU+NX];
	real_t wR[n];
	real_t Qw[BS*NU+NX];

	/* initialize Q with identity matrix */
	for ( j = 0 ;  j < BS*NU+NX ; j++ )
	{
		Q[j*(BS*NU+NX)+j] = 1.0;
	}

	for ( j = 0 ; j < n ; j++ )
	{
		nrm = norm(&R[j*(BS*NU+NX)+j], BS*NU+NX-j);

		if ( R[j*(BS*NU+NX) + j ] >= 0 )
		{
			s = -1;
		}
		else
		{
			s = 1;
		}

		u = R[j*(BS*NU+NX)+j] -s*nrm;

		for ( i = 0 ; i < BS*NU+NX-j ; i++ )
		{
			w[i] = R[j*(BS*NU+NX)+j+i]/u;
		}
		w[0] = 1.0;

		tau = -(s*u)/nrm;

	    /* calculate w'*R */
		for ( i = 0 ; i < n ; i++ )
		{
			wR[i] = 0.0;
            
			for( k = 0 ; k < BS*NU+NX-j ; k++ )
			{
				wR[i] += w[k]*R[i*(BS*NU+NX)+k+j];
			}
		}

		/* update R */ 
		for ( i = 0 ; i < n ; i++ )
		{
			for ( k = 0 ; k < BS*NU+NX-j ; k++ ) 
			{
				R[i*(BS*NU+NX)+j+k] -= tau*w[k]*wR[i];
			}
		}

		/* calculate Q*w */
		Qtimesw(Qw, &Q[j*(BS*NU+NX)], w, BS*NU+NX-j);

		/* update Q*/
		for ( k = 0 ; k < BS*NU+NX-j ; k++ )
		{
			for ( i = 0 ; i < BS*NU+NX ; i++ )
			{
				Q[(k+j)*(BS*NU+NX)+i] -= Qw[i]*tau*w[k];
			}
		}
	}
}


/* Helper function for QR decomposition: 
 * return 2-norm of a vector v (n x 1) */
real_t norm(const real_t *v, const int n)
{
	real_t res = 0.0;

	int i;

	for (i=0;i<n;i++)
	{
		res += v[i]*v[i];
	}
	res = sqrt(res);
    
	return res;
}


/* Helper function for QR decomposition: 
 * multiply matrix Q (BS*NU+NX x n) with a vector w (n x 1) */
void Qtimesw(real_t *res, const real_t *A, const real_t *x, int n)
{
    int i,j;
    
    /* first column seperately to initialize result*/
    for( i = 0 ; i < BS*NU+NX ; i++ )
    {
        res[i] = A[i]*x[0];
    }
    /* add all other columns to first one*/
    for( j = 1 ; j < n ; j++ )
    {
        for( i = 0 ; i < BS*NU+NX ; i++ )
        {
            res[i] += A[j*(BS*NU+NX)+i]*x[j];
        }
    } 
}


/* This function generates the elimination matrices P_i */
void form_P( ) 
{
    int i,j;
    int info, nz;
    real_t *Z, *H;
    real_t Ztmp[(BS*NU+NX)*(BS*NU+NX)];
    timer tmr;
    
    for ( i = 0; i < NB ; i++ )
    {
        tic(&tmr);
        
        /* Size of nullspace */
        nz = BS*NU+NX - (int)data.nact[i];

        if ( data.nact[i] > 0 )
        {

            /* Pointers to Z and H */
            Z = &data.RZ[i*(BS*NU+NX)*(BS*NU+NX)+(int)data.nact[i]*(BS*NU+NX)];
            H = &data.H[i*(BS*NU+NX)*(BS*NU+NX)];
            
            /* printMatrix("Hessian", H, BS*NU+NX,BS*NU+NX);
             * printMatrix("Z", Z, BS*NU+NX,BS*NU+NX-(int)data.nact[i]); */

            /* Form reduced Hessian of block i */
            ZHZ(Z, H, &data.Hred[i*(BS*NU+NX)*(BS*NU+NX)], BS*NU+NX, nz);
            /* printMatrix("reduced Hessian", Hred, BS*NU+NX-(int)data.nact[i],BS*NU+NX-(int)data.nact[i]); */
        
            /* Factorize reduced Hessian (Hred <-- L) */
            choleskyFactorization(&nz, &data.Hred[i*(BS*NU+NX)*(BS*NU+NX)], &nz, &info);
       
            /* Avoid overwritting Z */
            copyMatrix(Ztmp, Z, BS*NU+NX, nz);

            /* Two customized matrix solves, Ztmp <-- (inv(Htilde)*Z')' */
            matrixSubstitutionUpper(nz, &data.Hred[i*(BS*NU+NX)*(BS*NU+NX)], Ztmp);
            matrixSubstitutionLower(nz, &data.Hred[i*(BS*NU+NX)*(BS*NU+NX)], Ztmp);
            
            /* Calculate P_i */
            matMatMult(&data.P_NEW[i*(BS*NU+NX)*(BS*NU+NX)], Ztmp, Z, BS*NU+NX, nz, BS*NU+NX, 0, 1);
   
        }
        else
        {            
            /* Reduced Hessian is identical to actual Hessian */
            copyMatrix(&data.Hred[i*(BS*NU+NX)*(BS*NU+NX)],&data.H[i*(BS*NU+NX)*(BS*NU+NX)],BS*NU+NX,BS*NU+NX);
            
            /* Factorize reduced Hessian (Hred <-- L) */
            choleskyFactorization(&nz, &data.Hred[i*(BS*NU+NX)*(BS*NU+NX)], &nz, &info);
       
            /* Perform solves with identity matrix to compute inverse */
            copyMatrix(Ztmp, data.eye, BS*NU+NX, nz);

            /* Two matrix solves */
            matrixSubstitutionUpper(nz, &data.Hred[i*(BS*NU+NX)*(BS*NU+NX)], Ztmp);
            matrixSubstitutionLower(nz, &data.Hred[i*(BS*NU+NX)*(BS*NU+NX)], Ztmp);
            
            /* Copy result to P */
            copyMatrix(&data.P_NEW[i*(BS*NU+NX)*(BS*NU+NX)],Ztmp,BS*NU+NX,BS*NU+NX);
            
        }
        
        data.PTimes[i] = toc(&tmr);
        
        /* printMatrix("P",  &data.P_NEW[i*(BS*NU+NX)*(BS*NU+NX)], BS*NU+NX,BS*NU+NX);
         * printMatrix("Real P", &data.P[i*(BS*NU+NX)*(BS*NU+NX)], BS*NU+NX,BS*NU+NX); */
        
    }
}


/* Helper function for form_P :
 * Form reduced Hessians H_tilde */
void ZHZ(const real_t *Z, const real_t *H, real_t *res, int m, int n)
{
	
	int i,j,k,l;
	real_t Mtmp[m*n];

	/* Build first column, that also computes intermediate results */
	for ( i=0; i<n; i++ )
	{
		for (k = 0; k < m; k++)
		{
			Mtmp[i*m+k] = 0.0;
			for (l = 0; l < m; l++)
			{
				Mtmp[i*m+k] += H[k*m+l]*Z[i*m+l];
			}
			res[i] += Z[k]*Mtmp[i*m+k];
		}
	}

	for ( j=1; j<n; j++ )
	{
		/* Copy symmetric part */
		for ( i=0; i<j; i++)
		{
			res[j*n+i] = res[i*n+j];
		}	
		/* Compute rest of the jth column */
		for ( i=j; i<n; i++ )
		{
			for (k = 0; k < m; k++)
			{
				res[j*n+i] += Z[j*m+k]*Mtmp[i*m+k];
			}
		}
	}
}


/* Helper function for form_P :
 * Matrix substitution A = L*B (solve with respect to B)
 * - for efficiency, we solve instead A' = B'*R (with respect to B') 
 * - Function reads in L and but uses its transpose (aka R)
 * - in our case A' = Z 
 * - Bt has A' (aka Z) at input and the solution of A = L*B transposed at output */
void matrixSubstitutionUpper(const int nz, const real_t *L, real_t *Bt)
{
    int i,j,k;
    
    for ( j = 0; j < nz; j++ )
    {
        for ( i = 0; i < BS*NU+NX; i++ )
        {
            for ( k = 0; k < j; k++ )
            {
                /* Bt[j*(BS*NU+NX)+i] -= Bt[k*(BS*NU+NX)+i]*R[j*nz+k]; */
                Bt[j*(BS*NU+NX)+i] -= Bt[k*(BS*NU+NX)+i]*L[k*nz+j];
            }
            /* Bt[j*(BS*NU+NX)+i] /= R[j*nz+j]; */
            Bt[j*(BS*NU+NX)+i] /= L[j*nz+j];
        }
    }
    
}


/* Helper function for form_P :
 * Matrix substitution A = L'*B (solve with respect to B)
 * - for efficiency we solve instead A' = B'*L (with respect to B') 
 * - in our case A' = Bt from previous solve 
 * - Bt has A' (aka previous Bt) at input and the solution of A = L'*B transposed at output */
void matrixSubstitutionLower(const int nz, const real_t *L, real_t *Bt)
{
    int i,j,k;
    
    for ( j = nz-1 ; j >= 0; j-- )
    {
        for ( i = 0; i < BS*NU+NX; i++ )
        {
            for ( k = nz-1; k > j; k-- )
            {
                Bt[j*(BS*NU+NX)+i] -= Bt[k*(BS*NU+NX)+i]*L[j*nz+k];
            }
            Bt[j*(BS*NU+NX)+i] /= L[j*nz+j];
        }
    }
    
}


/* This function forms the diagonal and off-diagonal blocks of the 
 * dual Hessian based on the coupling and elimination matrices */
void form_dual_hessian( )
{
    timer tmr;

    int i;
    
    /* Build diagonal blocks Md: W_{k} depends on C_{k-1}, P_{k-1}, and P_{k} */
    for(i=0;i<NB-1;i++)
    {
        tic( &tmr );
        calculateDiagonalBlock(&data.Md[i*NX*NX], &data.CP[i*NX*(BS*NU+NX)], &data.C[i*NX*(BS*NU+NX)], &data.P_NEW[i*(BS*NU+NX)*(BS*NU+NX)]);
        data.formTimes[i] = toc( &tmr );
    }
    /* NOTE: cpu time of last operation in calculateDiagonalBlock should be 
     * part of the timings in the next block. However, the operation
     * E_{k}*P_{k}*E_{k}' is simply selecting upper block P_(1,1) and 
     * therefore its contribution is negligible */
    
    /* Build off diagonal blocks Ml: U_{k} depends on P_{k}, C_{k} */
    for(i=0;i<NB-2;i++)
    {
        tic( &tmr );
        couplingTimesEliminationTruncated(&data.Ml[i*NX*NX], &data.C[(i+1)*NX*(BS*NU+NX)], &data.P_NEW[(i+1)*(BS*NU+NX)*(BS*NU+NX)]);
        data.formTimes[i+1] += toc( &tmr );
    }
    
}


/* Helper function for form_dual_hessian */
void calculateDiagonalBlock(real_t *Md, real_t *CP, const real_t *C, const real_t *P) 
{
    int i,j;
 
    /* Result C*P (its diagonal blocks) will be reused */
    
    /* C_{k-1}*P_{k-1} */
    matMatMult(&CP[0], &C[0], &P[0], NX, (BS*NU+NX), (BS*NU+NX), 0, 0); 
        
    /* C_{k-1}*P_{k-1}*C_{k-1}' */
    matMatMult(&Md[0], &CP[0], &C[0], NX, (BS*NU+NX), NX, 0, 1);
        
    /* C_{k-1}*P_{k-1}*C_{k-1}' + E_{k}*P_{k}*E_{k}' */
    for(j=0;j<NX;j++)
    {
        for(i=0;i<NX;i++)   
        {
            /* adding upper block of next P */
            Md[j*NX+i] += P[(BS*NU+NX)*(BS*NU+NX)+j*(BS*NU+NX)+i];
        }
    }
}


/* Helper function for form_dual_hessian:
 * Calculates only the first NX columns of the multiplication C_{k}*P_{k} */
void couplingTimesEliminationTruncated(real_t *res, const real_t *C, const real_t *P)
{
    int i,j,k;

    for(k=0;k<NX;k++) /* instead of going until BS*NU+NX */
    {
        for(j=0;j<BS*NU+NX;j++)
        {
            for(i=0;i<NX;i++)
            {
                res[k*NX+i] += C[j*NX+i]*P[k*(BS*NU+NX)+j];
            }
        }
    }
}


/* Multiply coupling matrix C with a vector x 
 * C: (NB-1)*NX x NB*(BS*NU+NX) 
 * x: NB*(BS*NU+NX) x 1 */
void CtimesVector(real_t *res, const real_t *x )
{
    int i,j;
    
    for(i=0;i<NB-1;i++)
    {
        matVecMult_NX(&res[i*NX], &data.C[i*NX*(BS*NU+NX)], &x[i*(BS*NU+NX)]);
        for(j=0;j<NX;j++)
        {
            res[i*NX+j] += x[(i+1)*(BS*NU+NX)+j];
        }
    }
}


/* Multiply transpose of coupling matrix C with a vector x 
 * C': NB*(BS*NU+NX) x (NB-1)*NX 
 * x:  (NB-1)*NX  x 1 */
void CTtimesVector(real_t *res, const real_t *x )
{
    int i,j;
    
    matTransVecMult_NX(&res[0], &data.C[0], &x[0]);

    for(i=1;i<NB-1;i++)
    {
        matTransVecMult_NX(&res[i*(BS*NU+NX)], &data.C[i*NX*(BS*NU+NX)], &x[i*NX]);
        for(j=0;j<NX;j++)
        {
            res[i*(BS*NU+NX)+j] += x[(i-1)*NX+j];
        }
    }
    
    for(j=0;j<NX;j++)
    {
        res[(NB-1)*(BS*NU+NX)+j] += x[(NB-2)*NX+j];
    }
}


/* This function calculates the right hand side of the linear system */
void form_rhs( ) 
{   
    int i,j;
    int totact = 0;
        
    /* Calculate vector: RDRd - P*(H*RDRd + g) */
    
    for(i=0;i<NB;i++)
    {
        if ( data.nact[i] > 0 )
        {
            form_RDRd(&data.RZ[i*(BS*NU+NX)*(BS*NU+NX)], &data.DT[totact*(BS*NU+NX)], &data.d[totact], data.nact[i], i);
        }
        matVecMult_NV(&data.HRDRdg[i*(BS*NU+NX)], &data.H[i*(BS*NU+NX)*(BS*NU+NX)], &data.RDRd_NEW[i*(BS*NU+NX)]);
        for(j=0;j<BS*NU+NX;j++)
        {
            data.HRDRdg[i*(BS*NU+NX)+j] += data.g[i*(BS*NU+NX)+j];
        }        
        matVecMult_NV(&data.PHRDRdg[i*(BS*NU+NX)], &data.P_NEW[i*(BS*NU+NX)*(BS*NU+NX)], &data.HRDRdg[i*(BS*NU+NX)]);
        
        totact += data.nact[i];
    }
    
    for(i=0;i<NB*(BS*NU+NX); i++)
    {
        data.RHS_temp[i] = data.RDRd_NEW[i] - data.PHRDRdg[i];
    }
   
    /* multiply RHS_temp with coupling matrix C */
    CtimesVector(data.RHS, data.RHS_temp);
    
    /* subtract c */
    for(i=0;i<(NB-1)*NX;i++)
    {
        data.RHS[i] -= data.c[i];
    }
        
}


/* Helper function for form_rhs:
 * Factorize the term DR and then calculate R*((D*R)^-1)*d  */
void form_RDRd(const real_t *R, const real_t *Dt, real_t *d, const int dim, const int nblock)
{
    int info;
    int i;
    real_t maxerr = 0;
    real_t err;

    /* D*R */
    matMatMult(&data.DR[nblock*(BS*NU+NX)*(BS*NU+NX)], Dt, R, dim, BS*NU+NX, dim, 1, 0);

    /* positive definite (D*R)'*(D*R) */
    matMatMult(&data.DRDR[nblock*(BS*NU+NX)*(BS*NU+NX)], &data.DR[nblock*(BS*NU+NX)*(BS*NU+NX)], &data.DR[nblock*(BS*NU+NX)*(BS*NU+NX)], dim, dim, dim, 1, 0);

    /* Adjust d <-- (D*R)'*d */
    matVecMult(&data.DRd[nblock*(BS*NU+NX)], &data.DR[nblock*(BS*NU+NX)*(BS*NU+NX)], d, dim, dim, 1);     

    /* factorization */
    choleskyFactorization(&dim, &data.DRDR[nblock*(BS*NU+NX)*(BS*NU+NX)], &dim, &info);

    /* two solves to get inv(DR)*d */
    forwardSubstitution(dim, &data.DRDR[nblock*(BS*NU+NX)*(BS*NU+NX)], &data.DRd[nblock*(BS*NU+NX)]);
    backwardSubstitution(dim, &data.DRDR[nblock*(BS*NU+NX)*(BS*NU+NX)], &data.DRd[nblock*(BS*NU+NX)]);
    
    /* multiply with R */
    matVecMult(&data.RDRd_NEW[nblock*(BS*NU+NX)], R, &data.DRd[nblock*(BS*NU+NX)], BS*NU+NX, dim, 0);

}


/* Helper function for from_RDRd:
 * Vector forward substitution L*X = V (solve with respect to X) 
 * - Equivalent to X'*R = V' where elements in R are accessed by column 
 * - Results overwrites v 
 * - Here we read L instead of R (less efficient) to avoid transposing the Cholesky factor */
void forwardSubstitution(const int dim, const real_t *L, real_t *v)
{
    int i,j;
    
    for ( j = 0; j < dim; j++ )
    {
        for ( i = 0; i < j; i++ )
        {
            v[j] -=  v[i]*L[i*dim+j];
        }
        v[j] /= L[j*dim+j]; 
    }
}


/* Helper function for from_RDRd:
 * Vector backward substitution R*X = V (solve with respect to X) 
 * - Equivalent to X'*L = V' where elements in L are accessed by column 
 * - Results overwrites v */
void backwardSubstitution(const int dim, const real_t *L, real_t *v)
{
    int i,j;
    
    for ( j = dim-1; j >= 0; j-- )
    {
        for ( i = dim-1; i > j; i-- )
        {
            v[j] -=  v[i]*L[j*dim+i];
        }
        v[j] /= L[j*dim+j]; 
    }
}


/* This function solves the linear system, using a block cholesky factorization */
void solve_dual_system( ) 
{
    blockCholeskyFactorization(data.Md, data.Ml);
    blockForwardSubstitution(data.Md, data.Ml, data.RHS);
    blockBackwardSubstitution(data.Md, data.Ml, data.RHS);
    
}


/* This function calculates the primal step, given the current multipliers */
void update_primal( ) 
{
    int i;
    
    CTtimesVector(data.primalStep, data.RHS); /* lambdas are stored in RHS at this point */
  
    for(i=0;i<NB*(BS*NU+NX);i++)
    {
        data.primalStep[i] += data.HRDRdg[i];
    }
    
    
    for(i=0;i<NB;i++)
    {
        matVecMult_NV(&data.primal_temp[i*(BS*NU+NX)], &data.P_NEW[i*(BS*NU+NX)*(BS*NU+NX)], &data.primalStep[i*(BS*NU+NX)]);
    }
       
    for(i=0;i<NB*(BS*NU+NX);i++)
    {
        data.primalStep[i] = data.RDRd_NEW[i] - data.primal_temp[i];
    }
 
}


void blockCholeskyFactorization(real_t *Md, real_t *Ml)
{
    int nblk = NX;
    int info;
    
    int i,n;
    
    /* first block */
    choleskyFactorization(&nblk,&Md[0],&nblk,&info);
    
    for(n=1;n<=NB-2;n++)
    {
        /* Matrix substitution */
        myblas_dtrsm(NX, NX, &Md[(n-1)*NX*NX], &Ml[(n-1)*NX*NX]);
        
        /* Substitution of diagonal block
         * --> this operation may also add invalid data to upper triangular
         * part that should NOT be referenced */
        myblas_dgemm(NX, NX, NX, &Ml[(n-1)*NX*NX], &Ml[(n-1)*NX*NX], &Md[n*NX*NX]);

        /* Cholesky factorization */
        choleskyFactorization(&nblk,&Md[n*NX*NX],&nblk,&info);
    }
    
}


void blockForwardSubstitution(real_t *Md, real_t *Ml, real_t *x)
{
    int row, col, n;
    
    /* Substitute first diagonal block */
    for (row = 0; row < NX; row++)
    {
        for (col = 0; col <= row-1; col++)
        {
            x[row] -= Md[row + col*NX] * x[col];
        }
        x[row] = x[row] / Md[row + row*NX];
    }
    
    /* Substitute remaining diagonal and off-diagonal blocks */
    for (n = 1; n <= NB-2; n++)
    {
        for (row = 0; row < NX; row++)
        {
            for (col = 0; col < NX; col++)
            {
                x[n*NX + row] -= Ml[(n-1)*NX*NX + row + col*NX] * x[(n-1)*NX + col];
            }
            for (col = 0; col <= row-1; col++)
            {
                x[n*NX + row] -= Md[n*NX*NX + row + col*NX] * x[n*NX + col];
            }
            x[n*NX + row] = x[n*NX + row] / Md[n*NX*NX + row + row*NX];
        }
    }
}


void blockBackwardSubstitution(real_t *Md, real_t *Ml, real_t *x)
{
    int row, col, n;
    
    /* Substitute last diagonal block */
    for (row = NX-1; row >= 0; row--)
    {
        for (col = NX-1; col > row; col--)
        {
            x[(NB-2)*NX + row] -= Md[(NB-2)*NX*NX + col + row*NX] * x[(NB-2)*NX + col];
        }
        x[(NB-2)*NX + row] = x[(NB-2)*NX + row] / Md[(NB-2)*NX*NX + row + row*NX];
    }
    
    /* Substitute remaining diagonal and off-diagonal blocks */
    for (n = NB-3; n >= 0; n--)
    {
        for (row = NX-1; row >= 0; row--)
        {
            for (col = NX-1; col >= 0; col--)
            {
                x[n*NX + row] -= Ml[(n)*NX*NX + col + row*NX] * x[(n+1)*NX + col];
            }
            for (col = NX-1; col > row; col--)
            {
                x[n*NX + row] -= Md[n*NX*NX + col + row*NX] * x[n*NX + col];
            }
            x[n*NX + row] = x[n*NX + row] / Md[n*NX*NX + row + row*NX];
        }
    }
}


/* Matrix substitution */
void myblas_dtrsm(const int M, const int N, const real_t *A, real_t *B)
{
    /* NOTE: this is tailored to our application only */
    
    /* cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, M, N, 1.0, A, NX, B, NX); */
    
    int row, col, n;
    
    for (n = 0; n < N; n++)
    {
        for (col = 0; col < M; col++)
        {
            for (row = 0; row < col; row++)
            {
                B[col*N + n] -= A[col + row*M] * B[row*N + n];
            }
            B[col*N + n] = B[col*N + n] / A[col + col*M];
        }
        
    }
    
}


/* Tailored matrix matrix multiplication */
void myblas_dgemm(const int M, const int N, const int K, const real_t *A,
        const real_t *B, real_t *C)
{
    
    /* NOTE: this is tailored to our application only */
    
    int i, j, k;
    
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < K; k++)
                C[i + N*j] -= A[i + k*M] * B[k*K + j];
    
}


/* Dense Cholesky factorization (code adapted from qpOASES)
 *
 *  License:
 *
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
void choleskyFactorization(const int *_n, real_t *a, const int *_lda, int *info)
{
    
    /* char low = 'L';
     * dpotrf_ (&low, _n, a, _lda, info); */
    
    real_t sum;
    int i, j, k;
    int n = (int)(*_n);
    int lda = (int)(*_lda);
    real_t tmp;
    
    for( i=0; i<n; ++i )
    {
        /* j == i */
        sum = a[i + lda*i];
        
        for( k=(i-1); k>=0; --k )
            sum -= a[i+lda*k] * a[i+lda*k];
        
        if ( sum > 0.0 )
            a[i+lda*i] = sqrt( sum );
        else
        {
            a[0] = sum; /* tunnel negative diagonal element to caller */
            if (info != 0)
                *info = (int)i+1;
            return;
        }
        
        for( j=(i+1); j<n; ++j )
        {
            sum = a[i*lda + j];
            
            for( k=(i-1); k>=0; --k )
                sum -= a[i+lda*k] * a[j+lda*k];
            
            a[j+lda*i] = sum / a[i+lda*i];
        }
    }
    if (info != 0)
        *info = 0;
}