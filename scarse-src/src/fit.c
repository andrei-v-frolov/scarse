/* $Id: fit.c,v 1.1 2001/01/26 22:45:32 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Data fitting and approximation routines.
 * 
 * Copyright (C) 2000 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include <math.h>

#include "util.h"
#include "spaces.h"



/******************* Some linear algebra first ************************/

/* Gauss-Jordan elimination with full pivoting */
/* returns inverse of the matrix in A, and solution to Ax=B in B */
void gaussj(double **A, int n, double **B, int m)
{
	double pivot, Aji;
	int i, j, k, pr, pc, temp;
	double **U = matrix(n, n+m);
	int *c = ivector(n), *r = ivector(n);
	
	#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
	
	/* Initialize */
	for (i = 0; i < n; i++) {
		c[i] = r[i] = i;
		for (j = 0; j < n; j++) U[i][j] = 0.0; U[i][i] = 1.0;
		for (j = 0; j < m; j++) U[i][n+j] = B[i][j];
	}
	
	/* Gauss-Jordan elimination loop */
	for (i = 0; i < n; i++) {
		/* Find next pivot element */
		pivot = 0.0; pr = pc = i;
		
		for (j = i; j < n; j++) for (k = i; k < n; k++) {
			if (fabs(A[r[j]][c[k]]) > fabs(pivot)) {
				pivot = A[r[j]][c[k]];
				pr = j; pc = k;
			}
		}
		
		if (pivot == 0.0) error("Singular matrix in gaussj()");
		
		/* Bring it into top left corner */
		SWAP(r[i], r[pr]); SWAP(c[i], c[pc]);
		
		/* Eliminate all elements in a column except pivot one */
		A[r[i]][c[i]] = 1.0;
		for (j = i+1; j < n; j++) A[r[i]][c[j]] /= pivot;
		for (j = 0; j < n+m; j++) U[r[i]][j] /= pivot;
		
		for (j = 0; j < n; j++) if (j != i) {
			Aji = A[r[j]][c[i]]; A[r[j]][c[i]] = 0.0;
			for (k = i+1; k < n; k++) A[r[j]][c[k]] -= Aji*A[r[i]][c[k]];
			for (k = 0; k < n+m; k++) U[r[j]][k] -= Aji*U[r[i]][k];
		}
		
	}
	
	/* recover original order of the rows and columns */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) A[c[i]][j] = U[r[i]][j];
		for (j = 0; j < m; j++) B[c[i]][j] = U[r[i]][n+j];
	}
	
	free_matrix(U);
	free_vector(c);
	free_vector(r);
	
	#undef SWAP
}



/******************* Generalized least square fit *********************/

/* Evaluate approximation y = A * F(x) */
void approx(double **A, void (*basis)(double [], double []), int d, double x[], double y[])
{
	int i, j;
	double *F = vector(d);
	
	(*basis)(x, F);
	
	for (i = 0; i < 3; i++) {
		double *a = A[i], P = 0.0;
		
		for (j = d-1; j >= 0; j--) P += a[j]*F[j];
		
		y[i] = P;
	}
	
	free_vector(F);
}

/* Find best fit for x[3..5] = A * F(x[0..2]) wrt perceptual Lab metric */
double **best_fit(double **x, int n, void (*basis)(double [], double []), int d)
{
	int i, j, k, N = 3*d;
	double xk[3], yk[3], g[3][3];
	double *F = vector(d), **A = matrix(3,d);
	double **S = matrix(N,N), **SY = matrix(N,1);
	
	/* Zero covariance matrix */
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) S[i][j] = 0.0;
		
		SY[i][0] = 0.0;
	}
	
	/* Accumulate data */
	for (k = 0; k < n; k++) {
		for (i = 0; i < 3; i++) {
			xk[i] = x[i][k];
			yk[i] = x[i+3][k];
		}
		
		(*basis)(xk, F);
		XYZ_metric(yk, g);
		
		for (i = 0; i < N; i++) {
			int a = i/3, alpha = i%3;
			
			for (j = 0; j < N; j++) {
				int b = j/3, beta = j%3;
				
				S[i][j] += g[alpha][beta]*F[a]*F[b];
			}
			
			SY[i][0] += (g[alpha][0]*yk[0]+g[alpha][1]*yk[1]+g[alpha][2]*yk[2])*F[a];
		}
	}
	
	/* Find coefficients */
	gaussj(S, N, SY, 1);
	
	for (i = 0; i < N; i++) {
		int a = i/3, alpha = i%3;
		
		A[alpha][a] = SY[i][0];
	}
	
	/* Cleanup and exit */
	free_vector(F);
	free_matrix(S);
	free_matrix(SY);
	
	return A;
}


/******************* Wrapper: Linear fit ******************************/

static void linear(double x[], double F[])
{
	F[0] = x[0]; F[1] = x[1]; F[2] = x[2];
}

void best_linear_fit(double **x, int n, double M[3][3])
{
	int i, j;
	double **A = best_fit(x, n, linear, 3);
	
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) M[i][j] = A[i][j];
	
	free_matrix(A);
}


/******************* Wrapper: Polynomial fit **************************/

static void powers(double x[], double F[])
{
	/* x^1 */
	F[ 0] = x[0];
	F[ 1] = x[1];
	F[ 2] = x[2];
	
	/* x^2 */
	F[ 3] = x[0]*F[0];
	F[ 4] = x[0]*F[1];
	F[ 5] = x[0]*F[2];
	F[ 6] = x[1]*F[1];
	F[ 7] = x[1]*F[2];
	F[ 8] = x[2]*F[2];
	
	/* x^3 */
	F[ 9] = x[0]*F[3];
	F[10] = x[0]*F[4];
	F[11] = x[0]*F[5];
	F[12] = x[0]*F[6];
	F[13] = x[0]*F[7];
	F[14] = x[0]*F[8];
	F[15] = x[1]*F[6];
	F[16] = x[1]*F[7];
	F[17] = x[1]*F[8];
	F[18] = x[2]*F[8];
}

void poly_approx(double **A, double x[], double y[])
{
	approx(A, powers, 19, x, y);
}

double **best_poly_fit(double **x, int n)
{
	return best_fit(x, n, powers, 19);
}
