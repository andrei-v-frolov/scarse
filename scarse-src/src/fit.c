/* $Id: fit.c,v 1.2 2001/02/05 01:21:43 frolov Exp $ */

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

#include <stdlib.h>
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



/******************* Simulated annealing ******************************/

/* n-dimensional amoeba matrix layout:
 *
 *                    |   [0..n-1]   |  [n]  
 *             -------+--------------+-------
 *             [0..n] |simplex coords| f(x)  
 *             -------+--------------+-------
 *  S[i][j] =   [n+1] |best pt coords| f(b)  
 *              [n+2] |vertex sum    |   0   
 *              [n+3] |tmp           |   .   
 *             -------+--------------+-------
 */

/* Restart amoeba around best point */
void restart_amoeba(double **S, int n, double (*func)(double []), double lambda)
{
	int i, j;
	double *best = S[n+1], *sum = S[n+2];
	
	for (i = 0; i <= n; i++)
		for (j = 0; j <= n; j++)
			S[i][j] = best[j];
	
	for (i = 0; i < n; i++) {
		S[i][i] += lambda;
		S[i][n] = (*func)(S[i]);
	}
	S[n][n] = (*func)(S[n]);
	
	for (i = 0; i < n; i++)
		sum[i] = (n+1) * best[j] + lambda;
	sum[n] = 0.0;
}

/* Create new amoeba */
double **new_amoeba(double x[], int n, double (*func)(double []), double lambda)
{
	int i;
	double **S = matrix(n+4, n+1), *best = S[n+1];
	
	for (i = 0; i < n; i++)
		best[i] = x[i];
	best[n] = HUGE;
	
	restart_amoeba(S, n, func, lambda);
	
	return S;
}


/* Try to crawl in given direction */
static double crawl(double **S, int n, double (*func)(double []), double T, int dir, double amount)
{
	int i;
	double y, *x = S[dir], *best = S[n+1], *sum = S[n+2], *try = S[n+3];
	
	/* Give it a try... */
	for (i = 0; i < n; i++)
		try[i] = (1.0-amount) * (sum[i] - x[i])/n + amount * x[i];
	try[n] = (*func)(try);
	
	/* Best move ever? */
	if (try[n] < best[n])
		for (i = 0; i <= n; i++) best[i] = try[i];
	
	y = try[n] - T * log(RAND_MAX/rand());
	
	/* Favourable move? */
	if (y < x[n]) {
		for (i = 0; i < n; i++) {
			sum[i] += try[i] - x[i];
			x[i] = try[i];
		}
		
		x[n] = try[n];
	}
	
	return y;
}

/* Shrink the amoeba around given vertex */
static void shrink(double **S, int n, double (*func)(double []), double T, int dir, double amount)
{
	int i, j;
	double *x = S[dir], *best = S[n+1], *sum = S[n+2];
	
	/* Shrink the amoeba */
	for (i = 0; i <= n; i++) if (i != dir) {
		for (j = 0; j < n; j++)
			S[i][j] = (1.0-amount) * x[j] + amount * S[i][j];
		S[i][n] = (*func)(S[i]);
		
		if (S[i][n] < best[n])
			for (j = 0; j <= n; j++) best[j] = S[i][j];
	}
	
	/* Update vertex sum */
	for (i = 0; i < n; i++) {
		sum[i] = 0.0;
		
		for (j = 0; j <= n; j++)
			sum[i] += S[j][i];
	}
}

/* Minimize the function by simulated annealing */
void anneal(double **S, int n, double (*func)(double []), double T0, int maxsteps, double tol)
{
	int i, j, k, lo, hi;
	double T, y, ylo, yhi, yhi2;
	
	for (k = 0; k < maxsteps; ) {
		/* Cooling schedule */
		T = T0 * exp(-log(100.0)*k/maxsteps);
		
		/* Rank simplex vertices */
		lo = hi = 0;
		ylo = HUGE; yhi = yhi2 = -HUGE;
		
		for (i = 0; i <= n; i++) {
			y = S[i][n] + T * log(RAND_MAX/rand());
			
			if (y < ylo) { lo = i; ylo = y; }
			if (y > yhi) { yhi2 = yhi; hi = i; yhi = y; }
			else if (y > yhi2) { yhi2 = y; }
		}
		
		/* Are we done yet? */
		if (2.0*fabs(S[hi][n]-S[lo][n])/(fabs(S[hi][n])+fabs(S[lo][n])) < tol) break;
		
		/* Make a move: try reflect first */
		y = crawl(S, n, func, T, hi, -1.0); k++;
		
		if (y <= ylo) {
			/* was good, try expanding */
			y = crawl(S, n, func, T, hi, 2.0); k++;
		} else if (y >= yhi2 ) {
			/* no good, try contracting */
			y = crawl(S, n, func, T, hi, 0.5); k++;
			
			/* if that didn't work, try shrinking */
			if (y >= yhi2) { shrink(S, n, func, T, lo, 0.5); k += n; }
		}
	}
		printf("\nT = %g\n", T);
		for (i = 0; i <= n+2; i++) {
			for (j = 0; j <= n; j++) printf("%12.8g\t", S[i][j]);
			printf("\n");
		}
		
		
}



/******************* Non-linear curve model ***************************/

/* Curve data is passed as global variable */
static int _curve_pts_;
static double **_curve_data_;


/* Curve model is power law with possible saturation */
static double model(double p[], double x)
{
	double y;
	
	y = p[0]*pow(x, p[2]) + p[1];
	
	if (y < p[3]) y = p[3];
	if (y > p[4]) y = p[4];
	
	return y;
}

/* Minimization criterion is least square */
static double chi2(double p[])
{
	double t, s = 0.0;
	int i, n = _curve_pts_;
	double *x = _curve_data_[0], *y = _curve_data_[1];
	
	for (i = 0; i < n; i++) {
		t = y[i] - model(p, x[i]); s += t*t;
	}
	
	/* Penalty for saturation points wandering outside range */
	if (p[3] < 0.0) s += 1.0e6 * fabs(p[3]);
	if (p[4] > 1.0) s += 1.0e6 * fabs(p[4]-1.0);
	
	return s/n;
}

/* Fit parametric curve model y=f(x) to the data using simulated annealing */
double *fit_curve(double **data, int n)
{
	int i;
	#define D 5
	double *p = vector(D+1), **S;
	
	_curve_data_ = data; _curve_pts_ = n;
	p[0] = p[2] = p[4] = 1.0; p[1] = p[3] = 0.0;
	
	S = new_amoeba(p, D, chi2, 1.0); anneal(S, D, chi2, 1.0, 1000, 1.0e-6);
	restart_amoeba(S, D, chi2, 1.0); anneal(S, D, chi2, 0.1, 1000, 1.0e-6);
	restart_amoeba(S, D, chi2, 1.0); anneal(S, D, chi2, 0.0, 1000, 1.0e-6);
	
	for (i = 0; i < D; i++)
		p[i] = S[D+1][i];
	
	free_matrix(S);
	return p;
	
	#undef D
}


/* Check if the value is within curve range */
int within_range(double p[], double y)
{
	return (y > p[3]) && (y < p[4]);
}

/* Curve lookup x = f^(-1)(y) */
double curve(double p[], double y)
{
	register double x = (y - p[1])/p[0];
	
	return ppow(x, 1.0/p[2]);
}

/* Reverse lookup y = f(x) */
double curve_1(double p[], double x)
{
	return p[0]*pow(x, p[2]) + p[1];
}

/* Curve slope */
double curve_dydx(double p[], double x)
{
	return p[0]*p[2]*pow(x, p[2]-1.0);
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
