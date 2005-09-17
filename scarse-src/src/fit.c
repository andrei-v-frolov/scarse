/* $Id: fit.c,v 1.6 2005/09/17 03:13:38 afrolov Exp $ */

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



/******************* Multidimensional minimization ********************/

/* directional minimization along a vector xi in n dimensions */
double vmin(double p[], double xi[], int n, double (*f)(double []), double eps)
{
	double a, b, c, u, v, w, x, fa, fb, fc, fu, fv, fw, fx;
	double q, r, s, t, tol, e = 0.0, d = 0.0;
	
	int i, maxiter = 84; // maximal number of iterations
	
	#define GOLD  1.61803398874989484820458683436563811772030918
	#define CGOLD 0.38196601125010515179541316563436188227969082
	
	#define EVAL(X,F) { double t[n]; for (i = 0; i < n; i++) t[i] = p[i] + (X)*xi[i]; (F) = (*f)(t); }
	#define SWAP(A,B) { double T = (A); (A) = (B); (B) = T; } 
	#define SIGN(A,B) ((B) >= 0.0 ? fabs(A) : -fabs(A))
	
	/* initial bracketing of a minimum */
	a = 0.0; b = 1.0; EVAL(a,fa); EVAL(b,fb);
	if (fb > fa) { SWAP(a,b); SWAP(fa,fb); }
	c = b + GOLD*(b-a); EVAL(c,fc);
	
	while (fb > fc) {
		a = b; b = c; fa = fb; fb = fc;
		c = b + GOLD*(b-a); EVAL(c,fc);
	}
	
	/* Brent's minimization */
	x = w = v = b; fx = fw = fv = fb; if (a > c) SWAP(a,c);
	
	while (maxiter--) {
		b = (a+c)/2.0; tol = eps * sqrt(1.0 + x*x);
		if (fabs(x-b) <= (2.0*tol - (c-a)/2.0)) break;
		
		if (fabs(e) > tol) {
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			s = (x-v)*q-(x-w)*r;
			
			q = 2.0*(q-r); if (q > 0.0) s = -s; else q = -q;
			
			t = e; e = d;
			
			if (fabs(s) >= fabs(0.5*q*t) || s <= q*(a-x) || s >= q*(c-x)) {
				e = (x >= b ? a-x : c-x); d = CGOLD * e;
			} else {
				d = s/q; u = x+d; if (u-a < 2.0*tol || c-u < 2.0*tol) d = SIGN(tol,b-x);
			}
		} else { e = (x >= b ? a-x : c-x); d = CGOLD * e; }
		
		u = (fabs(d) >= tol ? x+d : x+SIGN(tol,d)); EVAL(u,fu);
				
		if (fu <= fx) {
			if (u >= x) a = x; else c = x;
			v = w; w = x; x = u;
			fv = fw; fw = fx; fx = fu;
		} else {
			if (u < x) a = u; else c = u;
			if (fu <= fw || w == x) { v = w; w = u; fv = fw; fw = fu; }
			else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
		}
	}
	
	/* update direction vectors */
	if (x != 0.0) for (i = 0; i < n; i++) { xi[i] *= x; p[i] = p[i] + xi[i]; }
	
	return fx;
}

/* n-dimensional minimization (using Powell's method) */
double mmin(double p[], double **e, int n, double (*f)(double []), double eps)
{
	int i, j, s; double fp, fo, fx, sd, t; double po[n], px[n], xi[n];
	
	int maxiter = 1000; // maximal number of iterations
	
	if (n == 1) { xi[0] = e[0][0]; return vmin(p, xi, n, f, eps); }
	
	fp = (*f)(p);
	
	while (maxiter--) {
		for (j = 0; j < n; j++) po[j] = p[j]; fo = fp; s = 0; sd = 0.0;
		
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) xi[j] = e[j][i];
			
			t = fp; fp = vmin(p, xi, n, f, eps);
			if (fabs(fp-t) > sd) { sd = fabs(fp-t); s = i; }
		}
		
		if (2.0*fabs(fo-fp) <= eps*eps*(fabs(fo)+fabs(fp))) return fp;
		
		for (j = 0; j < n; j++) { px[j] = 2.0*p[j]-po[j]; xi[j] = p[j]-po[j]; }
		
		fx = (*f)(px);
		
		if (fx < fo && 2.0*(fo-2.0*fp+fx)*(fo-fp-sd)*(fo-fp-sd) < sd*(fo-fx)*(fo-fx)) {
			fp = vmin(p, xi, n, f, eps);
			
			for (j = 0; j < n; j++) {
				for (i = s; i < n-1; i++) e[j][i] = e[j][i+1]; e[j][n-1] = xi[j];
			}
		}
	}
	
	return fp;
}



/******************* Non-linear curve model ***************************/

/* Curve data is passed as global variable */
static int _curve_pts_;
static double **_curve_data_;


/* Curve model is power law with possible clipping */
static double curve_model(double p[], double y)
{
	double x;
	
	x = p[0]*pow(y, p[2]) + p[1];
	
	if (x < p[3]) x = p[3];
	if (x > p[4]) x = p[4];
	
	return x;
}

/* Minimization criterion is least square */
static double curve_chi2(double p[])
{
	int i, n = _curve_pts_;
	double t, s = 0.0, min = 1.0, max = 0.0;
	double *y = _curve_data_[0], *x = _curve_data_[1], *dx = _curve_data_[2];
	
	for (i = 0; i < n; i++) {
		if (x[i] < min) min = x[i];
		if (x[i] > max) max = x[i];
		
		t = (x[i] - curve_model(p, y[i]))/dx[i]; s += t*t;
	}
	
	/* Penalty for clipping bounds wandering outside range */
	if (p[3] < min) s *= 1.0 + 100.0*(p[3]-min)*(p[3]-min);
	if (p[4] > max) s *= 1.0 + 100.0*(p[4]-max)*(p[4]-max);
	
	return s/n;
}

/* Fit parametric curve model y=f(x) to the data using simulated annealing */
double *fit_curve(double **data, int n)
{
	int i, j;
	#define D 5
	double *p = vector(D+1), **e = matrix(D,D);
	
	/* initialize curve data and model */
	_curve_data_ = data; _curve_pts_ = n;
	p[0] = p[2] = p[4] = 1.0; p[1] = p[3] = 0.0;
	
	/* initialize direction set */
	for (i = 0; i < D; i++) { for (j = 0; j < D; j++) e[i][j] = 0.0; e[i][i] = 1.0; }
	
	/* optimize */
	p[D] = mmin(p, e, D, curve_chi2, 1.0e-6);
	
	free_matrix(e); return p;
	
	#undef D
}


/* Curve lookup y = f(x) */
double lu_curve(double p[], double x)
{
	register double y = (x - p[1])/p[0];
	
	return ppow(y, 1.0/p[2]);
}

/* Reverse lookup x = f^(-1)(y) */
double lu_curve_1(double p[], double y)
{
	return p[0]*pow(y, p[2]) + p[1];
}

/* Curve slope */
double curve_dydx(double p[], double y)
{
	return 1.0/(p[0]*p[2]*pow(y, p[2]-1.0));
}


/******************* Robust matrix fit estimate ***********************/
#warning Done up to here...

/* LUT data is passed as global variable */
static int _lut_pts_;
static int _lut_channel_;
static double **_lut_data_;


/* LUT model is linear fit with possible clipping */
static double lut_model(double p[], double XYZ[])
{
	double x;
	
	x = p[0]*XYZ[0] + p[1]*XYZ[1] + p[2]*XYZ[2];
	
	if (x < p[3]) x = p[3];
	if (x > p[4]) x = p[4];
	
	return x;
}

/* Minimization criterion is least square */
static double lut_chi2(double p[])
{
	int i, n = _lut_pts_;
	double t, s = 0.0, min = 1.0, max = 0.0;
	
	for (i = 0; i < n; i++) {
		double *XYZ = _lut_data_[i], x = _lut_data_[i][_lut_channel_+3], dx = _lut_data_[i][_lut_channel_+6];
		
		if (x < min) min = x;
		if (x > max) max = x;
		
		t = (x - lut_model(p, XYZ))/dx; s += t*t;
	}
	
	/* Penalty for clipping bounds wandering outside range */
	if (p[3] < min) s *= 1.0 + 100.0*(p[3]-min)*(p[3]-min);
	if (p[4] > max) s *= 1.0 + 100.0*(p[4]-max)*(p[4]-max);
	
	return s/n;
}

/* Fit linear model RGB=M*XYZ to the data using simulated annealing */
void robust_linear_fit(double **data, int n, double M[3][3])
{
	int i, j, k;
	#define D 5
	
	for (k = 0; k < 3; k++) {
		double *p = vector(D+1), **e = matrix(D,D);
		
		/* initialize LUT data and model */
		_lut_data_ = data; _lut_pts_ = n; _lut_channel_ = k;
		p[0] = p[1] = p[2] = p[3] = 0.0; p[i] = p[4] = 1.0;
		
		/* initialize direction set */
		for (i = 0; i < D; i++) { for (j = 0; j < D; j++) e[i][j] = 0.0; e[i][i] = 1.0; }
		
		/* optimize */
		p[D] = mmin(p, e, D, lut_chi2, 1.0e-6);
		
		for (j = 0; j < 3; j++) M[k][j] = p[j];
		for (j = 0; j <= D; j++) fprintf(stderr, "%12.10g\t", p[j]);
		fprintf(stderr, "\n");
		
		free_matrix(e); free_vector(p);
	}
	
	#undef D
}



/**********************************************************************/

/*
 * Linear model fitting:
 *   chi^2 optimization via simulated annealing
 * 
 * Data models:
 *   - power-law curve with saturation
 *   - 3x3 matrix transform with saturation
 * 
 */


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
