/* $Id: interp.c,v 1.1.1.1 2001/01/26 22:45:31 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Interpolation in one and three dimensions.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include "util.h"


/******************* R^1 -> R^1 interpolation *************************/

/* Matrix M[0..3][0..n-1] contains:
 *   M[0] = x	-  coordinates of points
 *   M[1] = z	-  values of function in these points
 *   M[2] = z'	-  first derivatives of function in these ponts
 *   M[3] = z''	-  second derivatives of function in these ponts
 */


/* Calculate derivatives for smooth interpolation */
void derivs(double **M, unsigned long n)
{
	unsigned long i;
	double *u = vector(n-1);
	double *X = M[0], *Z = M[1], *A = M[2], *B = M[3];
	
	#define a(i) 2.0*(X[i+1]-X[i-1])
	#define h(i) (X[i+1]-X[i])
	#define p(i) ((Z[i+1]-Z[i])/(X[i+1]-X[i]))
	#define l(i) (h(i-1)/u[i-1])
	
	
	sort2(n, X, Z);
	for (i = 0; i < n-1; i++)
		if (h(i) == 0.0) error("Bad input to derivs()");
	
	u[1] = a(1);
	B[0] = B[n-1] = 0.0;
	B[1] = 3.0*(p(1)-p(0));
	for (i = 2; i < n-1; i++) {
		u[i] = a(i) - l(i)*h(i);
		B[i] = 3.0*(p(i)-p(i-1)) - l(i)*B[i-1];
	}
	
	B[n-2] = B[n-2]/u[n-2];
	for (i = n-3; i > 0; i--) {
		B[i] = (B[i] - h(i+1)*B[i+1])/u[i];
		A[i] = p(i) - (2.0*B[i] + B[i+1])*h(i)/3.0;
	}
	
	A[0] = p(0) - B[1]*h(0)/3.0;
	A[n-2] = p(n-2) - 2.0/3.0*B[n-2]*h(i);
	A[n-1] = A[n-2] + B[n-2]*h(n-2);
	
	free_vector(u);
	
	#undef a
	#undef h
	#undef p
	#undef l
}


/* Interpolate function in 1d using cubic splines */
double interp1d(double **M, unsigned long n, double x)
{
	unsigned long k = 0, l = n-1, m;
	double dx, *X = M[0], *Z = M[1], *A = M[2], *B = M[3], Ck;
	
	while (l-k > 1) {
		m = (l+k) >> 1;
		if (X[m] > x) l = m; else k = m;
	}
	
	dx = x - X[k];
	Ck = (B[l]-B[k])/(X[l]-X[k])/3.0;
	
	return Z[k] + dx*(A[k] + dx*(B[k] + dx*Ck));
}



/******************* R^3 -> R^3 interpolation *************************/

/* Implementation of choice is simple linear simplex interpolation.
 * It works better than fancier higher orders methods on our data,
 * because it is immune to overshooting caused by excessive stiffness.
 */

#define XMAX 10.0
#define YMAX 10.0
#define ZMAX 10.0


typedef struct {
	int vertex[4];
	double X0[3], M[3][3];
} simplex;

typedef struct {
	int n;
	double **m;
	int ns;
	simplex *simplx;
} partition;


/* Matrix m[0..5][0..n+3] contains:
 *   m[0..2] = x  -  coordinates of points
 *   m[3..5] = z  -  value of function in these points
 */

/* Create new simplex with given vertices */
static void mksimplx(simplex *s, double **m, int v0, int v1, int v2, int v3)
{
	int i, j;
	double T[3][3];
	
	s->vertex[0] = v0;
	s->vertex[1] = v1;
	s->vertex[2] = v2;
	s->vertex[3] = v3;
	
	s->X0[0] = m[0][v0];
	s->X0[1] = m[1][v0];
	s->X0[2] = m[2][v0];
	
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
		T[i][j] = m[i][s->vertex[j+1]] - s->X0[i];
	
	inv33(T, s->M);
}

/* Is point within a simplex? */
static int within(double x[], simplex *s)
{
	double *X0 = s->X0, X[3], Y[3];
	
	X[0] = x[0] - X0[0];
	X[1] = x[1] - X0[1];
	X[2] = x[2] - X0[2];
	
	apply33(s->M, X, Y);
	
	return  (Y[0] > -TINY) && (Y[1] > -TINY) &&
		(Y[2] > -TINY) && (Y[0] + Y[1] + Y[2] < 1.0+TINY);
}

/* Linear interpolation within a simplex */
static double interp_simplex(double x0[], simplex *s, double **m, int i)
{
	double x[4] = {m[0][s->vertex[0]]-x0[0], m[0][s->vertex[1]]-x0[0], m[0][s->vertex[2]]-x0[0], m[0][s->vertex[3]]-x0[0]};
	double y[4] = {m[1][s->vertex[0]]-x0[1], m[1][s->vertex[1]]-x0[1], m[1][s->vertex[2]]-x0[1], m[1][s->vertex[3]]-x0[1]};
	double z[4] = {m[2][s->vertex[0]]-x0[2], m[2][s->vertex[1]]-x0[2], m[2][s->vertex[2]]-x0[2], m[2][s->vertex[3]]-x0[2]};
	double F[4] = {m[3+i][s->vertex[0]], m[3+i][s->vertex[1]], m[3+i][s->vertex[2]], m[3+i][s->vertex[3]]};
	
	double  t5 = y[3]*z[2], t7 = z[2]*F[3], t8 = x[1]*y[0], t13 = F[2]*z[3],
		t16 = x[3]*y[0], t18 = x[2]*F[1], t19 = z[0]*y[3], t22 = y[3]*F[2],
		t24 = x[3]*z[0], t27 = x[2]*z[1], t32 = x[3]*F[0], t44 = y[0]*z[3],
		t46 = F[1]*x[0], t52 = z[1]*x[0], t58 = y[1]*x[0], t65 = z[2]*x[0],
		t70 = y[2]*x[0], t73 = y[2]*x[1], t75 = y[2]*x[3], t77 = x[2]*y[0],
		t86 = z[2]*x[1], t94 = x[2]*z[0], t97 = x[3]*z[2];
	
	double  t35 = -F[3]*x[2]*z[0]*y[1] + x[1]*F[0]*t5 - t7*t8 + z[3]*x[2]*F[0]*y[1]
		    + t13*t8 + F[1]*z[2]*t16 + t18*t19 - x[1]*z[0]*t22 + t24*F[2]*y[1]
		    - t27*F[0]*y[3] - z[1]*F[2]*t16 - t32*z[2]*y[1],
		t61 = -F[0]*y[2]*x[1]*z[3] + z[0]*y[2]*x[1]*F[3] - F[1]*y[2]*t24
		    - t18*t44 - t46*t5 + t46*y[2]*z[3] + t27*y[0]*F[3] - t52*y[2]*F[3]
		    + t52*t22 + z[1]*y[2]*t32 + t58*t7 - t58*t13,
		t83 = z[1]*y[3]*x[0] + t65*y[1] - t65*y[3] - y[1]*z[3]*x[0]
		    - t70*z[1] + t70*z[3] - t73*z[3] + t75*z[1] + t77*z[1]
		    - t77*z[3] + y[1]*x[3]*z[0] - t75*z[0],
		t101 = -z[1]*x[3]*y[0] + t86*y[3] + z[3]*x[2]*y[1] + t73*z[0]
		     - t86*y[0] + t44*x[1] - t19*x[1] + t94*y[3] - t27*y[3]
		     + t97*y[0] - t97*y[1] - t94*y[1];
	
	return (t35+t61)/(t83+t101);
}


/* Subdivide space into simplices */
void *subdivide(double **m, int n)
{
	double x[3];
	int i, j, ns = 1;
	int v0, v1, v2, v3;
	partition *p = (partition *)xmalloc(sizeof(partition));
	simplex *s = (simplex *)xmalloc((3*n+1)*sizeof(simplex));
	
	static double S0[4][3] = {
		{0.0, 0.0, 0.0},
		{XMAX, 0.0, 0.0},
		{0.0, YMAX, 0.0},
		{0.0, 0.0, ZMAX}
	};
	
	/* Define simplex surrounding area of interest */
	for (i = 0; i < 4; i++)
	    for (j = 0; j < 3; j++) {
		m[j][n+i] = S0[i][j];
		m[j+3][n+i] = 0.0;
	}
	
	mksimplx(s, m, n, n+1, n+2, n+3);
	
	/* Subdivide space into simplices */
	for (i = 0; i < n; i++) {
		x[0] = m[0][i]; x[1] = m[1][i]; x[2] = m[2][i];
		
		for (j = 0; j < ns; j++) if (within(x, &(s[j]))) {
			v0 = s[j].vertex[0]; v1 = s[j].vertex[1];
			v2 = s[j].vertex[2]; v3 = s[j].vertex[3];
			
			mksimplx(&(s[j]), m, i, v1, v2, v3);
			mksimplx(&(s[ns++]), m, v0, i, v2, v3);
			mksimplx(&(s[ns++]), m, v0, v1, i, v3);
			mksimplx(&(s[ns++]), m, v0, v1, v2, i);
			
			break;
		}
		
		if (j >= ns) error("Point (%g, %g, %g) outside a simplex in subdivide()", m[0][i], m[1][i], m[2][i]);
	}
	
	/* Store info in partition table */
	p->n = n;
	p->m = m;
	p->ns = ns;
	p->simplx = s;
	
	return (void *)(p);
}


/* Interpolate function in 3d using linear simplex interpolation */
void interp3d(void *part, double x[], double z[])
{
	int i, j;
	partition *p = (partition *)(part);
	
	for (j = 0; j < p->ns; j++) if (within(x, &(p->simplx[j]))) {
		for (i = 0; i < 3; i++)
			z[i] = interp_simplex(x, &(p->simplx[j]), p->m, i);
		return;
	}
	
	if (j >= p->ns) /* We are outside of the interpolation range */
		z[0] = z[1] = z[2] = 0.0;
}
