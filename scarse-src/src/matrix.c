/* $Id: matrix.c,v 1.1 2001/01/26 22:45:31 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Matrix operations on 3x3 matrices.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include "util.h"


/* determinant of matrix:  det(M) */
double det33(double M[3][3])
{
      return M[0][0]*M[1][1]*M[2][2]-M[0][0]*M[1][2]*M[2][1]
              -M[1][0]*M[0][1]*M[2][2]+M[1][0]*M[0][2]*M[2][1]
	      +M[2][0]*M[0][1]*M[1][2]-M[2][0]*M[0][2]*M[1][1];
}

/* invert matrix:  N = M^-1 */
void inv33(double M[3][3], double N[3][3])
{
	double t4 = M[0][0]*M[1][1], t6 = M[0][0]*M[1][2], t8 = M[0][1]*M[1][0];
	double t10 = M[0][2]*M[1][0], t12 = M[0][1]*M[2][0], t14 = M[0][2]*M[2][0];
	double t17 = 1/(t4*M[2][2]-t6*M[2][1]-t8*M[2][2]+t10*M[2][1]+t12*M[1][2]-t14*M[1][1]);
	
	N[0][0] = (M[1][1]*M[2][2]-M[1][2]*M[2][1])*t17;
	N[0][1] = -(M[0][1]*M[2][2]-M[0][2]*M[2][1])*t17;
	N[0][2] = -(-M[0][1]*M[1][2]+M[0][2]*M[1][1])*t17;
	N[1][0] = -(M[1][0]*M[2][2]-M[1][2]*M[2][0])*t17;
	N[1][1] = (M[0][0]*M[2][2]-t14)*t17;
	N[1][2] = -(t6-t10)*t17;
	N[2][0] = -(-M[1][0]*M[2][1]+M[1][1]*M[2][0])*t17;
	N[2][1] = -(M[0][0]*M[2][1]-t12)*t17;
	N[2][2] = (t4-t8)*t17;
}

/* multiply two matrices:  L = M*N */
void mult33(double M[3][3], double N[3][3], double L[3][3])
{
	L[0][0] = M[0][0]*N[0][0]+M[0][1]*N[1][0]+M[0][2]*N[2][0];
	L[0][1] = M[0][0]*N[0][1]+M[0][1]*N[1][1]+M[0][2]*N[2][1];
	L[0][2] = M[0][0]*N[0][2]+M[0][1]*N[1][2]+M[0][2]*N[2][2];
	L[1][0] = M[1][0]*N[0][0]+M[1][1]*N[1][0]+M[1][2]*N[2][0];
	L[1][1] = M[1][0]*N[0][1]+M[1][1]*N[1][1]+M[1][2]*N[2][1];
	L[1][2] = M[1][0]*N[0][2]+M[1][1]*N[1][2]+M[1][2]*N[2][2];
	L[2][0] = M[2][0]*N[0][0]+M[2][1]*N[1][0]+M[2][2]*N[2][0];
	L[2][1] = M[2][0]*N[0][1]+M[2][1]*N[1][1]+M[2][2]*N[2][1];
	L[2][2] = M[2][0]*N[0][2]+M[2][1]*N[1][2]+M[2][2]*N[2][2];
}

/* apply matrix to vector:  B = M*A */
void apply33(double M[3][3], double A[3], double B[3])
{
	B[0] = M[0][0]*A[0]+M[0][1]*A[1]+M[0][2]*A[2];
	B[1] = M[1][0]*A[0]+M[1][1]*A[1]+M[1][2]*A[2];
	B[2] = M[2][0]*A[0]+M[2][1]*A[1]+M[2][2]*A[2];
}

/* calculate quadratic form q = At*(M*B) */
double qform33(double M[3][3], double A[3], double B[3])
{
	return	(A[0]*M[0][0]+A[1]*M[1][0]+A[2]*M[2][0])*B[0] +
		(A[0]*M[0][1]+A[1]*M[1][1]+A[2]*M[2][1])*B[1] +
		(A[0]*M[0][2]+A[1]*M[1][2]+A[2]*M[2][2])*B[2];
}

/* make diagonal matrix:  M = diag(A) */
void diag33(double A[3], double M[3][3])
{
	M[0][0] = A[0]; M[0][1] = 0.0;  M[0][2] = 0.0;
	M[1][0] = 0.0;  M[1][1] = A[1]; M[1][2] = 0.0;
	M[2][0] = 0.0;  M[2][1] = 0.0;  M[2][2] = A[2];
}
