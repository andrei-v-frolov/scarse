/* $Id: util.h,v 1.4 2001/02/05 03:59:41 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Numerical and utility functions declarations.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include <stdio.h>
#include <stddef.h>


/**********************************************************************/

#ifndef __UTIL_H__
#define __UTIL_H__

#define TINY 1.0e-8

#define ppow(x,g) ((x > 0.0) ? pow(x, g) : -pow(-x, g))


/* Utility routines (util.c) */

extern char *program_name;
extern char *usage_msg[];

void usage();
void warning(char *fmt, ...);
void error(char *fmt, ...);
void fatal(char *msg);

void *xmalloc(size_t size);
void *xrealloc(void *addr, size_t size);
char *xstrdup(const char *s);

FILE *zfopen(const char *file, const char *mode);
FILE *xfopen(const char *file, const char *mode);
FILE *xfetch(const char *prefix, const char *file, const char *mode);

void expandopt(int *argcp, char ***argvp, char *options);
void readopt(int *argcp, char ***argvp, char *file, char *name);

double *vector(unsigned long n);
double *grow_vector(double *v, unsigned long n);
int *ivector(unsigned long n);
int *grow_ivector(int *v, unsigned long n);
unsigned long *uvector(unsigned long n);
unsigned long *grow_uvector(unsigned long *v, unsigned long n);
void free_vector(void *v);

double **matrix(unsigned long nr, unsigned long nc);
double **grow_matrix(double **m, unsigned long nr, unsigned long nc);
void free_matrix(double **m);


/* Quick sort and selection routines (sort.c) */

void sort(unsigned long n, double arr[]);
void sort2(unsigned long n, double arr[], double brr[]);
void indexx(unsigned long n, double arr[], unsigned long indx[]);
void rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
double seln(unsigned long k, unsigned long n, double arr[]);

double avg(unsigned long n, double arr[]);
double mean(unsigned long n, double arr[]);
double median(unsigned long n, double arr[]);
double stddev(unsigned long n, double arr[], double x);
double absdev(unsigned long n, double arr[], double x);


/* 3x3 matrix operations (matrix.c) */

double det33(double M[3][3]);
void inv33(double M[3][3], double N[3][3]);
void mult33(double M[3][3], double N[3][3], double L[3][3]);
void apply33(double M[3][3], double A[3], double B[3]);
double qform33(double M[3][3], double A[3], double B[3]);
void diag33(double A[3], double M[3][3]);


/* Data fitting and approximation routines (fit.c) */

void gaussj(double **A, int n, double **B, int m);

double **new_amoeba(double x[], int n, double (*func)(double []), double lambda);
void restart_amoeba(double **S, int n, double (*func)(double []), double lambda);
void anneal(double **S, int n, double (*func)(double []), double T0, int maxsteps, double tol);

double **best_fit(double **x, int n, void (*basis)(double [], double []), int d);
void approx(double **A, void (*basis)(double [], double []), int d, double x[], double y[]);

double *fit_curve(double **data, int n);
int within_range(double p[], double y);
double lu_curve(double p[], double y);
double lu_curve_1(double p[], double x);
double curve_dydx(double p[], double x);

void best_linear_fit(double **x, int n, double M[3][3]);

double **best_poly_fit(double **x, int n);
void poly_approx(double **A, double x[], double y[]);


/* Multi-dimensional interpolation routines (interp.c) */

void derivs(double **m, unsigned long n);
double interp1d(double **m, unsigned long n, double x);

void *subdivide(double **m, int n);
void interp3d(void *part, double x[], double z[]);


#endif /* __UTIL_H__ */
