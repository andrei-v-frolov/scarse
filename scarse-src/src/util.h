/* $Id: util.h,v 1.10 2005/09/24 01:20:45 afrolov Exp $ */

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


/* Shared definitions */

#ifndef TINY /* machine precision */
#define TINY 2.2204460492503131e-16
#endif /* TINY */

#ifndef HUGE /* largest number */
#define HUGE 1.7976931348623158e308
#endif /* HUGE */

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

void vcopy(double src[], double dest[], unsigned long n);
void vgamma(double src[], double dest[], unsigned long n, double gamma);

double **matrix(unsigned long nr, unsigned long nc);
double **grow_matrix(double **m, unsigned long nr, unsigned long nc);
void free_matrix(double **m);

void sort(unsigned long n, double arr[]);

double avg(unsigned long n, double arr[]);
double stddev(unsigned long n, double arr[], double x);


/* 3x3 matrix operations (matrix.c) */

void vcopy3(double A[3], double B[3]);
void diag33(double A[3], double M[3][3]);

double det33(double M[3][3]);
void inv33(double M[3][3], double N[3][3]);
void mult33(double M[3][3], double N[3][3], double L[3][3]);
void apply33(double M[3][3], double A[3], double B[3]);
double qform33(double M[3][3], double A[3], double B[3]);
void biscale33(double A[3], double M[3][3], double B[3], double N[3][3]);


/* Data fitting and approximation routines (fit.c) */

double vmin(double p[], double xi[], int n, double (*f)(double []), double eps);
double mmin(double p[], double **e, int n, double (*f)(double []), double eps);

double lu_curve(double p[], double x);
double lu_curve_1(double p[], double y);
double *fit_curve(double **data, int n);

void fit_matrix(double **data, int n, double M[3][3]);

void gaussj(double **A, int n, double **B, int m);

void evalf(double **A, void (*basis)(double [], double []), int d, double x[], double y[]);
double **lsq_fit(double **x, int n, void (*basis)(double [], double []), int d);

double **fit_poly(double **x, int n);
void evalf_poly(double **A, double x[], double y[]);


#endif /* __UTIL_H__ */
