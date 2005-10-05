/* $Id: scarse.h,v 1.7 2005/10/05 06:29:26 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Function & external variables declarations.
 * 
 * Copyright (C) 1999-2005 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <afrolov@stanford.edu>
 * 
 */


/**********************************************************************/

#ifndef __SCARSE_H__
#define __SCARSE_H__

/* configuration */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* standard headers */
#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include <unistd.h>
#include <stdarg.h>
#include <fnmatch.h>

/* required libraries */
#include <icc.h>
#include <tiffio.h>

/* our own headers */
#include "spaces.h"
#include "imageio/imageio.h"


/**********************************************************************/

/* Shared definitions */

#ifndef TINY /* machine precision */
#define TINY 2.2204460492503131e-16
#endif /* TINY */

#ifndef HUGE /* largest number */
#define HUGE 1.7976931348623158e308
#endif /* HUGE */

#define ppow(x,g) ((x >= 0.0) ? pow(x, g) : -pow(-x, g))


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


/* Data fitting and approximation routines (fit.c) */

double vmin(double p[], double xi[], int n, double (*f)(double []), double eps);
double mmin(double p[], double **e, int n, double (*f)(double []), double eps);

double lu_curve(double p[], double x);
double lu_curve_1(double p[], double y);
double *fit_curve(double **data, int n);

int *fit_matrix(double **data, int n, double M[3][3], double **C);

void gaussj(int n, double **A, double *B, double *X);

void evalf(double **A, void (*basis)(double [], double []), int d, double x[], double y[]);
double **lsq_fit(double **x, int n, void (*basis)(double [], double []), int d);

double **fit_poly(double **x, int n);
void evalf_poly(double **A, double x[], double y[]);

void prepint3d(double **m, int n, double Q[3][3]);
void interp3d(double **m, int n, double Q[3][3], double x[], double z[]);


/**********************************************************************/

/* Calibration targets (targets.c) */

typedef struct {
	int pts;
	struct _target_data { char *label; double RGB[6], XYZ[3]; } *data;
	int rows, cols, **layout;
	int subrows, *subridx;
	int subcols, *subcidx;
	int graypts, *grayscale;
} target;

void list_targets(char *type);
void find_target(char *type, char *batch, char **data, char **layout);

void parse_IT87_target(target *tg, char *data_file, char *layout_file);
void render_IT87_target(target *tg, char *file, char *geometry);
void read_IT87_target(target *tg, char *file, char *geometry);


#endif /* __SCARSE_H__ */
