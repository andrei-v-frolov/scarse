/* $Id: sort.c,v 1.2 2001/01/28 02:23:20 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Sorting, indexing, ranking, selection, and averages.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include <math.h>
#include "util.h"


/**********************************************************************/

#define WORTHY 16

#define PUSH(i,j) {if (sp>=ss) stack=grow_uvector(stack, ss<<=1); stack[sp++]=i; stack[sp++]=j;}
#define POP(i,j) {j=stack[--sp]; i=stack[--sp];}


/* Sort array into ascending numerical order */
void sort(unsigned long n, double arr[])
{
	register double a, t;
	register unsigned long i, j;
	
	unsigned long k = 0, l = n-1, m;
	int sp = 0, ss = 8*sizeof(unsigned long);
	unsigned long *stack = uvector(ss);
	
	#define XCHG(i,j) {t=arr[i]; arr[i]=arr[j]; arr[j]=t;}
	
	
	for (;;) if (l-k > WORTHY) {
		/* interval [k .. l] is worthy of preordering */
		
		m = (k+l) >> 1;
		
		if (arr[k] > arr[m]) XCHG(k, m)
		if (arr[k] > arr[l]) XCHG(k, l)
		if (arr[l] > arr[m]) XCHG(l, m)
		XCHG(m, l-1)
		
		/* Inner ordering loop */
		i = k;
		j = l-1;
		a = arr[l];
		
		for (;;) {
			while (arr[++i] < a);
			while (arr[--j] > a);
			if (j < i) break;
			XCHG(i, j)
		}
		
		XCHG(i, l)
		
		PUSH(i, l)	/* order [i .. l] later */
		l = j;		/* order [k .. j] next */
	} else {
		if (!sp) break; POP(k, l)
	}
	
	/* Do final insertion sort */
	for (j = 1; j < n; j++) {
		i = j;
		a = arr[j];
		
		while (i-- > 0 && arr[i] > a) arr[i+1] = arr[i];
		
		arr[i+1] = a;
	}
	
	free_vector(stack);
	
	#undef XCHG
}


/* Sort array into ascending numerical order, rearranging the second one the same way */
void sort2(unsigned long n, double arr[], double brr[])
{
	register double a, t;
	register unsigned long i, j;
	
	unsigned long k = 0, l = n-1, m;
	int sp = 0, ss = 8*sizeof(unsigned long);
	unsigned long *stack = uvector(ss);
	
	#define XCHG(i,j) {t=arr[i]; arr[i]=arr[j]; arr[j]=t;\
			   t=brr[i]; brr[i]=brr[j]; brr[j]=t;}
	
	
	for (;;) if (l-k > WORTHY) {
		/* interval [k .. l] is worthy of preordering */
		
		m = (k+l) >> 1;
		
		if (arr[k] > arr[m]) XCHG(k, m)
		if (arr[k] > arr[l]) XCHG(k, l)
		if (arr[l] > arr[m]) XCHG(l, m)
		XCHG(m, l-1)
		
		/* Inner ordering loop */
		i = k;
		j = l-1;
		a = arr[l];
		
		for (;;) {
			while (arr[++i] < a);
			while (arr[--j] > a);
			if (j < i) break;
			XCHG(i, j)
		}
		
		XCHG(i, l)
		
		PUSH(i, l)	/* order [i .. l] later */
		l = j;		/* order [k .. j] next */
	} else {
		if (!sp) break; POP(k, l)
	}
	
	/* Do final insertion sort */
	for (j = 1; j < n; j++) {
		i = j;
		a = arr[j];
		t = brr[j];
		
		while (i-- > 0 && arr[i] > a) {
			arr[i+1] = arr[i];
			brr[i+1] = brr[i];
		}
		
		arr[i+1] = a;
		brr[i+1] = t;
	}
	
	free_vector(stack);
	
	#undef XCHG
}


/* Index an array into ascending numerical order */
void indexx(unsigned long n, double arr[], unsigned long indx[])
{
	register double a;
	register unsigned long i, j, t;
	
	unsigned long k = 0, l = n-1, m;
	int sp = 0, ss = 8*sizeof(unsigned long);
	unsigned long *stack = uvector(ss);
	
	#define XCHG(i,j) {t=indx[i]; indx[i]=indx[j]; indx[j]=t;}
	
	
	/* Initialize index table */
	for (i = 0; i < n; i++) indx[i] = i;
	
	for (;;) if (l-k > WORTHY) {
		/* interval [k .. l] is worthy of preordering */
		
		m = (k+l) >> 1;
		
		if (arr[indx[k]] > arr[indx[m]]) XCHG(k, m)
		if (arr[indx[k]] > arr[indx[l]]) XCHG(k, l)
		if (arr[indx[l]] > arr[indx[m]]) XCHG(l, m)
		XCHG(m, l-1)
		
		/* Inner ordering loop */
		i = k;
		j = l-1;
		a = arr[indx[l]];
		
		for (;;) {
			while (arr[indx[++i]] < a);
			while (arr[indx[--j]] > a);
			if (j < i) break;
			XCHG(i, j)
		}
		
		XCHG(i, l)
		
		PUSH(i, l)	/* order [i .. l] later */
		l = j;		/* order [k .. j] next */
	} else {
		if (!sp) break; POP(k, l)
	}
	
	/* Do final insertion sort */
	for (j = 1; j < n; j++) {
		i = j;
		t = indx[j];
		a = arr[t];
		
		while (i-- > 0 && arr[indx[i]] > a) indx[i+1] = indx[i];
		
		indx[i+1] = t;
	}
	
	free_vector(stack);
	
	#undef XCHG
}


/* Construct table of ranks from index table */
void rank(unsigned long n, unsigned long indx[], unsigned long irank[])
{
	register unsigned long j;
	
	for (j = 0; j < n; j++)
		irank[indx[j]] = j;
}


/* Return r-th smallest value in array (destructive) */
/* Gives the same answer as sort(arr); seln = arr[r]; */
double seln(unsigned long r, unsigned long n, double arr[])
{
	register double a, t;
	register unsigned long i, j;
	unsigned long k = 0, l = n-1, m;
	
	#define XCHG(i,j) {t=arr[i]; arr[i]=arr[j]; arr[j]=t;}
	
	
	for (;;) if (l-k > 1) {
		m = (k+l) >> 1;
		
		if (arr[k] > arr[m]) XCHG(k, m)
		if (arr[k] > arr[l]) XCHG(k, l)
		if (arr[l] > arr[m]) XCHG(l, m)
		XCHG(m, l-1)
		
		/* Inner ordering loop */
		i = k;
		j = l-1;
		a = arr[l];
		
		for (;;) {
			while (arr[++i] < a);
			while (arr[--j] > a);
			if (j < i) break;
			XCHG(i, j)
		}
		
		XCHG(i, l)
		
		if (i >= r) l = i;
		if (i <= r) k = i;
	} else {
		if (l-k == 1 && arr[l] < arr[k]) XCHG(l, k)
		return arr[r];
	}
	
	#undef XCHG
}



/**********************************************************************/

/* Return average of median-filtered values in array */
double avg(unsigned long n, double arr[])
{
	double S = 0.0;
	unsigned long i, t = n >> 4;
	
	sort(n, arr);
	
	/* Ignore tails of distribution when averaging */
	for (i = t; i < n-t; i++) S += arr[i];
	
	return S/(double)(n-2*t);
}

/* Return mean of values in array */
double mean(unsigned long n, double arr[])
{
	double S = 0.0;
	unsigned long i;
	
	for (i = 0; i < n; i++) S += arr[i];
	
	return S/(double)(n);
}

/* Return median of values in array */
double median(unsigned long n, double arr[])
{
	return seln(n/2, n, arr);
}

/* Return standard deviation */
double stddev(unsigned long n, double arr[], double x)
{
	double S = 0.0;
	unsigned long i;
	
	for (i = 0; i < n; i++) S += (arr[i] - x) * (arr[i] - x);
	
	return sqrt(S/(double)(n-1));
}

/* Return absolute deviation */
double absdev(unsigned long n, double arr[], double x)
{
	double S = 0.0;
	unsigned long i;
	
	for (i = 0; i < n; i++) S += fabs(arr[i] - x);
	
	return S/(double)(n);
}
