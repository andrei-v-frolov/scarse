/* $Id: util.c,v 1.4 2005/10/05 06:29:26 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Error handlers, vectors and matrices, other misc stuff.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include "scarse.h"



/******************* Basic error handlers *****************************/

void usage()
{
	int i = 0; char *s;
	
	while ((s = usage_msg[i++]))
		fprintf(stderr, "%s\n", s);
	
	exit(1);
}

void warning(char *fmt, ...)
{
	va_list args;
	
	fprintf(stderr,"%s: ", program_name);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}

void error(char *fmt, ...)
{
	va_list args;
	
	fprintf(stderr, "%s: ", program_name);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	
	exit(-1);
}

void fatal(char *msg)
{
	error("Fatal error: %s\nTerminating process...", msg);
}


/******************* Foo or die routines ******************************/

/* allocate memory or die */
void *xmalloc(size_t size)
{
	register void *p = malloc(size);
	
	if (!p) fatal("virtual memory exhausted in malloc()");
	return p;
}

/* reallocate memory or die */
void *xrealloc(void *addr, size_t size)
{
	register void *p = realloc(addr, size);
	
	if (!p) fatal("virtual memory exhausted in realloc()");
	return p;
}

/* duplicate string or die */
char *xstrdup(const char *s)
{
	register char *p = strdup(s);
	
	if (!p) fatal("virtual memory exhausted in strdup()");
	return p;
}


/******************* File access wrappers *****************************/

/* open file (will read compressed files) */
FILE *zfopen(const char *file, const char *mode)
{
	char *f, *p;
	FILE *fp = fopen(file, mode);
	
	/* try opening as uncompressed stream */
	if (*mode != 'r' || strchr(mode, '+')) return fp;
	if (fnmatch("*.{gz,Z}", file, 0) && fp) return fp;
	
	/* failed; try adding suffixes */
	f = xstrdup(file);
	if (!fp) { free(f); asprintf(&f, "%s.gz", file); fp = fopen(f, mode); }
	if (!fp) { free(f); asprintf(&f, "%s.Z", file); fp = fopen(f, mode); }
	
	/* reopen as compressed input stream */
	if (fp) { fclose(fp); asprintf(&p, "zcat %s", f); fp = popen(p, mode); free(p); }
	
	free(f);
	
	return fp;
}

/* open file or die */
FILE *xfopen(const char *file, const char *mode)
{
	register FILE *fp = zfopen(file, mode);
	
	if (!fp) error("Can't open file '%s'", file);
	return fp;
}

/* find & open file for input */
FILE *xfetch(const char *prefix, const char *file, const char *mode)
{
	int i;
	FILE *fp = strcmp(file, "-") ? zfopen(file, mode) : stdin;
	char *t, *path[] = { ".", getenv("CMS_DATADIR"), "~/.scarse", DATADIR, "/usr/share/cms" };
	
	if (*mode != 'r' || strchr(mode, '+'))
		error("Internal error: use xfetch() for input only");
	
	/* if not absolute filename, search the path */
	asprintf(&path[1], "%s/.scarse", getenv("HOME"));
	if (*file != '/') for (i = 0; !fp && i < sizeof(path)/sizeof(char *); i++) if (path[i]) {
		asprintf(&t, "%s/%s%s%s", path[i], prefix ? prefix : "", prefix ? "/" : "", file);
		fp = zfopen(t, mode);
		free(t);
	}
	
	if (!fp) error("Can't open file '%s'", file);
	return fp;
}


/******************* Option expansion *********************************/

/* Expand option string into ARGV, as if it was inserted on the command line */
void expandopt(int *argcp, char ***argvp, char *options)
{
	#define OPT_MAX 64
	
	FILE *fp;
	int i, n = 0, null = 0;
	char *buffer, *opt[OPT_MAX], **nargv;
	
	/* invoke shell to tokenize option string */
	asprintf(&buffer, "printf '%%s\\0' %s", options); fp = popen(buffer, "r");
	if (!fp) error("Could not invoke shell to parse options string '%s'", options);
	
	do { opt[n] = NULL; null = 0; } while (getdelim(&(opt[n]), &null, 0, fp) != -1 && ++n < OPT_MAX);
	if (n == OPT_MAX) warning("Option number limit (%i) reached in expandopt() - ignoring the rest", OPT_MAX);
	
	/* merge it into argv at current optind */
	nargv = xmalloc((*argcp + n + 1)*sizeof(char *));
	
	for (i = 0; i < optind; i++) nargv[i] = (*argvp)[i];
	for (i = 0; i < n; i++) nargv[optind+i] = opt[i];
	for (i = optind; i <= *argcp; i++) nargv[i+n] = (*argvp)[i];
	
	/* update argv and clean up */
	*argvp = nargv; *argcp += n;
	
	free(buffer);
	pclose(fp);
	
	#undef OPT_MAX
}

/* Read long option definition from file, and expand it */
void readopt(int *argcp, char ***argvp, char *file, char *name)
{
	char *p, *q;
	int found = 0;
	size_t bsize = 128;
	char *buffer = (char *)xmalloc(bsize);
	FILE *fp = xfetch("etc", file, "r");
	
	while (getline(&buffer, &bsize, fp) != -1) {
		/* skip over comments and empty lines */
		if (*buffer == '#' || *buffer == '\n' || *buffer == '\r') continue;
		
		/* options file format: name options... */
		if (sscanf(buffer, " %as %a[^\n\r]", &p, &q) != 2)
			error("%s: Error parsing options file", file);
		
		/* expand matching macros */
		if (!fnmatch(p, name, 0)) { found = 1; expandopt(argcp, argvp, q); }
		
		free(p);
		free(q);
	}
	
	if (!found) { warning("invalid long option -- %s", name); usage(); }
	
	free(buffer);
	fclose(fp);
}


/******************* Vectors and matrices *****************************/

/* allocate a vector with subscript range v[0..n] */
double *vector(unsigned long n)
{
	register double *v = (double *)xmalloc((size_t)(n*sizeof(double)));
	
	return v;
}

/* grow a vector to subscript range v[0..n] */
double *grow_vector(double *v, unsigned long n)
{
	v = (double *)xrealloc(v, (size_t)(n*sizeof(double)));
	
	return v;
}

/* allocate an int vector with subscript range v[0..n] */
int *ivector(unsigned long n)
{
	register int *v = (int *)xmalloc((size_t)(n*sizeof(int)));
	
	return v;
}

/* grow an int vector to subscript range v[0..n] */
int *grow_ivector(int *v, unsigned long n)
{
	v = (int *)xrealloc(v, (size_t)(n*sizeof(int)));
	
	return v;
}

/* allocate an unsigned long vector with subscript range v[0..n] */
unsigned long *uvector(unsigned long n)
{
	register unsigned long *v = (unsigned long *)xmalloc((size_t)(n*sizeof(unsigned long)));
	
	return v;
}

/* grow an unsigned long vector to subscript range v[0..n] */
unsigned long *grow_uvector(unsigned long *v, unsigned long n)
{
	v = (unsigned long *)xrealloc(v, (size_t)(n*sizeof(unsigned long)));
	
	return v;
}

/* free a vector allocated with vector() */
void free_vector(void *v)
{
	free((void *)(v));
}


/* copy vector data with subscript range [0..n] */
void vcopy(double src[], double dest[], unsigned long n)
{
	register unsigned long i;
	
	for (i = 0; i < n; i++) dest[i] = src[i];
}

/* apply gamma transformation to a vector */
void vgamma(double src[], double dest[], unsigned long n, double gamma)
{
	register unsigned long i;
	
	for (i = 0; i < n; i++) dest[i] = ppow(src[i], gamma);
}


/* allocate a matrix with subscript range m[0..nr][0..nc] */
double **matrix(unsigned long nr, unsigned long nc)
{
	register unsigned long i;
	register double **m = (double **)xmalloc((size_t)(nr*sizeof(double *)+sizeof(unsigned long)));
	
	*(((unsigned long *)(m))++) = nr;
	
	for (i = 0; i < nr; i++) m[i] = vector(nc);
	
	return m;
}

/* grow a matrix to subscript range m[0..nr][0..nc] */
double **grow_matrix(double **m, unsigned long nr, unsigned long nc)
{
	register unsigned long i;
	unsigned long old_nr = *(--((unsigned long *)(m)));
	
	/* Reallocate row index if necessary */
	if (nr != old_nr)
		m = (double **)xrealloc(m, (size_t)(nr*sizeof(double *)+sizeof(unsigned long)));
	
	*(((unsigned long *)(m))++) = nr;
	
	/* Reallocate rows */
	for (i = 0; i < old_nr; i++) m[i] = grow_vector(m[i], nc);
	for (i = old_nr; i < nr; i++) m[i] = vector(nc);
	
	return m;
}

/* free a matrix allocated by matrix() */
void free_matrix(double **m)
{
	register unsigned long i;
	unsigned long nr = *((unsigned long *)(m)-1);
	
	for (i = 0; i < nr; i++) free_vector(m[i]);
	
	free((void *)((unsigned long *)(m)-1));
}



/******************* Sorting and averages *****************************/

/* sort array into ascending numerical order using heap sort */
void sort(unsigned long n, double arr[])
{
	unsigned long i, j, k = n >> 1, l = n-1; double t;
	
	while (n > 1) {
		if (k > 0) { t = arr[--k]; } else {
			t = arr[l]; arr[l] = arr[0];
			if (--l == 0) { arr[0] = t; break; }
		}
		
		i = k; j = k+k+1;
		
		while (j <= l) {
			if (j < l && arr[j] < arr[j+1]) j++;
			
			if (t < arr[j]) { arr[i] = arr[j]; i = j; j = j+j+1; }
			else break;
		}
		
		arr[i] = t;
	}
}

/* return average of median-filtered values in array */
double avg(unsigned long n, double arr[])
{
	double S = 0.0;
	unsigned long i, t = n >> 4;
	
	sort(n, arr);
	
	/* Ignore tails of distribution when averaging */
	for (i = t; i < n-t; i++) S += arr[i];
	
	return S/(double)(n-2*t);
}

/* return standard deviation */
double stddev(unsigned long n, double arr[], double x)
{
	double S = 0.0;
	unsigned long i;
	
	for (i = 0; i < n; i++) S += (arr[i] - x) * (arr[i] - x);
	
	return sqrt(S/(double)(n-1));
}
