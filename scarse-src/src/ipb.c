/* $Id: ipb.c,v 1.1.1.1 2001/01/26 22:45:31 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * ICC profile builder - build profile based on measured data.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

/* CREDITS:
 *   Derived from icclib & examples by Graeme W. Gill
 */

/* TODO:
 *
 */

#define SELF "ipb"

#define _GNU_SOURCE
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <icc.h>
#include "scarse.h"



/* Usage */
char *program_name = SELF;
char *usage_msg[] = {
	"ICC profile builder, Version " VERSION,
	"Author: Andrei Frolov <andrei@phys.ualberta.ca>",
	"",
	"Usage: " SELF " [options] profile.icm",
	"",
	" General options:",
	"  -h		print this message",
	"  -v		verbosity (cumulative)",
	"  --longopt	expand option macro 'longopt' as defined in etc/" SELF ".options",
	"",
	" Color space options:",
	"  -c class	profile class (i = input, d = display, o = output)",
	"			      (c = color space conversion, a = abstract)",
	"  -[i|o] space[:gamma]",
	"		set input (device) and output (connection) color spaces",
	"		space = {XYZ, Lab, Luv, Yxy, RGB, GRAY, HSV, CMY, CMYK}",
	"		if followed by gamma, power law scaling will be applied",
	"  -p [white|red|green|blue]:x,y | standard",
	"		set (x,y) values of white point and RGB primaries",
	"		known standards are:",
	"		    illuminants: D50 (default), D55, D65, D75, D93, A, B, C, E",
	"		  RGB primaries: 700/525/450nm, EBU, HDTV, P22, Trinitron",
	"		     RGB spaces: Adobe (default), Apple, ColorMatch, CIE,",
	"		                 NTSC, PAL/SECAM, sRGB, SMPTE-C, WideGamut",
	"",
	" Color mapping algorithm options:",
	"  -C file	calibration data; '-' means read from stdin",
	"		required format of the data is described elsewhere",
	"  -U		generate unidirectional (device -> PCS) profile only",
	"  -E gamma	LUT shadow expansion exponent (none=1.0, default=3.0)",
	"  -M		generate matrix based profile only, do not output LUT",
	"  -L		ignore color correction LUT from calibration data",
	"",
	" ICC profile options:",
	"  -m manufacturer[:model] of device for which profile is created",
	"  -d string	description of device profile",
	"  -r string	copyright string of device profile",
	"",
	"  -b bits	LUT bits per sample (bits = 8 or 16)",
	"  -s 1D[:3D]	LUT table size, for 1D and 3D tables",
	NULL
};



/* Options */
static int verbose = 0;

static icProfileClassSignature class = icSigInputClass;
static icColorSpaceSignature     ins = icSigRgbData;
static icColorSpaceSignature    outs = icSigXYZData;

static double        primaries[4][2];

static double      media_white_pt[3] = {0.9642, 1.0000, 0.8249};
static double      media_black_pt[3] = {0.0, 0.0, 0.0};

static FILE        *calibration_data = NULL;
static int                invertible = 1;
static int                ignore_lut = 0;
static int               matrix_only = 0;
static double   lut_shadow_expansion = 3.0;

static char            *manufacturer = "none";
static char                   *model = "none";
static char             *description = "";
static char               *copyright = "";

static int             bitspersample = 16;
static int               lut_size_1d = 256;
static int               lut_size_3d = 33;

static char            *profile = NULL;



/******************* Transform functions ******************************/

/* Transform data */
static int              ins_channels = 0;
static double              ins_gamma = 1.0;
static transform             ins2XYZ = NULL,
                             XYZ2ins = NULL;

static int             outs_channels = 0;
static double             outs_gamma = 1.0;
static transform            outs2XYZ = NULL,
                            XYZ2outs = NULL;

static int	use_ins_curves = 0,
		use_outs_curves = 0;

static int	ins_curves_n[MAXCHANNELS],
		outs_curves_n[MAXCHANNELS];

static double	**ins_curves[MAXCHANNELS],
		**ins_curves_1[MAXCHANNELS],
		**outs_curves[MAXCHANNELS],
		**outs_curves_1[MAXCHANNELS];

static int          calibration_data_pts = 0;
static double	*calibration_data_vector = NULL;

static double	L[3][3]     = { {1.0, 0.0, 0.0},
				{0.0, 1.0, 0.0},
				{0.0, 0.0, 1.0} },
		L1[3][3]    = { {1.0, 0.0, 0.0},
				{0.0, 1.0, 0.0},
				{0.0, 0.0, 1.0} };

static int	use_higher_order_approximation = 0;
static double	**P = NULL, **P1 = NULL;

static void	*ELUT = NULL, *ELUT_1 = NULL;


/* Translate data through curve */
void xlate_curve(double in[], double out[], int channels, double **m[], int n[])
{
	int i;
	
	for (i = 0; i < channels; i++)
		out[i] = interp1d(m[i], n[i], in[i]);
}


/* Input curves */
void incurves(void *cntx, double out[], double in[])
{
	if (use_ins_curves)
		xlate_curve(in, out, ins_channels, ins_curves, ins_curves_n);
	else if (ins_gamma != 1.0)
		gamma_scale(in, out, ins_channels, ins_gamma);
	else
		identity(in, out, ins_channels);
}

/* Input curves for shadow-expanded LUT */
void incurves_expanded(void *cntx, double out[], double in[])
{
	incurves(cntx, out, in);
	
	if (lut_shadow_expansion != 1.0)
		gamma_scale(out, out, ins_channels, 1.0/lut_shadow_expansion);
}

/* Inverse of input curves */
void incurves_1(void *cntx, double out[], double in[])
{
	if (use_ins_curves)
		xlate_curve(in, out, ins_channels, ins_curves_1, ins_curves_n);
	else if (ins_gamma != 1.0)
		gamma_scale(in, out, ins_channels, 1.0/ins_gamma);
	else
		identity(in, out, ins_channels);
}

/* Output curves */
void outcurves(void *cntx, double out[], double in[])
{
	if (use_outs_curves)
		xlate_curve(in, out, outs_channels, outs_curves, outs_curves_n);
	else if (outs_gamma != 1.0)
		gamma_scale(in, out, outs_channels, 1.0/outs_gamma);
	else
		identity(in, out, outs_channels);
}

/* Inverse of output curves */
void outcurves_1(void *cntx, double out[], double in[])
{
	if (use_outs_curves)
		xlate_curve(in, out, outs_channels, outs_curves_1, outs_curves_n);
	else if (outs_gamma != 1.0)
		gamma_scale(in, out, outs_channels, outs_gamma);
	else
		identity(in, out, outs_channels);
}

/* Inverse of output curves for shadow-expanded LUT */
void outcurves_1_expanded(void *cntx, double out[], double in[])
{
	outcurves_1(cntx, out, in);
	
	if (lut_shadow_expansion != 1.0)
		gamma_scale(out, out, ins_channels, 1.0/lut_shadow_expansion);
}


/* 3D -> 3D color transformation */
void ctransform(void *cntx, double out[], double in[])
{
	int i;
	double t[3], XYZ[3], diff[3];
	
	if (lut_shadow_expansion != 1.0)
		gamma_scale(in, in, ins_channels, lut_shadow_expansion);
		/* this clobbers input, but we don't care... */
	
	(*ins2XYZ)(in, t);
	
	if (use_higher_order_approximation)
		poly_approx(P, t, XYZ); else apply33(L, t, XYZ);
	
	if (ELUT) {
		interp3d(ELUT, XYZ, diff);
		for (i = 0; i < 3; i++) XYZ[i] += diff[i];
	}
	
	(*XYZ2outs)(XYZ, out);
}

/* Inverse 3D -> 3D color transformation */
void ctransform_1(void *cntx, double out[], double in[])
{
	int i;
	double t[3], XYZ[3], diff[3];
	
	if (lut_shadow_expansion != 1.0)
		gamma_scale(in, in, ins_channels, lut_shadow_expansion);
		/* this clobbers input, but we don't care... */
	
	(*outs2XYZ)(in, XYZ);
	
	if (use_higher_order_approximation)
		poly_approx(P1, XYZ, t); else apply33(L1, XYZ, t);
	
	if (ELUT_1) {
		interp3d(ELUT_1, XYZ, diff);
		for (i = 0; i < 3; i++) XYZ[i] += diff[i];
	}
	
	(*XYZ2ins)(t, out);
}



/******************* Gamut constraint checking ************************/

/* Device gamut constrains */
static int	ins_gamut_constrained = 0,
		outs_gamut_constrained = 0;

static double	ins_gamut[2], outs_gamut[2];


/* Initialize gamut constraint to slightly bigger than default range */
void init_gamut_constraint(icColorSpaceSignature space, double *gmt)
{
	if (space != icSigHsvData)
		range(space, 0, gmt);
	else
		range(space, 2, gmt);
	
	gmt[0] = 1.1*gmt[0] - 0.1;
	gmt[1] = 1.1*gmt[1] + 0.1;
}

/* Return gamut constraints on a channel of color space */
double *gamut_constraint(icColorSpaceSignature space, int channel, double *gmt)
{
	switch (space) {
		case icSigGrayData:
		case icSigXYZData:
		case icSigRgbData:
		case icSigCmyData:
		case icSigCmykData:
			return gmt;
		case icSigYxyData:
		case icSigLabData:
		case icSigLuvData:
			return channel == 0 ? gmt : NULL;
		case icSigHsvData:
			return channel == 2 ? gmt : NULL;
		default:
			error("[%s]: Color space not supported",
				ColorSpaceSignature2str(space));
	}
	
	return NULL;
}

/* Check if pixel color satisfies gamut constraint */
int in_gamut(icColorSpaceSignature space, double gmt[], double p[])
{
	switch (space) {
		case icSigGrayData:
		case icSigYxyData:
		case icSigLabData:
		case icSigLuvData:
			return (p[0] > gmt[0]) && (p[0] < gmt[1]);
		case icSigXYZData:
		case icSigRgbData:
		case icSigCmyData:
			return  (p[0] > gmt[0]) && (p[0] < gmt[1]) &&
				(p[1] > gmt[0]) && (p[1] < gmt[1]) &&
				(p[2] > gmt[0]) && (p[2] < gmt[1]);
		case icSigCmykData:
			return  (p[0] > gmt[0]) && (p[0] < gmt[1]) &&
				(p[1] > gmt[0]) && (p[1] < gmt[1]) &&
				(p[2] > gmt[0]) && (p[2] < gmt[1]) &&
				(p[3] > gmt[0]) && (p[3] < gmt[1]);
		case icSigHsvData:
			return (p[2] > gmt[0]) && (p[2] < gmt[1]);
		default:
			error("[%s]: Color space not supported",
				ColorSpaceSignature2str(space));
	}
	
	return -1;
}



/******************* Calibration data handling ************************/

/* Handle possible degeneracy in curve data */
int degenerate(int *nptr, double **m, int io, double gmt[])
{
	int j, k, kl = -1, ku = -1, n = *nptr;
	int *eq = ivector(n), *dg = ivector(n);
	double *x = vector(n), *y = vector(n);
	double a = gmt ? gmt[0] : -0.1, b = gmt ? gmt[1] : 1.1, c;
	
	#define EEQUIV(a,b) (fabs(a-b) < 2.0e-4/n)
	#define EQUIV(i,j) (EEQUIV(m[0][i],m[0][j]) || EEQUIV(m[1][i],m[1][j]))
	
	/* Build equivalence class table */
	for (j = 0; j < n; j++) {
		eq[j] = j;
		
		for (k = 0; k < j; k++) {
			eq[k] = eq[eq[k]];
			if (EQUIV(j,k)) eq[eq[eq[k]]] = j;
		}
	}
	
	for (j = 0; j < n; j++) eq[j] = eq[eq[j]];
	
	/* Combine points within the same equivalence class */
	for (j = 0; j < n; j++) {
		dg[j] = 0;
		x[j] = 0.0;
		y[j] = 0.0;
	}
	
	for (j = 0; j < n; j++) {
		dg[eq[j]]++;
		x[eq[j]] += m[0][j];
		y[eq[j]] += m[1][j];
	}
	
	/* Remove degeneracy and write out corrected data */
	k = 0;
	for (j = 0; j < n; j++) if (dg[j]) {
		if (dg[j] == 1) {
			m[0][k] = x[j];
			m[1][k] = y[j];
		} else {
			m[0][k] = x[j]/dg[j];
			m[1][k] = y[j]/dg[j];
			
			if (gmt) {
				c = m[io][k];
				
				if ((c > a) && (c-a < b-c)) { a = c; kl = k; }
				if ((c < b) && (b-c < c-a)) { b = c; ku = k; }
			}
		}
		
		k++;
	}
	*nptr = k;
	
	/* Return gamut bounds of non-singular data if requested */
	if (gmt && k != n) {
		unsigned long *idx = uvector(k), *rnk = uvector(k);
		
		indexx(k, m[io], idx); rank(k, idx, rnk);
		
		if (kl != -1) { j = rnk[kl]; gmt[0] = (m[io][kl] + m[io][idx[j<k-1 ? j+1 : j]])/2.0; }
		if (ku != -1) { j = rnk[ku]; gmt[1] = (m[io][ku] + m[io][idx[j>0 ? j-1 : j]])/2.0; }
		
		if (verbose > 3)
			fprintf(stderr, "\t\t(%12.10g, %12.10g) regular gamut range\n", gmt[0], gmt[1]);
		
		free_vector(idx);
		free_vector(rnk);
	}
	
	free_vector(eq);
	free_vector(dg);
	free_vector(x);
	free_vector(y);
	
	return k != n;
	
	#undef EEQUIV
	#undef EQUIV
}

/* Read curve data */
int read_curve(FILE *fp, int io, int curve, double ***m, double ***im)
{
	int i;
	char *p, *np;
	int n = 0, size = 64;
	size_t bsize = 128;
	char *buffer = (char *)xmalloc(bsize);
	
	
	*m = matrix(4, size);
	
	/* Read curve data */
	while (getline(&buffer, &bsize, fp) != -1) {
		if (*buffer == '#') continue;
		if (*buffer == ';') break;
		
		if (n >= size)
			*m = grow_matrix(*m, 4, size<<=1);
		
		(*m)[0][n] = strtod(buffer, &p);
		(*m)[1][n] = strtod(p, &np);
		
		if ((p == buffer) || (np == p))
			error("syntax error in curve data");
		
		n++;
	}
	
	/* Check if we got any data */
	if (!n) error("No curve data supplied");
	
	/* Check for degeneracy */
	if (degenerate(&n, *m, io, gamut_constraint(io == 0 ? ins : outs, curve,
						    io == 0 ? ins_gamut : outs_gamut))) {
		warning("Degenerate calibration data, doing my best to cope...");
		
		(io == 0 ? ins_gamut_constrained : outs_gamut_constrained) = 1;
	}
	
	/* Prepare curve data for interpolation */
	*im = matrix(4, size);
	
	for (i = 0; i < n; i++) {
		(*im)[0][i] = (*m)[1][i];
		(*im)[1][i] = (*m)[0][i];
	}
	
	if (verbose > 1) {
		if (verbose > 3)
		    for (i = 0; i < n; i++)
			fprintf(stderr, "\t\t%12.10g %12.10g\n", (*m)[0][i], (*m)[1][i]);
		
		fprintf(stderr, "\t%i points, prepping for interpolation...", n);
		fflush(stderr);
	}
	
	derivs(*m, n); derivs(*im, n);
	
	if (verbose > 1) fprintf(stderr, " done.\n");
	free(buffer);
	
	return n;
}

/* Read LUT data */
void read_lut(FILE *fp, void **part, void **ipart)
{
	int i, j;
	int n = 0, size = 64, cdp = 0, cdsize = 64;
	double **m = matrix(6, size), **im, *cdv = vector(cdsize);
	
	size_t bsize = 128;
	char *label, *p, *np;
	char *buffer = (char *)xmalloc(bsize);
	double ip[MAXCHANNELS], op[MAXCHANNELS], tp[MAXCHANNELS];
	double e = 0.0, dE[3] = {0.0 /* max */, 1.0e12 /* min */, 0.0 /* avg */};
	
	
	/* Read LUT data */
	while (getline(&buffer, &bsize, fp) != -1) {
		if (*buffer == '#') continue;
		if (*buffer == ';') break;
		if (ignore_lut) continue;
		
		/* Read data */
		sscanf(buffer, "%as %n", &label, &i);
		p = buffer + i; np = NULL;
		
		for (i = 0; i < ins_channels; i++) {
			ip[i] = strtod(p, &np);
			if (p == np)
				error("syntax error in lut data");
			p = np;
			
			/* Tuck a copy away for later */
			if (cdp >= cdsize)
				cdv = grow_vector(cdv, cdsize <<= 1);
			cdv[cdp++] = ip[i];
		}
		
		for (i = 0; i < outs_channels; i++) {
			op[i] = strtod(p, &np);
			if (p == np)
				error("syntax error in lut data");
			p = np;
			
			/* Tuck a copy away for later */
			if (cdp >= cdsize)
				cdv = grow_vector(cdv, cdsize <<= 1);
			cdv[cdp++] = op[i];
		}
		
		/* Check that the data point is inside gamut constraints */
		if ((ins_gamut_constrained && !in_gamut(ins, ins_gamut, ip)) ||
		    (outs_gamut_constrained && !in_gamut(outs, outs_gamut, op))) {
			if (verbose > 3)
				fprintf(stderr, "%20s  outside of gamut bounds, skipping!\n", label);
			continue;
		}
		
		/* Process LUT data point */
		if (n+4 >= size)
			m = grow_matrix(m, 6, size <<= 1);
		
		/* push input data forward to XYZ */
		incurves(NULL, ip, ip); (*ins2XYZ)(ip, ip);
		m[0][n] = ip[0]; m[1][n] = ip[1]; m[2][n] = ip[2];
		
		/* pull output data back to XYZ */
		outcurves_1(NULL, op, op); (*outs2XYZ)(op, op);
		m[3][n] = op[0]; m[4][n] = op[1]; m[5][n] = op[2];
		
		if (verbose > 3)
			fprintf(stderr, "%20s %12.10g %12.10g %12.10g -> %12.10g %12.10g %12.10g\n",
				label, m[0][n], m[1][n], m[2][n], m[3][n], m[4][n], m[5][n]);
		
		free(label);
		
		n++;
	}
	
	/* Build inverse data if required */
	if (invertible) {
		im = matrix(6, size);
		
		for (j = 0; j < n; j++) {
			im[0][j] = m[3][j];
			im[1][j] = m[4][j];
			im[2][j] = m[5][j];
			
			im[3][j] = m[0][j];
			im[4][j] = m[1][j];
			im[5][j] = m[2][j];
		}
	}
	
	/* Tuck a copy of calibration data away for later */
	calibration_data_pts = cdp;
	calibration_data_vector = cdv;
	
	/* Fit calibration data and prepare for interpolation */
	if (!ignore_lut && n) {
		/* Find best linear RGB transform */
		best_linear_fit(m, n, L); inv33(L, L1);
		
		if (verbose > 4) {
			fprintf(stderr, "\tLinear transform prefactor:\n");
			fprintf(stderr, "\t\t    [ %13.10g %13.10g %13.10g ]\n", L[0][0], L[0][1], L[0][2]);
			fprintf(stderr, "\t\tL = [ %13.10g %13.10g %13.10g ]\n", L[1][0], L[1][1], L[1][2]);
			fprintf(stderr, "\t\t    [ %13.10g %13.10g %13.10g ]\n", L[2][0], L[2][1], L[2][2]);
		}
		
		/* If enough data, find better non-linear fit */
		if (n > 64) {
			use_higher_order_approximation = 1;
			P = best_poly_fit(m, n);
			if (invertible) P1 = best_poly_fit(im, n);
			
			if (verbose > 3)
				fprintf(stderr, "\tHave enough data, using better non-linear fit...\n");
			
			if (verbose > 4) {
				fprintf(stderr, "\tNon-linear approximation coefficients:\n");
				for (i = 0; i < 19; i++)
					fprintf(stderr, "\t\t    [ %13.10g %13.10g %13.10g ]\n", P[0][i], P[1][i], P[2][i]);
			}
		}
		
		/* Prepare remainder for interpolation */
		for (j = 0; j < n; j++) {
			tp[0] = m[0][j];
			tp[1] = m[1][j];
			tp[2] = m[2][j];
			
			op[0] = m[3][j];
			op[1] = m[4][j];
			op[2] = m[5][j];
			
			if (use_higher_order_approximation)
				poly_approx(P, tp, ip); else apply33(L, tp, ip);
			
			m[0][j] = ip[0];
			m[1][j] = ip[1];
			m[2][j] = ip[2];
			
			m[3][j] -= m[0][j];
			m[4][j] -= m[1][j];
			m[5][j] -= m[2][j];
			
			e = XYZ_dE(ip, op);
			if (e > dE[0]) dE[0] = e;
			if (e < dE[1]) dE[1] = e;
			dE[2] += e;
		}
		
		if (verbose > 4)
			fprintf(stderr, "\tResidual error to correct by multi-dimensional interpolation:\n"
					"\t\tdE = (%6.4g max, %6.4g min, %6.4g avg)\n", dE[0], dE[1], dE[2]/n);
		
		if (verbose > 1) {
			fprintf(stderr, "\t%i points, prepping for interpolation...", n);
			fflush(stderr);
		}
		
		*part = subdivide(m, n);
		
		/* Do the same for inverted data if required */
		if (invertible) {
			for (j = 0; j < n; j++) {
				tp[0] = im[0][j];
				tp[1] = im[1][j];
				tp[2] = im[2][j];
				
				if (use_higher_order_approximation)
					poly_approx(P1, tp, ip); else apply33(L1, tp, ip);
				
				im[0][j] = ip[0];
				im[1][j] = ip[1];
				im[2][j] = ip[2];
				
				im[3][j] -= im[0][j];
				im[4][j] -= im[1][j];
				im[5][j] -= im[2][j];
			}
			
			*ipart = subdivide(im, n);
		}
		
		if (verbose > 1) fprintf(stderr, " done.\n");
	} else {
		free_matrix(m); if (invertible) free_matrix(im);
		if (verbose > 1) fprintf(stderr, " ignored.\n");
	}
	
	free(buffer);
}


/* Read calibration data */
void read_calibration_data(FILE *fp)
{
	int i;
	size_t bsize = 128;
	char *p, *buffer = (char *)xmalloc(bsize);
	
	if (verbose) fprintf(stderr, "Reading calibration data...\n");
	
	init_gamut_constraint(ins, ins_gamut);
	init_gamut_constraint(outs, outs_gamut);
	
	while (getline(&buffer, &bsize, fp) != -1) {
		buffer[strlen(buffer)-1] = 0;
		
		if (*buffer == '#') continue;
		
		if ((p = strstr(buffer, "WHITEPOINT:"))) {
			char *p1 = p + 11, *p2, *p3, *np;
			
			if (verbose > 1) fprintf(stderr, "\t%s\n", buffer);
			
			media_white_pt[0] = strtod(p1, &p2);
			media_white_pt[1] = strtod(p2, &p3);
			media_white_pt[2] = strtod(p3, &np);
			
			if ((p2 == p1) || (p3 == p2) || (np == p3))
				error("syntax error in whitepoint data");
			
		} else if ((p = strstr(buffer, "BLACKPOINT:"))) {
			char *p1 = p + 11, *p2, *p3, *np;
			
			if (verbose > 1) fprintf(stderr, "\t%s\n", buffer);
			
			media_black_pt[0] = strtod(p1, &p2);
			media_black_pt[1] = strtod(p2, &p3);
			media_black_pt[2] = strtod(p3, &np);
			
			if ((p2 == p1) || (p3 == p2) || (np == p3))
				error("syntax error in blackpoint data");
			
		} else if (sscanf(buffer, "CURVE IN%i", &i) == 1) {
			if (verbose > 1) {
				fprintf(stderr, "\t%s", buffer);
				if (verbose > 3) fprintf(stderr, "\n");
				else fflush(stderr);
			}
			
			ins_curves_n[i] = read_curve(fp, 0 /*in*/, i, &(ins_curves[i]), &(ins_curves_1[i]));
			use_ins_curves |= 1 << i;
			
		} else if (sscanf(buffer, "CURVE OUT%i", &i) == 1) {
			if (verbose > 1) {
				fprintf(stderr, "\t%s", buffer); 
				if (verbose > 3) fprintf(stderr, "\n");
				else fflush(stderr);
			}
			
			outs_curves_n[i] = read_curve(fp, 1 /*out*/, i, &(outs_curves[i]), &(outs_curves_1[i]));
			use_outs_curves |= 1 << i;
			
		} else if (strstr(buffer, "LUT:")) {
			if ((use_ins_curves && (use_ins_curves != (1 << ins_channels) - 1)) ||
			    (use_outs_curves && (use_outs_curves != (1 << outs_channels) - 1)) )
				error("All curves must be specified before trying to read in LUT");
			
			if (verbose > 1) {
				fprintf(stderr, "\t%s", buffer);
				if (verbose > 3 && !ignore_lut) fprintf(stderr, "\n");
				else fflush(stderr);
			}
			
			read_lut(fp, &ELUT, &ELUT_1);
		}
	}
	
	if ((use_ins_curves && (use_ins_curves != (1 << ins_channels) - 1)) ||
	    (use_outs_curves && (use_outs_curves != (1 << outs_channels) - 1)) )
		error("Some curves are missing from calibration data");
	
	if (use_ins_curves && (ins_gamma != 1.0))
		warning("Choosing supplied input curves in favor of gamma correction");
	if (use_outs_curves && (outs_gamma != 1.0))
		warning("Choosing supplied output curves in favor of gamma correction");
	
	free(buffer);
}



/******************* Profile builder **********************************/

/* Clip value to specified range */
static double clip(double v, double min, double max, int *clipped)
{
	if (v < min) { *clipped |= 0x01; return min; }
	if (v > max) { *clipped |= 0x01; return max; }
	
	return v;
}


/* Create a matrix or LUT-based profile */
void build_profile(char *file)
{
	icc *icco = new_icc();
	FILE *fp = xfopen(file, "wb");
	
	if (verbose) fprintf(stderr, "Creating [%s]->[%s] %s-based %s device profile '%s'\n",
		ColorSpaceSignature2str(ins), ColorSpaceSignature2str(outs),
		matrix_only ? "matrix" : "LUT", class == icSigInputClass ? "input" :
		(class == icSigDisplayClass ? "display" : "output"), file);
	
	if (!icco) error ("Creation of ICC object failed");
	
	
	/* Add all the tags required */
	
	/* The header: */
	{
		icmHeader *h = icco->header;
		
		/* Profile creator */
		h->cmmId = h->creator = str2tag("scrs");
		
		/* Values that must be set before writing */
		h->deviceClass     = class;
    		h->colorSpace      = ins;
    		h->pcs             = outs;
		
		/* icclib seems somewhat confused on the relative vs. absolute issue,
		   so to be on the safe side, rendering intent is set to perceptual */
    		h->renderingIntent = icPerceptual;
		
		/* Values that should be set before writing */
		h->manufacturer = str2tag(manufacturer);
    		h->model        = str2tag(model);
	}
	
	/* Profile Description Tag: */
	{
		icmTextDescription *dt = (icmTextDescription *)icco->add_tag(
			icco, icSigProfileDescriptionTag, icSigTextDescriptionType);
		
		if (!dt) error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		dt->size = strlen(description)+1;	/* Allocated and used size of desc, inc null */
		dt->allocate((icmBase *)dt);		/* Allocate space */
		strcpy(dt->desc, description);		/* Copy the string in */
	}
	
	/* Copyright Tag: */
	{
		icmText *ct = (icmText *)icco->add_tag(icco, icSigCopyrightTag, icSigTextType);
		
		if (!ct) error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		ct->size = strlen(copyright)+1;		/* Allocated and used size of text, inc null */
		ct->allocate((icmBase *)ct);		/* Allocate space */
		strcpy(ct->data, copyright);		/* Copy the text in */
	}
	
	/* Media White Point Tag: */
	{
		icmXYZArray *wpt = (icmXYZArray *)icco->add_tag(
			icco, icSigMediaWhitePointTag, icSigXYZArrayType);
		/* Note that tag types icSigXYZType and icSigXYZArrayType are identical */
		
		if (!wpt) error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		wpt->size = 1;
		wpt->allocate((icmBase *)wpt);	/* Allocate space */
		wpt->data[0].X = media_white_pt[0];
		wpt->data[0].Y = media_white_pt[1];
		wpt->data[0].Z = media_white_pt[2];
	}
	
	/* Media Black Point Tag: */
	{
		icmXYZArray *bpt = (icmXYZArray *)icco->add_tag(
			icco, icSigMediaBlackPointTag, icSigXYZArrayType);
		
		if (!bpt) error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		bpt->size = 1;
		bpt->allocate((icmBase *)bpt);	/* Allocate space */
		bpt->data[0].X = media_black_pt[0];
		bpt->data[0].Y = media_black_pt[1];
		bpt->data[0].Z = media_black_pt[2];
	}
	
	
	/******************* Matrix-based profile *********************/
	
	/* Red, Green and Blue Colorant Tags: */
	{
		double M[3][3];
		
		icmXYZArray *r = (icmXYZArray *)icco->add_tag(
			icco, icSigRedColorantTag, icSigXYZArrayType);
		icmXYZArray *g = (icmXYZArray *)icco->add_tag(
			icco, icSigGreenColorantTag, icSigXYZArrayType);
		icmXYZArray *b = (icmXYZArray *)icco->add_tag(
			icco, icSigBlueColorantTag, icSigXYZArrayType);
		
		if (!r || !g || !b)
			error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		r->size = g->size = b->size = 1;
		
		r->allocate((icmBase *)r);	/* Allocate space */
		g->allocate((icmBase *)g);
		b->allocate((icmBase *)b);
		
		mult33(L, M_RGB2XYZ, M);
		
		r->data[0].X = M[0][0]; r->data[0].Y = M[1][0]; r->data[0].Z = M[2][0];
		g->data[0].X = M[0][1]; g->data[0].Y = M[1][1]; g->data[0].Z = M[2][1];
		b->data[0].X = M[0][2]; b->data[0].Y = M[1][2]; b->data[0].Z = M[2][2];
	}
	
	/* Red, Green and Blue Tone Reproduction Curve Tags: */
	{
		int i, rv = 0;
		double p[MAXCHANNELS];
		
		icmCurve *r = (icmCurve *)icco->add_tag(
				icco, icSigRedTRCTag, icSigCurveType);
		icmCurve *g = (icmCurve *)icco->add_tag(
				icco, icSigGreenTRCTag, icSigCurveType);
		icmCurve *b = (icmCurve *)icco->add_tag(
				icco, icSigBlueTRCTag, icSigCurveType);
		
		if (!r || !g || !b)
			error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		r->flag = g->flag = b->flag = icmCurveSpec; 	/* Specified version */
		r->size = g->size = b->size = lut_size_1d;	/* Number of entries (min must be 2!) */
		
		r->allocate((icmBase *)r);	/* Allocate space */
		g->allocate((icmBase *)g);
		b->allocate((icmBase *)b);
		
		for (i = 0; i < lut_size_1d; i++) {
			p[0] = p[1] = p[2] = i/(lut_size_1d-1.0);
			
			incurves(NULL, p, p);	/* Transfer function */
			r->data[i] = clip(p[0], 0.0, 1.0, &rv);
			g->data[i] = clip(p[1], 0.0, 1.0, &rv);
			b->data[i] = clip(p[2], 0.0, 1.0, &rv);
		}
		
		if (verbose > 2 && rv)
			warning("%s: warning: RGB curves were clipped", file);
	}
	
	
	/********************** LUT-based profile *********************/
	
	/* dev -> pcs lut: */
	if (!matrix_only) {
		/* Intent 0 = perceptual */
		icmLut *lut = (icmLut *)icco->add_tag(icco, icSigAToB0Tag,
			bitspersample == 16 ? icSigLut16Type : icSigLut8Type);
		
		if (!lut) error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		lut->inputChan = channels(ins); lut->outputChan = channels(outs);
		lut->inputEnt = lut->outputEnt = lut_size_1d;
		
		lut->clutPoints = lut_size_3d;
		
		lut->allocate((icmBase *)lut);	/* Allocate space */
		
		/* The matrix is only applicable to XYZ input space, */
		/* so it is not used here. */
		
		/* Use helper function to do the hard work. */
		if (lut->set_tables(lut, NULL,
				ins,			/* Input color space */
				outs,			/* Output color space */
				incurves_expanded,	/* Input transfer function (NULL = default) */
				NULL, NULL,		/* Use default Maximum range of values */
				ctransform,		/* LUT transform function */
				NULL, NULL,		/* Use default Maximum range of values */
				outcurves		/* Output transfer function */
		) != 0)
			error("Setting %i bit [%s]->[%s] Lut failed: %d, %s",
				bitspersample, ColorSpaceSignature2str(ins),  ColorSpaceSignature2str(outs),
				icco->errc, icco->err);
	}
	
	/* dev -> pcs lut: link intent 1 to intent 0 */
	if (!matrix_only) {
		/* Intent 1 = relative colorimetric */
		icmLut *lut = (icmLut *)icco->link_tag(icco, icSigAToB1Tag, icSigAToB0Tag);
		
		if (!lut) error("link_tag failed: %d, %s", icco->errc, icco->err);
	}
	
	/* dev -> pcs lut: link intent 2 to intent 0 */
	if (!matrix_only) {
		/* Intent 2 = saturation */
		icmLut *lut = (icmLut *)icco->link_tag(icco, icSigAToB2Tag, icSigAToB0Tag);
		
		if (!lut) error("link_tag failed: %d, %s", icco->errc, icco->err);
	}
	
	/* pcs -> dev lut: */
	if (!matrix_only && invertible) {
		/* Intent 0 = perceptual */
		icmLut *lut = (icmLut *)icco->add_tag(icco, icSigBToA0Tag,
			bitspersample == 16 ? icSigLut16Type : icSigLut8Type);
		
		if (!lut) error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		lut->inputChan = channels(outs); lut->outputChan = channels(ins);
		lut->inputEnt = lut->outputEnt = lut_size_1d;
		
		lut->clutPoints = lut_size_3d;
		
		lut->allocate((icmBase *)lut);	/* Allocate space */
		
		/* The matrix is only applicable to XYZ input space, */
		/* so it is not used here. */
		
		/* Use helper function to do the hard work. */
		if (lut->set_tables(lut, NULL,
				outs,			/* Input color space */
				ins,			/* Output color space */
				outcurves_1_expanded,	/* Input transfer function */
				NULL, NULL,		/* Use default range */
				ctransform_1,		/* LUT transform function */
				NULL, NULL,		/* Use default range */
				incurves_1		/* Output transfer function (NULL = deflt) */
		) != 0)
			error("Setting %i bit [%s]->[%s] Lut failed: %d, %s",
				bitspersample, ColorSpaceSignature2str(outs),  ColorSpaceSignature2str(ins),
				icco->errc, icco->err);
	}
	
	/* pcs -> dev lut: link intent 1 to intent 0 */
	if (!matrix_only && invertible) {
		/* Intent 1 = relative colorimetric */
		icmLut *lut = (icmLut *)icco->link_tag(icco, icSigBToA1Tag, icSigBToA0Tag);
		
		if (!lut) error("link_tag failed: %d, %s", icco->errc, icco->err);
	}
	
	/* pcs -> dev lut: link intent 2 to intent 0 */
	if (!matrix_only && invertible) {
		/* Intent 2 = saturation */
		icmLut *lut = (icmLut *)icco->link_tag(icco, icSigBToA2Tag, icSigBToA0Tag);
		
		if (!lut) error("link_tag failed: %d, %s", icco->errc, icco->err);
	}
	
	
	/* Gamut tag is not generated (yet?) - AF */
	
	
	/* Set profile size in the header */
	icco->header->size = icco->get_size(icco);
	
	/* Dump ICC profile header */
	if (verbose > 2) {
		fprintf(stderr, "%s ", file);
		icco->header->dump(icco->header, stderr, 1);
	}
	
	/* Write the file out */
	if (icco->write(icco, fp, 0) != 0)
		error("Write file failed: %s", icco->err);
	
	icco->free(icco);
	fclose(fp);
}


/* Test profile */
void test_profile(char *file)
{
	icc *icco = new_icc();
	FILE *fp = xfopen(file, "rb");
	
	
	if (verbose) fprintf(stderr, "Testing profile... ");
	
	if (!icco) error ("Creation of ICC object failed");
	
	if (icco->read(icco, fp, 0))
		error("Read file failed: %s", icco->err);
	
	
	/* Test forward lookup */
	if (calibration_data_pts) {
		int i = 0, n = 0;
		double ip[3], op[3], tp[MAXCHANNELS];
		double e = 0.0, dE[3] = {-HUGE /* max */, HUGE /* min */, 0.0 /* avg */};
		
		/* Get a conversion object */
		icmLuBase *luo = icco->get_luobj(icco, icmFwd, icmDefaultIntent, icmLuOrdNorm);
		if (!luo) error("%s: Error %d, %s", file, icco->errc, icco->err);
		
		for (i = 0; i < calibration_data_pts; i += ins_channels+outs_channels) {
			if (luo->lookup(luo, tp, &(calibration_data_vector[i])) > 1)
				error("%s: Error %d, %s", file, icco->errc, icco->err);
			
			(*outs2XYZ)(tp, ip);
			(*outs2XYZ)(&(calibration_data_vector[i+ins_channels]), op);
			
			e = XYZ_dE(ip, op);
			if (e > dE[0]) dE[0] = e;
			if (e < dE[1]) dE[1] = e;
			dE[2] += e; n++;
		}
		
		if (verbose) fprintf(stderr, "\n\t Forward lookup dE = (%6.4g max, %6.4g min, %6.4g avg)\n", dE[0], dE[1], dE[2]/n);
		
		luo->free(luo);
	}
	
	/* Test backward lookup */
	if (calibration_data_pts && invertible) {
		int i = 0, n = 0;
		double ip[3], op[3], tp[MAXCHANNELS];
		double e = 0.0, dE[3] = {-HUGE /* max */, HUGE /* min */, 0.0 /* avg */};
		
		/* Get a conversion object */
		icmLuBase *luo = icco->get_luobj(icco, icmBwd, icmDefaultIntent, icmLuOrdNorm);
		if (!luo) error("%s: Error %d, %s", file, icco->errc, icco->err);
		
		for (i = 0; i < calibration_data_pts; i += ins_channels+outs_channels) {
			if (luo->lookup(luo, tp, &(calibration_data_vector[i+ins_channels])) > 1)
				error("%s: Error %d, %s", file, icco->errc, icco->err);
			
			(*ins2XYZ)(tp, ip);
			(*ins2XYZ)(&(calibration_data_vector[i]), op);
			
			e = XYZ_dE(ip, op);
			if (e > dE[0]) dE[0] = e;
			if (e < dE[1]) dE[1] = e;
			dE[2] += e; n++;
		}
		
		if (verbose) fprintf(stderr, "\tBackward lookup dE = (%6.4g max, %6.4g min, %6.4g avg)\n", dE[0], dE[1], dE[2]/n);
		
		luo->free(luo);
	}
	
	if (verbose) fprintf(stderr, "OK.\n");
	
	
	icco->free(icco);
	fclose(fp);
}



/******************* Main routine *************************************/

int main(int argc, char *argv[])
{
	char c;
	
	/* Set default primaries */
	LookupPrimaries("Adobe", primaries);
	LookupPrimaries("D50", primaries);
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv-:c:i:o:p:C:UE:MLm:d:r:b:s:")) != -1)
	switch (c) {
	/* General options */
		case 'h':				/* Help message */
			usage(); break;
		case 'v':				/* Verbosity */
			verbose++; break;
		case '-':				/* Long options */
			readopt(&argc, &argv, SELF ".options", optarg);
			break;
	
	/* Color space options */
		case 'c':				/* Profile class */
			switch (optarg[0]) {
				case 'i':
					class = icSigInputClass;
					ins = icSigRgbData;
					outs = icSigXYZData;
					break;
				case 'd':
					class = icSigDisplayClass;
					ins = icSigRgbData;
					outs = icSigXYZData;
					break;
				case 'o':
					class = icSigOutputClass;
					ins = icSigRgbData;
					outs = icSigXYZData;
					break;
				case 'c':
					class = icSigColorSpaceClass;
					ins = icSigXYZData;
					outs = icSigLabData;
					break;
				case 'a':
					class = icSigAbstractClass;
					ins = icSigXYZData;
					outs = icSigXYZData;
					break;
				default:
					usage();
			}
			break;
		case 'i':				/* Input color space */
			{
				char *s = xstrdup(optarg);
				char *c = strchr(s, ':');
				
				if (c) { *c = 0; ins_gamma = strtod(++c, NULL); }
				ins = str2ColorSpaceSignature(s);
				
				free(s);
			}
			break;
		case 'o':				/* Output color space */
			{
				char *s = xstrdup(optarg);
				char *c = strchr(s, ':');
				
				if (c) { *c = 0; outs_gamma = strtod(++c, NULL); }
				outs = str2ColorSpaceSignature(s);
				
				free(s);
			}
			break;
		case 'p':				/* Primaries */
			if (!LookupPrimaries(optarg, primaries)) {
				double x, y;
				char *s = strchr(optarg, ':');
				
				if (!s) usage();
				x = strtod(++s, &s);
				if (*s != ',') usage();
				y = strtod(++s, NULL);
				
				switch (optarg[0]) {
					case 'w':
						primaries[0][0] = x;
						primaries[0][1] = y;
						break;
					case 'r':
						primaries[1][0] = x;
						primaries[1][1] = y;
						break;
					case 'g':
						primaries[2][0] = x;
						primaries[2][1] = y;
						break;
					case 'b':
						primaries[3][0] = x;
						primaries[3][1] = y;
						break;
					default:
						usage();
				}
			}
			break;
	
	/* Color mapping algorithm options */
		case 'C':				/* Calibration data */
			calibration_data = xfetch(NULL, optarg, "r");
			break;
		case 'U':				/* Unidirectional */
			invertible = 0;
			break;
		case 'E':				/* LUT shadow expansion */
			lut_shadow_expansion = strtod(optarg, NULL);
			if (lut_shadow_expansion == 0.0) usage();
			break;
		case 'M':				/* Matrix profile only */
			matrix_only = 1;
			break;
		case 'L':				/* Ignore LUT */
			ignore_lut = 1;
			break;
	
	/* ICC profile options */
		case 'm':				/* Manufacturer & model */
			{
				char *s = xstrdup(optarg);
				char *c = strchr(s, ':');
				
				if (c) { *c = 0; model = ++c; }
				manufacturer = s;
			}
			break;
		case 'd':				/* Description */
			description = optarg;
			break;
		case 'r':				/* Copyright */
			copyright = optarg;
			break;
		case 'b':				/* LUT bits per sample */
			bitspersample = atoi(optarg);
			if (bitspersample != 8 && bitspersample != 16) usage();
			break;
		case 's':				/* LUT size */
			if (strchr(optarg, ':'))
				sscanf(optarg, "%i:%i", &lut_size_1d, &lut_size_3d);
			else
				lut_size_1d = lut_size_3d = atoi(optarg);
			
			if (lut_size_1d < 3 || lut_size_3d < 3) usage ();
			break;
	
	/* Default response */
		default:
			usage();
	}
	
	if (argc != optind+1) usage();
	profile = argv[optind++];
	
	
	/* Initialize transform data */
	if (verbose) fprintf(stderr, "%s\n", usage_msg[0]);
	
	SetPrimaries(primaries);
	
	ins_channels = channels(ins);
	outs_channels = channels(outs);
	
	ins2XYZ = toXYZ(ins); XYZ2ins = fromXYZ(ins);
	outs2XYZ = toXYZ(outs); XYZ2outs = fromXYZ(outs);
	
	
	/* Read in calibration data if available */
	if (calibration_data)
		read_calibration_data(calibration_data);
	
	
	/* Build profile */
	build_profile(profile);
	
	/* Test profile accuracy */
	test_profile(profile);
	
	return 0;
}
