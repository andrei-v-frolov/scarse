/* $Id: ipb.c,v 1.26 2005/10/20 06:14:58 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * ICC profile builder - build profile based on measured data.
 * 
 * Copyright (C) 1999-2005 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <frolov@cita.utoronto.ca>
 * 
 */

/* CREDITS:
 *   Derived from icclib & examples by Graeme W. Gill
 */

/* TODO:
 *   - Convert to robust fitting (work in progress...)
 *   - Generalize variance handling (input class only now)
 *   - Sort out the mess of input/output vs device/PCS notation
 *   - Do something about gamut boundaries and tables...
 */

#include "scarse.h"
#define SELF "ipb"



/******************* Options and defaults *****************************/

/* Usage */
char *program_name = SELF;
char *usage_msg[] = {
	"ICC profile builder, Version " PACKAGE_VERSION,
	"Author: " PACKAGE_BUGREPORT,
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
	"		     RGB spaces: Adobe (default), Apple, Best, Beta, Bruce, CIE,",
	"		                 ColorMatch, ECI, EktaSpace, NTSC, PAL/SECAM,",
	"		                 ProPhoto, SMPTE-C, sRGB (simplified), WideGamut",
	"",
	" Color mapping algorithm options:",
	"  -C file	calibration data; '-' means read from stdin",
	"		required format of the data is described elsewhere",
	"  -U		generate unidirectional (device -> PCS) profile only",
	"  -E gamma	LUT shadow expansion exponent (none=1.0, default=3.0)",
	"  -M		generate matrix based profile only, do not output LUT",
	"  -I pattern	ignore patches matching pattern; 'none' forces all to be used",
	"",
	" ICC profile options:",
	"  -m manufacturer[:model] of device for which profile is created",
	"  -d string	description of device profile (required by Photoshop)",
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

static FILE *calibration_data_stream = NULL;
static int                invertible = 1;
static int                ignore_lut = 0;
static int               ignore_none = 0;
static char          *ignore_pattern = NULL;
static int               matrix_only = 0;
static double   lut_shadow_expansion = 3.0;

static char            *manufacturer = "none";
static char                   *model = "none";
static char             *description = NULL;
static char               *copyright = "";

static int             bitspersample = 16;
static int               lut_size_1d = 1024;
static int               lut_size_3d = 33;

static char                 *profile = NULL;



/******************* Transform data model *****************************/

/* General transform info */
static int              ins_channels = 0;
static double              ins_gamma = 1.0;
static transform             ins2XYZ = NULL,
                             XYZ2ins = NULL;

static int             outs_channels = 0;
static double             outs_gamma = 1.0;
static transform            outs2XYZ = NULL,
                            XYZ2outs = NULL;


/* Curve data and models */
typedef struct {
	int n;		/* data points */
	double **data;	/* [y,x,dx][n] */
	double *fit;	/* curve fit */
} curve;

static int	      use_ins_curves = 0,
		     use_outs_curves = 0;

static curve ins_curves[MAXCHANNELS],
            outs_curves[MAXCHANNELS];


/* LUT data and models */
typedef struct {
	int flag;			/* outlier? */
	char *label;			/* patch label */
	/* Measured values */
	double  in[MAXCHANNELS];	/* input */
	double out[MAXCHANNELS];	/* output */
	double var[MAXCHANNELS];	/* variance */
	/* Curve pullbacks */
	double DEV[MAXCHANNELS];	/* (linearized) device */
	double XYZ[3];			/* PCS */
} datapt;

static int           calibration_pts = 0;
static datapt      *calibration_data = NULL;

static int                  use_poly = 0;
static double                    **P = NULL,
                                **P1 = NULL;

static int                use_interp = 0;
static int                       PTS = 0;
static double                    **T = NULL,  Q[3][3];
static double                   **T1 = NULL, Q1[3][3];



/******************* Transform functions ******************************/

/* Soft clip transfer function */
static double sclip(double x, double min, double max)
{
	if (x < min) { return min*min/(2.0*min - x); }
	if (x > max) { return (x - max*max)/(x - 2.0*max + 1.0); }
	
	return x;
}

/* Translate data through curves */
void xlate_curve(double in[], double out[], int channels, curve c[])
{
	int i;
	
	for (i = 0; i < channels; i++)
		out[i] = sclip(lu_curve(c[i].fit, in[i]), RMIN, RMAX);
}

/* Translate data through reverse curves */
void xlate_curve_1(double in[], double out[], int channels, curve c[])
{
	int i;
	
	for (i = 0; i < channels; i++)
		out[i] = lu_curve_1(c[i].fit, in[i]);
}


/* Input curves */
void incurves(void *cntx, double out[], double in[])
{
	if (use_ins_curves)
		xlate_curve(in, out, ins_channels, ins_curves);
	else if (ins_gamma != 1.0)
		vgamma(in, out, ins_channels, ins_gamma);
	else
		vcopy(in, out, ins_channels);
}

/* Input curves for shadow-expanded LUT */
void incurves_expanded(void *cntx, double out[], double in[])
{
	incurves(cntx, out, in);
	
	if (lut_shadow_expansion != 1.0)
		vgamma(out, out, ins_channels, 1.0/lut_shadow_expansion);
}

/* Inverse of input curves */
void incurves_1(void *cntx, double out[], double in[])
{
	if (use_ins_curves)
		xlate_curve_1(in, out, ins_channels, ins_curves);
	else if (ins_gamma != 1.0)
		vgamma(in, out, ins_channels, 1.0/ins_gamma);
	else
		vcopy(in, out, ins_channels);
}

/* Output curves */
void outcurves(void *cntx, double out[], double in[])
{
	if (use_outs_curves)
		xlate_curve_1(in, out, outs_channels, outs_curves);
	else if (outs_gamma != 1.0)
		vgamma(in, out, outs_channels, 1.0/outs_gamma);
	else
		vcopy(in, out, outs_channels);
}

/* Inverse of output curves */
void outcurves_1(void *cntx, double out[], double in[])
{
	if (use_outs_curves)
		xlate_curve(in, out, outs_channels, outs_curves);
	else if (outs_gamma != 1.0)
		vgamma(in, out, outs_channels, outs_gamma);
	else
		vcopy(in, out, outs_channels);
}

/* Inverse of output curves for shadow-expanded LUT */
void outcurves_1_expanded(void *cntx, double out[], double in[])
{
	outcurves_1(cntx, out, in);
	
	if (lut_shadow_expansion != 1.0)
		vgamma(out, out, ins_channels, 1.0/lut_shadow_expansion);
}


/* 3D -> 3D color transformation */
void ctransform(void *cntx, double out[], double in[])
{
	double XYZ[3], Lab[3], t[3];
	
	if (lut_shadow_expansion != 1.0)
		vgamma(in, in, ins_channels, lut_shadow_expansion);
		/* this clobbers input, but we don't care... */
	
	(*ins2XYZ)(in, XYZ);
	
	if (use_poly) evalf_poly(P, XYZ, XYZ);
	
	if (use_interp) {
		XYZ2Lab(XYZ, Lab); interp3d(T, PTS, Q, Lab, t);
		XYZ[0] += t[0]; XYZ[1] += t[1]; XYZ[2] += t[2];
	}
	
	XYZ[0] = sclip(XYZ[0], RMIN*XYZ_ILLUM[0], RMAX*XYZ_ILLUM[0]);
	XYZ[1] = sclip(XYZ[1], RMIN*XYZ_ILLUM[1], RMAX*XYZ_ILLUM[1]);
	XYZ[2] = sclip(XYZ[2], RMIN*XYZ_ILLUM[2], RMAX*XYZ_ILLUM[2]);
	
	(*XYZ2outs)(XYZ, out);
}

/* Inverse 3D -> 3D color transformation */
void ctransform_1(void *cntx, double out[], double in[])
{
	double XYZ[3], Lab[3], t[3];
	
	if (lut_shadow_expansion != 1.0)
		vgamma(in, in, ins_channels, lut_shadow_expansion);
		/* this clobbers input, but we don't care... */
	
	(*outs2XYZ)(in, XYZ);
	
	if (use_poly) evalf_poly(P1, XYZ, XYZ);
	
	if (use_interp) {
		XYZ2Lab(XYZ, Lab); interp3d(T1, PTS, Q1, Lab, t);
		XYZ[0] += t[0]; XYZ[1] += t[1]; XYZ[2] += t[2];
	}
	
	XYZ[0] = sclip(XYZ[0], RMIN*XYZ_ILLUM[0], RMAX*XYZ_ILLUM[0]);
	XYZ[1] = sclip(XYZ[1], RMIN*XYZ_ILLUM[1], RMAX*XYZ_ILLUM[1]);
	XYZ[2] = sclip(XYZ[2], RMIN*XYZ_ILLUM[2], RMAX*XYZ_ILLUM[2]);
	
	(*XYZ2ins)(XYZ, out);
}



/******************* Calibration data handling ************************/

/* small deviations are suspicious; limit them from below */
#define LIMDEV(x,dx) (sqrt(1.0e-6 + 1.0e-4*(x)*(x) + (dx)*(dx)))

/* Read curve data */
void read_curve(FILE *fp, int io, curve *c, char *label)
{
	int i; int n = 0, size = 64; double **m = matrix(3, size);
	size_t bsize = 128; char *buffer = (char *)xmalloc(bsize);
	
	/* Read curve data */
	while (getline(&buffer, &bsize, fp) != -1) {
		if (*buffer == '#') continue;
		if (*buffer == ';') break;
		
		if (n >= size) m = grow_matrix(m, 3, size<<=1);
		
		m[0][n] = m[1][n] = m[2][n] = 0.0;
		if (sscanf(buffer, "%lf %lf %lf", &(m[io?1:0][n]), &(m[io?0:1][n]), &(m[2][n])) < 2)
			error("syntax error in curve data");
		
		m[2][n] = LIMDEV(m[1][n], m[2][n]);
		
		n++;
	}
	free(buffer);
	
	
	/* Fit curve data */
	if (!n) error("No curve data supplied");
	
	if (verbose > 3) {
		for (i = 0; i < n; i++)
			fprintf(stderr, "\t\t%12.10g %12.10g \t(%12.10g)\n", m[0][i], m[1][i], m[2][i]);
	}
	
	if (verbose > 1) { fprintf(stderr, "\t%i points, fitting curve...", n); fflush(stderr); }
	
	c->n = n;
	c->data = m;
	c->fit = fit_curve(m, n);
	
	if (verbose > 1) {
		fprintf(stderr, " done.\n");
		
		if (verbose > 2)
			fprintf(stderr, "\t(Best fit: y = %9.7f * x^%9.7f + %9.7f; Chi2 = %.7f)\n",
					c->fit[0], c->fit[1], c->fit[2], c->fit[3]);
	}
	
	
	/* Plot curve fit and data using gnuplot */
	#ifdef DEBUG
	#define PTS 257
	#define DMIN 1.0e-3
	{
		double x, y; FILE *gp = popen("gnuplot", "w");
		
		fprintf(gp, "set terminal postscript eps enhanced color; set title '%s'; set key left\n", label);
		fprintf(gp, "set size 0.75,1.50; set lmargin 10; set multiplot; set size 0.75,1.03; set origin 0,0.47\n");
		fprintf(gp, "set ylabel 'Y (device)'; set xrange [%g:1]; set yrange [%g:1]; set logscale\n", DMIN, DMIN);
		fprintf(gp, "plot '-' title 'Y->X lookup' with lines lt 2, '-' title 'X->Y lookup' with lines lt 3, '-' title 'data' with yerrorbars lt 1 pt 6\n");
		
		/* forward lookup */
		for (i = 0; i < PTS; i++) {
			y = pow(DMIN, i/(PTS-1.0)); x = lu_curve(c->fit, y);
			fprintf(gp, "%12.10g\t%12.10g\n", x, y);
		}
		fprintf(gp, "e\n");
		
		/* backward lookup */
		for (i = 0; i < PTS; i++) {
			x = pow(DMIN, i/(PTS-1.0)); y = lu_curve_1(c->fit, x);
			fprintf(gp, "%12.10g\t%12.10g\n", x, y);
		}
		fprintf(gp, "e\n");
		
		/* data points and error bars */
		for (i = 0; i < n; i++) {
			fprintf(gp, "%12.10g\t%12.10g\t%12.10g\n", m[0][i], m[1][i], m[2][i]);
		}
		fprintf(gp, "e\n");
		
		fprintf(gp, "set size 0.75,0.50; set origin 0,0; set title; set nokey; set nologscale y\n");
		fprintf(gp, "set xlabel 'X (target)'; set ylabel 'fit residual'; set yrange [*:*]\n");
		fprintf(gp, "plot 1 lt 0, '-' title 'data' with yerrorbars lt 1 pt 6\n");
		
		/* fit residuals */
		for (i = 0; i < n; i++) {
			x = m[0][i]; y = lu_curve_1(c->fit, x);
			fprintf(gp, "%g %g %g\n", x, m[1][i]/y, m[2][i]/y);
		}
		fprintf(gp, "e\n");
		
		fprintf(gp, "set nomultiplot\n"); pclose(gp);
	}
	#undef DMIN
	#undef PTS
	#endif /* DEBUG */
}


/* Read LUT data */
void read_lut(FILE *fp)
{
	int i, j; char *p, *q; int n = 0, size = 64;
	datapt *data = (datapt *)xmalloc(size*sizeof(datapt));
	size_t bsize = 128; char *buffer = (char *)xmalloc(bsize);
	
	/* Read LUT data */
	while (getline(&buffer, &bsize, fp) != -1) {
		if (*buffer == '#') continue;
		if (*buffer == ';') break;
		if (ignore_lut) continue;
		
		/* initialize data structure */
		if (n >= size) data = (datapt *)xrealloc(data, (size <<= 1)*sizeof(datapt));
		
		data[n].flag = 0; data[n].label = NULL;
		
		for (i = 0; i < MAXCHANNELS; i++)
			data[n].in[i] = data[n].out[i] = data[n].var[i] = 0.0;
		
		/* read data */
		if (sscanf(buffer, "%as %n", &(data[n].label), &i) < 1)
			error("syntax error in lut data");
		p = buffer + i; q = NULL;
		
		for (i = 0; i < ins_channels; i++) {
			data[n].in[i] = strtod(p, &q);
			if (p == q) error("syntax error in lut data"); else p = q;
		}
		
		for (i = 0; i < outs_channels; i++) {
			data[n].out[i] = strtod(p, &q);
			if (p == q) error("syntax error in lut data"); else p = q;
		}
		
		for (i = 0; i < ins_channels; i++) {
			data[n].var[i] = strtod(p, &q); p = q;
		}
		
		/* flag points that match ignored pattern */
		if (ignore_pattern && !fnmatch(ignore_pattern, data[n].label, 0) && !ignore_none) data[n].flag = 1;
		
		/* pull output data back to XYZ profile connection space */
		outcurves_1(NULL, data[n].XYZ, data[n].out);
		(*outs2XYZ)(data[n].XYZ, data[n].XYZ);
		
		n++;
	}
	
	/* Tuck a copy of calibration data away for later */
	calibration_pts = n; calibration_data = data;
	
	/* Linear fit to calibration data */
	if (!ignore_lut && n) {
		int *f, idx[n]; double **m = matrix(9, n);
		
		for (i = 0, n = 0; i < calibration_pts; i++) if (!data[i].flag) {
			for (j = 0; j < 3; j++) {
				m[j  ][n] = data[i].XYZ[j];
				m[j+3][n] = data[i].in[j];
				m[j+6][n] = LIMDEV(data[i].in[j], data[i].var[j]);
			}
			
			idx[n++] = i;
		}
		
		/* Find best-fitting linear RGB transform, while adjusting curves */
		if (n) {
			double *C[] = { ins_curves[0].fit, ins_curves[1].fit, ins_curves[2].fit };
			
			f = fit_matrix(m, n, M_XYZ2RGB, C); inv33(M_XYZ2RGB, M_RGB2XYZ);
		}
		
		for (i = 0; i < n; i++) if (f[i] && !ignore_none) data[idx[i]].flag = f[i];
		
		/* Flag potentially clipped points */
		for (i = 0, n = 0; i < calibration_pts; i++) {
			/* push input data forward to linearized device space */
			incurves(NULL, data[i].DEV, data[i].in);
			
			if (verbose > 3) fprintf(stderr,
				"%c%18s: %12.10g %12.10g %12.10g -> %12.10g %12.10g %12.10g\n",
				(data[i].flag ? '#' : ' '), data[i].label,
				data[i].DEV[0], data[i].DEV[1], data[i].DEV[2],
				data[i].XYZ[0], data[i].XYZ[1], data[i].XYZ[2]);
			
			if (!data[i].flag) n++;
		}
		
		if (verbose > 1) { fprintf(stderr, "\t%i points (%i ignored), matrix fit...", calibration_pts, calibration_pts-n); fflush(stderr); }
		
		if (verbose > 2) {
			fprintf(stderr, "\n\tXYZ -> RGB conversion matrix:\n");
			fprintf(stderr, "\t\t[R]   [ %13.10g %13.10g %13.10g ] [X]\n", M_XYZ2RGB[0][0], M_XYZ2RGB[0][1], M_XYZ2RGB[0][2]);
			fprintf(stderr, "\t\t[G] = [ %13.10g %13.10g %13.10g ] [Y]\n", M_XYZ2RGB[1][0], M_XYZ2RGB[1][1], M_XYZ2RGB[1][2]);
			fprintf(stderr, "\t\t[B]   [ %13.10g %13.10g %13.10g ] [Z]\n", M_XYZ2RGB[2][0], M_XYZ2RGB[2][1], M_XYZ2RGB[2][2]);
			
			fprintf(stderr, "\tAdjusted curve fits:\n");
			for (j = 0; j < 3; j++) fprintf(stderr, "\t       y[%i] = %9.7f * x[%i]^%9.7f + %9.7f; Chi2 = %.7f\n",
					j, ins_curves[j].fit[0], j, ins_curves[j].fit[1], ins_curves[j].fit[2], ins_curves[j].fit[3]);
		}
		
		/* Plot linear lookup table fit scatter using gnuplot */
		#ifdef DEBUG
		{
			int k; double t[3]; FILE *gp = popen("gnuplot", "w");
			
			fprintf(gp, "set terminal postscript eps enhanced color; set nokey\n");
			fprintf(gp, "set size 0.75,1.50; set lmargin 10; set bmargin 1.5; set multiplot; set size 0.75,0.47\n");
			fprintf(gp, "set xrange [-0.1:1.1]; set yrange [0:1]\n");
			
			for (k = 0; k < 3; k++) {
				if (k == 2) {
					fprintf(gp, "set xlabel 'X (target pullback)'\n");
					fprintf(gp, "set label 'Scatter after linear regression' at screen 0.4,1.47 center\n");
				}
				fprintf(gp, "set ylabel 'Y (device, channel %i)'\n", k);
				fprintf(gp, "set origin 0,%f\n", 0.06 + 0.47*(2-k));
				fprintf(gp, "plot x with line lt 2, '-' title 'data' with yerrorbars lt 3 pt 6, '-' title 'data' with yerrorbars lt 1 pt 6\n");
				
				/* data points flagged as outliers */
				for (i = 0; i < calibration_pts; i++) if (data[i].flag) {
					(*XYZ2ins)(data[i].XYZ, t); incurves_1(NULL, t, t);
					fprintf(gp, "%12.10g\t%12.10g\t%12.10g\n", t[k], data[i].in[k], data[i].var[k]);
				}
				fprintf(gp, "e\n");
				
				/* data points and error bars */
				for (i = 0; i < calibration_pts; i++) if (!data[i].flag) {
					(*XYZ2ins)(data[i].XYZ, t); incurves_1(NULL, t, t);
					fprintf(gp, "%12.10g\t%12.10g\t%12.10g\n", t[k], data[i].in[k], data[i].var[k]);
				}
				fprintf(gp, "e\n");
			}
			
			fprintf(gp, "set nomultiplot\n"); pclose(gp);
		}
		#endif /* DEBUG */
		
		free_matrix(m); free_vector(f);
	}
	
	/* Polynomial regression and interpolation */
	if (!ignore_lut && !matrix_only && n > 64) {
		double **m = matrix(6, n), **im = matrix(6, n);
		
		for (i = 0, n = 0; i < calibration_pts; i++) if (!data[i].flag) {
			double t[3]; (*ins2XYZ)(data[i].DEV, t);
			
			for (j = 0; j < 3; j++) {
				m[j  ][n] = im[j+3][n] = t[j];
				m[j+3][n] = im[j  ][n] = data[i].XYZ[j];
			}
			
			n++;
		}
		
		if (use_poly = 1) { /* polynomial fit */
			if (verbose > 1) { fprintf(stderr, " poly fit..."); fflush(stderr); }
			
			P = fit_poly(m, n); if (invertible) P1 = fit_poly(im, n);
			
			for (i = 0; i < n; i++) {
				double t1[3] = {  m[0][i],  m[1][i],  m[2][i] };
				double t2[3] = { im[0][i], im[1][i], im[2][i] };
				
				evalf_poly(P, t1, t1); for (j = 0; j < 3; j++) m[j][i] = t1[j];
				if (invertible) { evalf_poly(P1, t2, t2); for (j = 0; j < 3; j++) im[j][i] = t2[j]; }
			}
		}
		
		if (use_interp = 1) { /* prepare for interpolation */
			if (verbose > 1) { fprintf(stderr, " interpolation..."); fflush(stderr); }
			
			for (i = 0; i < n; i++) {
				double t1[3] = {  m[0][i],  m[1][i],  m[2][i] };
				double t2[3] = { im[0][i], im[1][i], im[2][i] };
				
				XYZ2Lab(t1, t1); XYZ2Lab(t2, t2);
				
				for (j = 0; j < 3; j++) {
					m[j+3][i] -= m[j][i]; m[j][i] = t1[j];
					im[j+3][i] -= im[j][i]; im[j][i] = t2[j];
				}
			}
			
			prepint3d(T = m, PTS = n, Q); if (invertible) prepint3d(T1 = im, PTS = n, Q1);
		}
	}
	
	if (verbose > 1) fprintf(stderr, " %s.\n", (ignore_lut ? "ignored" : "done"));
	
	free(buffer);
}


/* Read calibration data */
void read_calibration_data(FILE *fp)
{
	int i;
	size_t bsize = 128;
	char *p, *buffer = (char *)xmalloc(bsize);
	
	if (verbose) fprintf(stderr, "Reading calibration data...\n");
	
	while (getline(&buffer, &bsize, fp) != -1) {
		buffer[strlen(buffer)-1] = 0;
		
		if (*buffer == '#') continue;
		
		if ((p = strstr(buffer, "WHITEPOINT:"))) {
			if (verbose > 1) fprintf(stderr, "\t%s\n", buffer);
			
			if (sscanf(p+11, "%lf %lf %lf", media_white_pt, media_white_pt+1, media_white_pt+2) != 3)
				error("syntax error in whitepoint data");
			
		} else if ((p = strstr(buffer, "BLACKPOINT:"))) {
			if (verbose > 1) fprintf(stderr, "\t%s\n", buffer);
			
			if (sscanf(p+11, "%lf %lf %lf", media_black_pt, media_black_pt+1, media_black_pt+2) != 3)
				error("syntax error in blackpoint data");
			
		} else if (sscanf(buffer, "CURVE IN%i", &i) == 1) {
			if (verbose > 1) {
				fprintf(stderr, "\t%s", buffer);
				if (verbose > 2) fprintf(stderr, "\n");
				else fflush(stderr);
			}
			
			read_curve(fp, 1 /*in*/, &(ins_curves[i]), buffer);
			use_ins_curves |= 1 << i;
			
		} else if (sscanf(buffer, "CURVE OUT%i", &i) == 1) {
			if (verbose > 1) {
				fprintf(stderr, "\t%s", buffer); 
				if (verbose > 2) fprintf(stderr, "\n");
				else fflush(stderr);
			}
			
			read_curve(fp, 0 /*out*/, &(outs_curves[i]), buffer);
			use_outs_curves |= 1 << i;
			
		} else if (strstr(buffer, "LUT:")) {
			if ((use_ins_curves && (use_ins_curves != (1 << ins_channels) - 1)) ||
			    (use_outs_curves && (use_outs_curves != (1 << outs_channels) - 1)) )
				error("All curves must be specified before trying to read in LUT");
			
			if (verbose > 1) {
				fprintf(stderr, "\t%s", buffer);
				if (verbose > 2 && !ignore_lut) fprintf(stderr, "\n");
				else fflush(stderr);
			}
			
			read_lut(fp);
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
		h->cmmId = h->creator = str2tag("SCRS");
		
		/* Values that must be set before writing */
		h->deviceClass     = class;
    		h->colorSpace      = ins;
    		h->pcs             = outs;
		
		/* I am somewhat confused on the relative vs. absolute profile issue,
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
		
		if (!description) asprintf(&description, "SCARSE: %s", profile);
		
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
		
		r->data[0].X = M_RGB2XYZ[0][0];
		r->data[0].Y = M_RGB2XYZ[1][0];
		r->data[0].Z = M_RGB2XYZ[2][0];
		g->data[0].X = M_RGB2XYZ[0][1];
		g->data[0].Y = M_RGB2XYZ[1][1];
		g->data[0].Z = M_RGB2XYZ[2][1];
		b->data[0].X = M_RGB2XYZ[0][2];
		b->data[0].Y = M_RGB2XYZ[1][2];
		b->data[0].Z = M_RGB2XYZ[2][2];
	}
	
	/* Red, Green and Blue Tone Reproduction Curve Tags: */
	{
		int i; double p[MAXCHANNELS];
		
		icmCurve *r = (icmCurve *)icco->add_tag(
				icco, icSigRedTRCTag, icSigCurveType);
		icmCurve *g = (icmCurve *)icco->add_tag(
				icco, icSigGreenTRCTag, icSigCurveType);
		icmCurve *b = (icmCurve *)icco->add_tag(
				icco, icSigBlueTRCTag, icSigCurveType);
		
		if (!r || !g || !b)
			error("add_tag failed: %d, %s", icco->errc, icco->err);
		
		if (use_ins_curves) { /* tabulated curve */
			r->flag = g->flag = b->flag = icmCurveSpec;
			r->size = g->size = b->size = lut_size_1d;
			
			/* Allocate space */
			r->allocate((icmBase *)r);
			g->allocate((icmBase *)g);
			b->allocate((icmBase *)b);
			
			for (i = 0; i < lut_size_1d; i++) {
				p[0] = p[1] = p[2] = i/(lut_size_1d-1.0);
				
				incurves(NULL, p, p);	/* Transfer function */
				
				r->data[i] = p[0];
				g->data[i] = p[1];
				b->data[i] = p[2];
			}
		} else if (ins_gamma != 1.0) { /* power law */
			r->flag = g->flag = b->flag = icmCurveGamma;
			
			/* Allocate space */
			r->allocate((icmBase *)r);
			g->allocate((icmBase *)g);
			b->allocate((icmBase *)b);
			
			r->data[0] = g->data[0] = b->data[0] = ins_gamma;
		} else { /* identity */
			r->flag = g->flag = b->flag = icmCurveLin;
			
			/* Allocate space */
			r->allocate((icmBase *)r);
			g->allocate((icmBase *)g);
			b->allocate((icmBase *)b);
		}
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
	if (verbose > 3) {
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
	if (calibration_pts) {
		int i = 0, n = 0; double ip[3], op[3], tp[MAXCHANNELS];
		double e = 0.0, dE[3] = {-HUGE /* max */, HUGE /* min */, 0.0 /* avg */};
		
		/* Get a conversion object */
		icmLuBase *luo = icco->get_luobj(icco, icmFwd, icmDefaultIntent, icmLuOrdNorm);
		if (!luo) error("%s: Error %d, %s", file, icco->errc, icco->err);
		
		for (i = 0; i < calibration_pts; i++) {
			if (luo->lookup(luo, tp, calibration_data[i].in) > 1)
				error("%s: Error %d, %s", file, icco->errc, icco->err);
			
			(*outs2XYZ)(tp, ip); (*outs2XYZ)(calibration_data[i].out, op);
			
			e = dE_XYZ(ip, op);
			
			if (!calibration_data[i].flag) {
				if (e > dE[0]) dE[0] = e;
				if (e < dE[1]) dE[1] = e;
				dE[2] += e; n++;
			}
			
			if (verbose > 4) fprintf(stderr, "\n%20s:\tdE = %9.6f %s",
				calibration_data[i].label, e, calibration_data[i].flag ? "(ignored)" : "");
		}
		
		if (verbose > 1) fprintf(stderr, "\n\t Forward lookup dE = (%6.4g max, %6.4g min, %6.4g avg)", dE[0], dE[1], dE[2]/n);
		
		luo->free(luo);
	}
	
	/* Test backward lookup */
	if (calibration_pts && invertible) {
		int i = 0, n = 0; double ip[3], op[3], tp[MAXCHANNELS];
		double e = 0.0, dE[3] = {-HUGE /* max */, HUGE /* min */, 0.0 /* avg */};
		
		/* Get a conversion object */
		icmLuBase *luo = icco->get_luobj(icco, icmBwd, icmDefaultIntent, icmLuOrdNorm);
		if (!luo) error("%s: Error %d, %s", file, icco->errc, icco->err);
		
		for (i = 0; i < calibration_pts; i++) {
			if (luo->lookup(luo, tp, calibration_data[i].out) > 1)
				error("%s: Error %d, %s", file, icco->errc, icco->err);
			
			(*ins2XYZ)(tp, ip); (*ins2XYZ)(calibration_data[i].in, op);
			
			e = dE_XYZ(ip, op);
			
			if (!calibration_data[i].flag) {
				if (e > dE[0]) dE[0] = e;
				if (e < dE[1]) dE[1] = e;
				dE[2] += e; n++;
			}
			
			if (verbose > 4) fprintf(stderr, "\n%20s:\tdE = %9.6f %s",
				calibration_data[i].label, e, calibration_data[i].flag ? "(ignored)" : "");
		}
		
		if (verbose > 1) fprintf(stderr, "\n\tBackward lookup dE = (%6.4g max, %6.4g min, %6.4g avg)\n", dE[0], dE[1], dE[2]/n);
		
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
	LookupPrimaries("Adobe", primaries, NULL);
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv-:c:i:o:p:C:UE:MI:m:d:r:b:s:")) != -1)
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
			if (!LookupPrimaries(optarg, primaries, &ins_gamma)) {
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
			calibration_data_stream = xfetch(NULL, optarg, "r");
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
		case 'I':				/* Ignore LUT points */
			if (!strcmp(optarg, "all"))  { ignore_lut = 1;  break; }
			if (!strcmp(optarg, "none")) { ignore_none = 1; break; }
			ignore_pattern = optarg;
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
				lut_size_1d = atoi(optarg);
			
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
	vcopy3(XYZ_WPT, media_white_pt);
	
	ins_channels = channels(ins);
	outs_channels = channels(outs);
	
	ins2XYZ = toXYZ(ins); XYZ2ins = fromXYZ(ins);
	outs2XYZ = toXYZ(outs); XYZ2outs = fromXYZ(outs);
	
	
	/* Read in calibration data if available */
	if (calibration_data_stream)
		read_calibration_data(calibration_data_stream);
	
	
	/* Build profile */
	build_profile(profile);
	
	/* Test profile accuracy */
	test_profile(profile);
	
	
	return 0;
}
