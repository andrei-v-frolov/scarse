/* $Id: cmap.c,v 1.3 2005/09/28 01:18:33 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Translate image through ICC profiles.
 * 
 * Copyright (C) 1999, 2000 Scarse Project
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

#define SELF "cmap"

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <icc.h>
#include "scarse.h"



/******************* Options and defaults *****************************/

/* Usage */
char *program_name = SELF;
char *usage_msg[] = {
	"Translate image through ICC profiles, Version " VERSION,
	"Author: Andrei Frolov <andrei@phys.ualberta.ca>",
	"",
	"Usage: " SELF " [options] [[-p|-r profile] ...] input.tif output.tif",
	"",
	" General options:",
	"  -h		print this message",
	"  -v		verbosity (cumulative)",
	"  -P		preview (scale image down to below 640x480 pixels)",
	"  --longopt	expand option macro 'longopt' as defined in etc/" SELF ".options",
	"",
	" ICC profile options:",
	"  -f function	f = forward, b = backwards, g = gamut, p = preview",
	"  -i intent	p = perceptual, r = relative colorometric,",
	"		s = saturation, a = absolute, x = absolute XYZ",
	"  -o order	n = normal (priority: lut > matrix > monochrome)",
	"		r = reverse (priority: monochrome > matrix > lut)",
	"  -p profile	ICC profile to apply",
	"  -r profile	reverse profile; equivalent to -fb -p profile",
	"  -e		embed last profile into the output image file",
	"",
	" Color adjustment options:",
	"  -O order	n = normal (color mapping order: profiles > adjustments)",
	"		r = reverse (color mapping order: adjustments > profiles)",
	"  -[B|W|G] [=]value[,value,value]",
	"		adjust (or set if =) [B]lack/[W]hite point or contrast [G]rade",
	"		if triplet of values is given, adjust each channel separately",
	"  -A flags	enable automatic adjustment; flag values can be:",
	"		w = white point (scalar), W = white point (vector)",
	"		b = black point (scalar), B = black point (vector)",
	"		g = contrast grade (scalar), G = contrast grade (vector)",
	"		L = auto levels (equivalent to WB)",
	"  -S split,shadow_exposure,highlight_exposure",
	"		split-contrast print; low and high contrast exposures are in %",
	"  -C length	apply center filter for lens of given focal length (35mm equiv)",
	"",
	" Data storage options: (partial support by formats other than TIFF)",
	"  -b bits	bits per sample (bits = 8 or 16)",
	"  -n		add sub-bit noise when quantizing",
	"  -w weave	c = contigious (e.g. RGBRGBRGB...),",
	"		i = interlaced (e.g. RRR...GGG...BBB...)",
	"  -c compress	l = lzw, z = zip, p = pack bits, n = none",
	"		l:2, z:2 = lzw, zip with horizontal differencing",
	NULL
};


/* ICC profile info */
typedef struct _iccProfile {
	const char *name;			/* Profile name */
	icc *icco;				/* ICC profile object */
	icmLuBase *luo;				/* Lookup object */
	icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
	int inn, outn;				/* Number of components */
	icmLuAlgType alg;			/* Type of lookup algorithm */
	
	struct _iccProfile *next;		/* Next profile in list */
} iccProfile;

/* Color adjustment flags */
enum {
	ADJ_WHITE_SCALAR	= 0x01,
	ADJ_WHITE_VECTOR	= 0x02,
	ADJ_BLACK_SCALAR	= 0x04,
	ADJ_BLACK_VECTOR	= 0x08,
	ADJ_CONTRAST_SCALAR	= 0x10,
	ADJ_CONTRAST_VECTOR	= 0x20,
	ADJ_CONTRAST_SPLIT	= 0x40,
	ADJ_CENTER		= 0x80,
	/* Composite flags */
	ADJ_WHITE		= ADJ_WHITE_SCALAR | ADJ_WHITE_VECTOR,
	ADJ_BLACK		= ADJ_BLACK_SCALAR | ADJ_BLACK_VECTOR,
	ADJ_LEVELS		= ADJ_WHITE | ADJ_BLACK,
	ADJ_CONTRAST		= ADJ_CONTRAST_SCALAR | ADJ_CONTRAST_VECTOR,
	ADJ_CURVES		= ADJ_CONTRAST | ADJ_CONTRAST_SPLIT,
	ADJ_LOCAL		= ADJ_LEVELS | ADJ_CURVES
};


/* Options */
static int              verbose = 0;
static int              preview = 0;
static unsigned long       skip = 1;

static icmLookupFunc       func = icmFwd;
static icRenderingIntent intent = icmDefaultIntent;
static icmLookupOrder     order = icmLuOrdNorm;
static iccProfile     *icc_pipe = NULL;
static int                embed = 0;

static int               adjust = 0;
static int          auto_adjust = 0;
static int            preadjust = 0;
static double               wpt = 0.0, wptv[3] = {0.0, 0.0, 0.0};
static double               bpt = 0.0, bptv[3] = {0.0, 0.0, 0.0};
static double             grade = 0.0, gradev[3] = {0.0, 0.0, 0.0};
static double             split = 0.0, shadow = 0.0, highlight = 0.0;
static double          cflength = 0.0;

static short      bitspersample = -1;
extern int              noisify;
static short       planarconfig = -1;
static short        compression = -1;
static short          predictor = 1;

static image                *in = NULL;
static image               *out = NULL;



/******************* Profile pipelines ********************************/

/* Add profile to conversion pipeline, open it if necessary */
static iccProfile *add_profile_to_pipe(const char *name, icc *icco, iccProfile *pipe, int tail)
{
	iccProfile *p = (iccProfile *)xmalloc(sizeof(iccProfile)), *t = pipe;
	
	p->name = name;
	
	/* Open ICC profile */
	if (!icco) {
		int rv = 0;
		FILE *file = xfopen(name, "rb");
		
		p->icco = new_icc();
		if (!p->icco) error("Creation of ICC object failed");
		
		/* Read profile */
		rv = p->icco->read(p->icco, file, 0);
		if (rv) error("%s: Error %d, %s", p->name, rv, p->icco->err);
	} else p->icco = icco;
	
	/* Get a conversion object */
	p->luo = p->icco->get_luobj(p->icco, func, intent, order);
	if (!p->luo) error("%s: Error %d, %s", p->name, p->icco->errc, p->icco->err);
	
	/* Get details of conversion (Arguments may be NULL if info not needed) */
	p->luo->spaces(p->luo, &(p->ins), &(p->inn), &(p->outs), &(p->outn), &(p->alg));
	
	p->next = NULL;
	
	if (!t) return p;
	
	/* Add new profile to the pipeline */
	if (tail) {	/* new profile goes last */
		while (t->next) t = t->next; t->next = p;
		
		if ((p->ins != t->outs) || (p->inn != t->outn))
			error("Input space [%s] of profile %s does not match output space [%s] of previous profile %s",
				ColorSpaceSignature2str(p->ins), p->name,
				ColorSpaceSignature2str(t->outs), t->name);
	} else {	/* new profile goes first */
		p->next = t;
		
		if ((p->outs != t->ins) || (p->outn != t->inn))
			error("Input space [%s] of profile %s does not match output space [%s] of previous profile %s",
				ColorSpaceSignature2str(t->outs), t->name,
				ColorSpaceSignature2str(p->ins), p->name);
	}
	
	return tail ? pipe : p;
}

/* Split profile pipe at specified PCS */
static void pcs_split(iccProfile *pipe, icColorSpaceSignature pcs, iccProfile **head, iccProfile **tail)
{
	iccProfile *p = pipe; *head = *tail = NULL;
	
	while (p) {
		if (p->outs == pcs) {
			*head = pipe;
			*tail = p->next;
			p->next = NULL;
			break;
		}
		
		if (!p->next) error("Cannot find %s profile connection space in the pipeline", ColorSpaceSignature2str(pcs));
		
		p = p->next;
	}
}

/* Join two profile pipes */
static iccProfile *join_pipes(iccProfile *head, iccProfile *tail)
{
	iccProfile *p = head;
	
	if (!head) return tail;
	if (!tail) return head;
	
	while (p->next) p = p->next; p->next = tail;
	
	if ((p->outs != tail->ins) || (p->outn != tail->inn))
			error("Input space [%s] of profile %s does not match output space [%s] of previous profile %s",
				ColorSpaceSignature2str(tail->ins), tail->name,
				ColorSpaceSignature2str(p->outs), p->name);
	
	return head;
}

/* Get input color space */
static icColorSpaceSignature input_space(iccProfile *pipe)
{
	if (pipe) return pipe->ins; else return -1;
}

/* Get output color space */
static icColorSpaceSignature output_space(iccProfile *pipe)
{
	iccProfile *p = pipe;
	
	while (p) { if (!p->next) return p->outs; p = p->next; } return -1;
}

/* Dump conversion pipeline info */
static void dump_pipe(iccProfile *pipe)
{
	iccProfile *p = pipe;
	
	if (!pipe) return;
	
	if (adjust && preadjust)
		fprintf(stderr, "(Adjustment) -> ");
	
	while (p) {
		fprintf(stderr, "[%s] -> %s (%s) -> ",
			ColorSpaceSignature2str(p->ins),
			p->name, LuAlg2str(p->alg));
		if (!p->next) fprintf(stderr, "[%s]",
			ColorSpaceSignature2str(p->outs));
		
		p = p->next;
	}
	
	if (adjust && !preadjust)
		fprintf(stderr, " -> (Adjustment)");
	
	fprintf(stderr, "\n");
	
	if (verbose > 2) {
		fprintf(stderr, "\n");
		
		p = pipe;
		while (p) {
			fprintf(stderr, "%s ", p->name);
			p->icco->header->dump(p->icco->header, stderr, 1);
			
			p = p->next;
		}
	}
}

/* Run pixel through ICC profile pipeline */
static int lu_pipe(iccProfile *pipe, double in[], double out[])
{
	iccProfile *p = pipe;
	int j = 0, rv = 0;
	double tmp[MAXCHANNELS], *alt[] = {tmp, out};
	
	if (!pipe) return 0;
	
	while (p) {
		rv |= p->luo->lookup(p->luo, p->next ? alt[j & 1]: out, j ? alt[(j-1) & 1] : in); j++;
		
		p = p->next;
	}
	
	if (rv & (~0x01)) error("%s: Error %d, %s",
		p->name, p->icco->errc, p->icco->err);
	
	return rv;
}

/* Free profiles in a conversion pipeline */
static iccProfile *free_profiles(iccProfile *pipe)
{
	iccProfile *p = pipe, *n;
	
	while (p) {
		fclose(p->icco->fp);
		p->luo->free(p->luo);
		p->icco->free(p->icco);
		n = p->next;
		free(p);
		p = n;
	}
	
	return NULL;
}



/******************* Embedded profiles ********************************/

#ifdef TIFFTAG_ICCPROFILE

/* Read embedded ICC profile from image */
static icc *read_embedded_profile(image *img)
{
	FILE *tmp;
	icc *icco;
	int rv = 0;
	void *buffer;
	unsigned long l;
	
	/* Extract embedded profile to temporary file */
	if (!(*img->get_field)(img, TIFFTAG_ICCPROFILE, &l, &buffer)) return NULL;
	if (!(tmp = tmpfile()) || fwrite(buffer, 1, l, tmp) != l)
		error("%s: Error extracting embedded profile", img->file);
	
	/* Read profile */
	if (!(icco = new_icc()))
		error("Creation of ICC object failed");
	if ((rv = icco->read(icco, tmp, 0)))
		error("%s: Error %d, %s", img->file, rv, icco->err);
	
	return icco;
}

/* Embed ICC profile into image */
static void write_embedded_profile(image *img, icc *icco)
{
	unsigned long l = icco->header->size;
	void *buffer = xmalloc((size_t)(l));
	
	/* Read profile data into memory */
	if (fseek(icco->fp, icco->of, 0) || fread(buffer, 1, l, icco->fp) != l)
		error("Error re-reading profile to embed");
	
	/* Embed profile */
	if (!(*img->set_field)(img, TIFFTAG_ICCPROFILE, l, buffer))
		warning("%s: Profile embedding is not supported by image format", img->file);
	
	free(buffer);
}

#else  /* Embedded profiles are not supported by libtiff prior to v3.5 */
static icc *read_embedded_profile(image *img)
	{ warning("Embedded profile support is not compiled in"); return NULL; }
static void write_embedded_profile(image *img, icc *icco)
	{ warning("Embedded profile support is not compiled in"); }

#endif /* TIFFTAG_ICCPROFILE */



/****************** Automatic color adjustment ************************/

#define MOMENTA 5	/* number of distribution momenta to collect */

/* Pixel value distribution and statistics */
typedef struct {
	double min, max, avg[MOMENTA];			/* statistics */
	unsigned long under, count[256], over, total;	/* histogram */
} distribution;

/* Recommended adjustments */
typedef struct {
	double wpt, bpt, grade;
} recommendation;


/* Initialize distribution structure */
static distribution *new_dist()
{
	int i;
	distribution *d = (distribution *)xmalloc(sizeof(distribution));
	
	d->min = +HUGE; d->max = -HUGE;
	for (i = 0; i < MOMENTA; i++) d->avg[i] = 0.0;
	
	d->under = d->over = d->total = 0;
	for (i = 0; i < 256; i++) d->count[i] = 0;
	
	return d;
}

/* Collect distribution data */
static void collect_dist(distribution *d, double x)
{
	int i;
	double X;
	
	if (x < d->min) d->min = x;
	if (x > d->max) d->max = x;
	
	for (i = 1, X = x; i < MOMENTA; i++, X *= x) d->avg[i] += X;
	
	if (x < 0.0) d->under++;
	else if (x > 1.0) d->over++;
	else d->count[(int)floor(255.0*x+0.5)]++;
	
	d->total++;
}

/* Analyze distribution data and recommend adjustments */
static void analyze_dist(distribution *d, const char *channel, recommendation *r)
{
	int i;
	double x;
	unsigned long c;
	
	if (!d->total) error("Internal error: no data to analyze in analyze_dist()");
	
	/* be more aggressive setting white level - necessary to handle dust on negs */
	i = 256; c = d->over; while (c < (d->total >> 12)) c+= d->count[--i];
	if (i < 255) d->max = (double)(i+1)/255.0;
	
	/* median */
	i = 0; c = d->under + d->count[0]; while (c < (d->total >> 1)) c += d->count[++i];
	x = ((double)(d->total)/2.0 - (double)(c-d->count[i]))/((double)(d->count[i]));
	d->avg[0] = (x + (double)(i-1))/255.0;
	
	/* mean and higher momenta */
	for (i = 1; i < MOMENTA; i++) d->avg[i] /= (double)(d->total);
	
	/* dump stats */
	if (verbose > 4 && channel) {
		fprintf(stderr, "\n%s channel data distribution:\n", channel);
		fprintf(stderr, "  %c_min = %10.6g\n", *channel, d->min);
		fprintf(stderr, "  %c_max = %10.6g\n", *channel, d->max);
		
		for (i = 0; i < MOMENTA; i++)
			fprintf(stderr, "  <%c^%i> = %10.6g\n", *channel, i, d->avg[i]);
		
		#if (MOMENTA > 2)
		{
			double var = d->avg[2] - d->avg[1]*d->avg[1];
			fprintf(stderr, "    Var = %10.6g\n", var);
		
		#if (MOMENTA > 3)
		{
			double skew = (d->avg[3] - 3.0*d->avg[2]*d->avg[1] + 2.0*pow(d->avg[1], 3.0))/pow(var, 1.5);
			fprintf(stderr, "   Skew = %10.6g\n", skew);
		}
		#endif
		
		#if (MOMENTA > 4)
		{
			double kurt = (d->avg[4] - 4.0*d->avg[3]*d->avg[1] + 6.0*d->avg[2]*pow(d->avg[1], 2.0) - 3.0*pow(d->avg[1], 4.0))/pow(var, 2.0);
			fprintf(stderr, "   Kurt = %10.6g\n", kurt);
		}
		#endif
		
		}
		#endif /* (MOMENTA > 2) */
	}
	
	/* make adjustment recommendation */
	{
		double min, max, a, b, g;
		
		min = (auto_adjust & ADJ_BLACK) ? d->min : 0.0;
		max = (auto_adjust & ADJ_WHITE) ? d->max : 1.0;
		
		a = 1.0/(max - min); b = - a*min;
		g = log(a * d->avg[0] + b)/log(0.18);
		
		r->wpt = 100.0 * log10(a);
		r->bpt = -100.0 * log10(1.0 - b);
		r->grade = -10.0 * log10(g);
		
		if (verbose > 4 && channel)
			printf("(recommend adjustments: W = %-.2g, B = %-.2g, C = %-.2g)\n", r->wpt, r->bpt, r->grade);
	}
}


/* Analyze image and make automatic adjustments */
static void analyze_image(image *in)
{
	int i;
	double p[MAXCHANNELS], L;
	unsigned long x, y, X, Y;
	iccProfile *icc_pipe1, *icc_pipe2;
	distribution *d[4]; recommendation r[4];
	char *channel[4] = {"Red", "Green", "Blue", "Y"};
	unsigned long pxls = in->w * in->h, s = pxls ? (pxls-1)/640000 + 1 : 1;
	
	/* initialize */
	for (i = 0; i < 4; i++) d[i] = new_dist();
	
	/* scan image */
	pcs_split(icc_pipe, icSigXYZData, &icc_pipe1, &icc_pipe2);
	
	for (y = 0, Y = 0; Y < in->h; y++, Y += s) {
		IMG_ReadRow(in, Y);
		
		for (x = 0, X = 0; X < in->w; x++, X += s) {
			IMG_ReadPixel(in, X, p);
			
			if (preadjust) for (i = 0; i < 3; i++) collect_dist(d[i], p[i]);
			
			lu_pipe(icc_pipe1, p, p);
			
			if (icc_pipe) L = p[1];
			else L = 0.2973*pow(p[0], DEFAULT_GAMMA) + 0.6274*pow(p[1], DEFAULT_GAMMA) + 0.0753*pow(p[2], DEFAULT_GAMMA);
			collect_dist(d[3], L);
			
			lu_pipe(icc_pipe2, p, p);
			
			if (!preadjust) for (i = 0; i < 3; i++) collect_dist(d[i], p[i]);
		}
	}
	
	icc_pipe = join_pipes(icc_pipe1, icc_pipe2);
	
	/* do analysis */
	for (i = 0; i < 4; i++)
		analyze_dist(d[i], channel[i], &(r[i]));
	
	/* white point */
	if (auto_adjust & ADJ_WHITE) {
		if (auto_adjust & ADJ_WHITE_VECTOR) {
			for (i = 0; i < 3; i++)
				wptv[i] += r[i].wpt;
		} else {
			double w = HUGE; /* scalar adjustment */
			
			for (i = 0; i < 3; i++)
				if (r[i].wpt < w) w = r[i].wpt;
			
			wpt += w;
		}
	}
	
	/* black point */
	if (auto_adjust & ADJ_BLACK) {
		if (auto_adjust & ADJ_BLACK_VECTOR) {
			for (i = 0; i < 3; i++)
				bptv[i] += r[i].bpt;
		} else {
			double b = -HUGE; /* scalar adjustment */
			
			for (i = 0; i < 3; i++)
				if (r[i].bpt > b) b = r[i].bpt;
			
			bpt += b;
		}
	}
	
	/* contrast grade */
	if (auto_adjust & ADJ_CONTRAST) {
		if (auto_adjust & ADJ_CONTRAST_VECTOR) {
			double g = (r[0].grade + r[1].grade + r[2].grade)/3.0;
			
			for (i = 0; i < 3; i++)
				gradev[i] += r[i].grade - g;
		}
		
		grade += r[3].grade;
	}
}


/* Print out color adjustment settings */
static void dump_adjustment()
{
	if (adjust) fprintf(stderr, "\n");
	
	if (adjust & ADJ_LEVELS) {
		fprintf(stderr, "Adjusting levels:\n");
		if (adjust & ADJ_WHITE)
			fprintf(stderr, "\t\t* white [%+4.2g, %+4.2g, %+4.2g]\n",
				wpt+wptv[0], wpt+wptv[1], wpt+wptv[2]);
		if (adjust & ADJ_BLACK)
			fprintf(stderr, "\t\t* black [%+4.2g, %+4.2g, %+4.2g]\n",
				bpt+bptv[0], bpt+bptv[1], bpt+bptv[2]);
	}
	
	if (adjust & ADJ_CURVES) {
		fprintf(stderr, "Adjusting curves:\n");
		if (adjust & ADJ_CONTRAST)
			fprintf(stderr, "\t\t* contrast grade [%+4.2g, %+4.2g, %+4.2g]\n",
				grade+gradev[0], grade+gradev[1], grade+gradev[2]);
		if (adjust & ADJ_CONTRAST_SPLIT)
			fprintf(stderr, "\t\t* split-contrast [%+4.2g, %4.2g, %4.2g]\n",
				split, shadow, highlight);
	}
	
	if (adjust & ADJ_CENTER)
		fprintf(stderr, "Applying center filter for %gmm lens (35mm equivalent)\n", cflength);
}



/****************** Color correction loops ****************************/

/* Copy colors identically */
static void copy_loop(image *in, image *out)
{
	double p[MAXCHANNELS];
	unsigned long x, y, X, Y;
	
	for (y = 0, Y = 0; y < out->h; y++, Y += skip) {
		IMG_ReadRow(in, Y);
		
		for (x = 0, X = 0; x < out->w; x++, X += skip) {
			IMG_ReadPixel(in, X, p);
			IMG_WritePixel(out, x, p);
		}
		
		IMG_WriteRow(out, y);
	}
}


/* Translate colors through profile pipeline */
static void calibration_loop(image *in, image *out)
{
	double p[MAXCHANNELS];
	unsigned long x, y, X, Y, clipped = 0;
	
	for (y = 0, Y = 0; y < out->h; y++, Y += skip) {
		IMG_ReadRow(in, Y);
		
		for (x = 0, X = 0; x < out->w; x++, X += skip) {
			IMG_ReadPixel(in, X, p);
			clipped += lu_pipe(icc_pipe, p, p) |
			IMG_WritePixel(out, x, p);
		}
		
		IMG_WriteRow(out, y);
	}
	
	if (verbose && clipped)
		warning("%s: %lu pixels were out of gamut and had been clipped (%4.2g%%)",
			out->file, clipped, 100.0*clipped/out->w/out->h);
}


/* Linear table interpolation accelerator */
static double ltinterp(double y[256], double x, int *clip)
{
	double sx = 255.0*x, ix = floor(sx); int i = (int)ix;
	
	if (i < 0) { *clip |= 1; return y[0]; }
	else if (i > 254) { *clip |= 1; return y[255]; }
	
	return y[i] + (sx-ix)*(y[i+1]-y[i]);
}

/* Translate colors through profile pipeline and adjust them */
static void adjustment_loop(image *in, image *out)
{
	int i, k, rv;
	double p[MAXCHANNELS];
	iccProfile *icc_pipe1, *icc_pipe2;
	double a[3], b[3], **c, f2, norm;
	unsigned long x, y, X, Y, clipped = 0;
	
	/* Levels */
	if (adjust & ADJ_LEVELS) {
		for (i = 0; i < 3; i++) {
			a[i] = pow(10.0, (wpt+wptv[i])/100.0);
			b[i] = 1.0 - pow(10.0, -(bpt+bptv[i])/100.0);
		}
	}
	
	/* Curves */
	if (adjust & ADJ_CURVES) {
		double g, s, v;
		
		c = matrix(3, 256);
		s = pow(10.0, split/10.0);
		
		for (i = 0; i < 3; i++) {
			g = pow(10.0, (grade+gradev[i])/10.0) - 1.0;
			
			for (k = 0; k < 256; k++) {
				v = ((double)(k))/255.0;
				c[i][k] = pow(v, g) * (shadow*pow(v, 1.0/s) + (100.0-shadow-highlight)*v + highlight*pow(v, s))/100.0;
			}
		}
	}
	
	/* Center filter */
	if (adjust & ADJ_CENTER) {
		int p = 0;
		double dx, dy, c2, c4;
		
		f2 = (out->h*out->h + out->w*out->w) * cflength*cflength/468.0;
		
		/* Calculate normalization factor */
		norm = 0.0;
		for (y = 0; y < out->h; y++)
			for (x = 0; x < out->w; x++) {
				dx = (double)(x << 1) - (double)(out->w);
				dy = (double)(y << 1) - (double)(out->h);
				c2 = 1.0 + (dx*dx + dy*dy)/f2; c4 = c2*c2;
				
				norm += c4; p++;
		}
		norm = (p > 1024) ? norm/p : 1.0;
	}
	
	
	/* Correction loop */
	pcs_split(icc_pipe, icSigXYZData, &icc_pipe1, &icc_pipe2);
	
	for (y = 0, Y = 0; y < out->h; y++, Y += skip) {
		IMG_ReadRow(in, Y);
		
		for (x = 0, X = 0; x < out->w; x++, X += skip) {
			IMG_ReadPixel(in, X, p); rv = 0;
			
			for (k = 0; k < 2; k++) if ((k+preadjust) & 0x01) {
				if (adjust & ADJ_LEVELS) {
					for (i = 0; i < 3; i++)
						p[i] = a[i]*p[i] + b[i];
				}
				
				if (adjust & ADJ_CURVES) {
					for (i = 0; i < 3; i++)
						p[i] = ltinterp(c[i], p[i], &rv);
				}
			} else {
				rv |= lu_pipe(icc_pipe1, p, p);
				
				if (adjust & ADJ_CENTER) {
					double dx = (double)(x << 1) - (double)(out->w);
					double dy = (double)(y << 1) - (double)(out->h);
					double c2 = 1.0 + (dx*dx + dy*dy)/f2, c4 = c2*c2/norm;
					
					/* We assume default gamma if not in XYZ */
					if (!icc_pipe) c4 = pow(c4, 1.0/DEFAULT_GAMMA);
					
					for (i = 0; i < 3; i++) p[i] *= c4;
				}
				
				rv |= lu_pipe(icc_pipe2, p, p);
			}
			
			clipped += rv | IMG_WritePixel(out, x, p);
		}
		
		IMG_WriteRow(out, y);
	}
	
	icc_pipe = join_pipes(icc_pipe1, icc_pipe2);
	
	if (verbose && clipped)
		warning("%s: %lu pixels were out of gamut and had been clipped (%4.2g%%)",
			out->file, clipped, 100.0*clipped/out->w/out->h);
	
	if (adjust & ADJ_CURVES) free_matrix(c);
}



/**********************************************************************/

int main(int argc, char *argv[])
{
	char c;
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv-:Pf:i:o:p:r:eO:W:B:G:A:S:C:b:nw:c:")) != -1)
	switch (c) {
	/* General options */
		case 'h':				/* Help message */
			usage(); break;
		case 'v':				/* Verbosity */
			verbose++; break;
		case '-':				/* Long options */
			readopt(&argc, &argv, SELF ".options", optarg);
			break;
		case 'P':				/* Preview */
			preview = 1; break;
	
	/* ICC profile options */
		case 'f':				/* Function */
    			switch (*optarg) {
				case 'f':
					func = icmFwd;
					break;
				case 'b':
					func = icmBwd;
					break;
				case 'g':
					func = icmGamut;
					break;
				case 'p':
					func = icmPreview;
					break;
				default:
					usage();
			}
			break;
		case 'i':				/* Intent */
    			switch (*optarg) {
				case 'p':
					intent = icPerceptual;
					break;
				case 'r':
					intent = icRelativeColorimetric;
					break;
				case 's':
					intent = icSaturation;
					break;
				case 'a':
					intent = icAbsoluteColorimetric;
					break;
				case 'x':
					intent = icmAbsoluteColorimetricXYZ;
					break;
				default:
					usage();
			}
			break;
		case 'o':				/* Search order */
    			switch (*optarg) {
				case 'n':
					order = icmLuOrdNorm;
					break;
				case 'r':
					order = icmLuOrdRev;
					break;
				default:
					usage();
			}
			break;
		case 'p':				/* Profile */
			icc_pipe = add_profile_to_pipe(optarg, NULL, icc_pipe, 1);
			break;
		case 'r':				/* Profile in reverse */
			{
				icmLookupFunc current_func = func;
				
				func = icmBwd;
				icc_pipe = add_profile_to_pipe(optarg, NULL, icc_pipe, 1);
				func = current_func;
			}
			break;
		case 'e':				/* Embed profile */
			embed = 1; break;
	
	/* Color adjustment options */
		case 'O':				/* Adjustment order */
    			switch (*optarg) {
				case 'n':
					preadjust = 0;
					break;
				case 'r':
					preadjust = 1;
					break;
				default:
					usage();
			}
			break;
		case 'W':				/* White point */
			{
				double t, *v;
				int i, set = 0, n;
				char *p = optarg, *q;
				
				if (*p == '=') { set = 1; p++; auto_adjust &= ~ADJ_WHITE; }
				
				if (!strchr(p, ',')) {
					n = 1; v = &wpt;
					adjust |= ADJ_WHITE_SCALAR;
				} else {
					n = 3; v = wptv;
					adjust |= ADJ_WHITE_VECTOR;
				}
				
				for (i = 0; i < n; i++) {
					t = strtod(p, &q);
					if (q == p) usage();
					if (*q == ',') q++; p = q;
					
					v[i] = set ? t : v[i]+t;
				}
			}
			break;
		case 'B':				/* Black point */
			{
				double t, *v;
				int i, set = 0, n;
				char *p = optarg, *q;
				
				if (*p == '=') { set = 1; p++; auto_adjust &= ~ADJ_BLACK; }
				
				if (!strchr(p, ',')) {
					n = 1; v = &bpt;
					adjust |= ADJ_BLACK_SCALAR;
				} else {
					n = 3; v = bptv;
					adjust |= ADJ_BLACK_VECTOR;
				}
				
				for (i = 0; i < n; i++) {
					t = strtod(p, &q);
					if (q == p) usage();
					if (*q == ',') q++; p = q;
					
					v[i] = set ? t : v[i]+t;
				}
			}
			break;
		case 'G':				/* Contrast grade */
			{
				double t, *v;
				int i, set = 0, n;
				char *p = optarg, *q;
				
				if (*p == '=') { set = 1; p++; auto_adjust &= ~ADJ_CONTRAST; }
				
				if (!strchr(p, ',')) {
					n = 1; v = &grade;
					adjust |= ADJ_CONTRAST_SCALAR;
				} else {
					n = 3; v = gradev;
					adjust |= ADJ_CONTRAST_VECTOR;
				}
				
				for (i = 0; i < n; i++) {
					t = strtod(p, &q);
					if (q == p) usage();
					if (*q == ',') q++; p = q;
					
					v[i] = set ? t : v[i]+t;
				}
			}
			break;
		case 'A':				/* Auto-adjust */
			while (*optarg) switch (*(optarg++)) {
				case 'w':
					auto_adjust |= ADJ_WHITE_SCALAR;
					break;
				case 'W':
					auto_adjust |= ADJ_WHITE_VECTOR;
					break;
				case 'b':
					auto_adjust |= ADJ_BLACK_SCALAR;
					break;
				case 'B':
					auto_adjust |= ADJ_BLACK_VECTOR;
					break;
				case 'g':
					auto_adjust |= ADJ_CONTRAST_SCALAR;
					break;
				case 'G':
					auto_adjust |= ADJ_CONTRAST_VECTOR;
					break;
				case 'L':
					auto_adjust |= ADJ_LEVELS;
					break;
				default:
					usage();
			}
			adjust |= auto_adjust;
			break;
		case 'S':				/* Split contrast */
			{
				char *p = optarg, *q;
				
				split = strtod(p, &q); if (q == p) usage();
				
				if (*q != ',') usage(); p = ++q;
				shadow = strtod(p, &q); if (q == p) usage();
				
				if (*q != ',') usage(); p = ++q;
				highlight = strtod(p, &q); if (q == p) usage();
				
				adjust |= ADJ_CONTRAST_SPLIT;
			}
			break;
		case 'C':				/* Center filter */
			cflength = strtod(optarg, NULL);
			if (cflength == 0.0) usage();
			adjust |= ADJ_CENTER;
			break;
	
	/* Data storage options */
		case 'b':				/* Bits per sample */
			bitspersample = atoi(optarg);
			if (bitspersample != 8 && bitspersample != 16) usage();
			break;
		case 'n':				/* Noise */
			noisify = 1; break;
		case 'w':				/* Weave */
			switch (*optarg) {
				case 'c':
					planarconfig = PLANARCONFIG_CONTIG;
					break;
				case 'i':
					planarconfig = PLANARCONFIG_SEPARATE;
					break;
				default:
					usage();
			}
			break;
		case 'c':				/* Compression */
    			switch (*optarg) {
				case 'l':
					compression = COMPRESSION_LZW;
					{
						char *cp = strchr(optarg, ':');
						if (cp) predictor = atoi(cp+1);
					}
					break;
				case 'z':
					compression = COMPRESSION_DEFLATE;
					{
						char *cp = strchr(optarg, ':');
						if (cp) predictor = atoi(cp+1);
					}
					break;
				case 'p':
					compression = COMPRESSION_PACKBITS;
					break;
				case 'n':
					compression = COMPRESSION_NONE;
					break;
				default:
					usage();
			}
			break;
	
	/* Default response */
		default:
			usage();
	}
	
	/* Open input and output images */
	if (argc != optind+2) usage();
	in = IMG_Open(argv[optind++], "r");
	out = IMG_Open(argv[optind++], "w");
	IMG_CpHeader(in, out);
	
	
	/* Check profile pipeline and setup colorspaces */
	if (icc_pipe) {
		unsigned int space;
		icc *embedded = read_embedded_profile(in);
		
		IMG_GetField(in, TIFFTAG_COLORSPACE, &space);
		
		/* Make sure that colorspaces match */
		if (space != input_space(icc_pipe)) {
			if (!embedded) error("%s: Color space does not match profile [%s]",
				in->file, ColorSpaceSignature2str(input_space(icc_pipe)));
			
			/* Maybe with embedded profile... */
			func = icmFwd; intent = icmDefaultIntent; order = icmLuOrdNorm;
			icc_pipe = add_profile_to_pipe(in->file, embedded, icc_pipe, 0);
			
			if (space != input_space(icc_pipe))
				error("%s: Color space does not match embedded profile [%s]",
					in->file, ColorSpaceSignature2str(input_space(icc_pipe)));
		} else if (embedded) warning("%s: Ignoring embedded profile", in->file);
		
		IMG_SetField(out, TIFFTAG_COLORSPACE, output_space(icc_pipe));
	}
	
	/* Embed profile if requested */
	if (embed) {
		iccProfile *p = icc_pipe;
		
		while (p && p->next) p = p->next;
		
		if (p) write_embedded_profile(out, p->icco);
		else warning("%s: Nothing to embed, no profiles specified", out->file);
	}
	
	/* Print translation info if requested */
	if (verbose) {
		fprintf(stderr, "%s\n", usage_msg[0]);
		dump_pipe(icc_pipe);
	}
	
	
	/* Scale image if requested */
	if (preview) {
		unsigned long w, h, s1, s2;
		
		IMG_GetField(in, TIFFTAG_IMAGEWIDTH, &w);
		IMG_GetField(in, TIFFTAG_IMAGELENGTH, &h);
		
		if (w != 0 && h != 0) {
			s1 = ((w > h ? w : h) - 1)/640 + 1;
			s2 = ((w < h ? w : h) - 1)/480 + 1;
			skip = s1 > s2 ? s1 : s2;
		}
		
		if (verbose && skip != 1)
			printf("%s: Scaling down by a factor of %lu\n", in->file, skip);
		
		IMG_SetField(out, TIFFTAG_IMAGEWIDTH, w/skip);
		IMG_SetField(out, TIFFTAG_IMAGELENGTH, h/skip);
	}
	
	/* Set image data storage options */
	if (bitspersample != -1) if (!IMG_SetField(out, TIFFTAG_BITSPERSAMPLE, bitspersample))
		warning("%s: Setting bits per sample is not supported by image format", out->file);
	if (planarconfig != -1) if (!IMG_SetField(out, TIFFTAG_PLANARCONFIG, planarconfig))
		warning("%s: Setting planar configuration is not supported by image format", out->file);
	if (compression != -1) if (!IMG_SetField(out, TIFFTAG_COMPRESSION, compression) ||
				!IMG_SetField(out, TIFFTAG_PREDICTOR, predictor))
		warning("%s: Setting compression is not supported by image format", out->file);
	
	
	/* Finalize settings and allocate buffers */
	IMG_Alloc(in); IMG_Alloc(out);
	
	if (verbose > 3) {
		IMG_PrintStats(in, stderr);
		IMG_PrintStats(out, stderr);
	}
	
	
	/* Make automatic adjustments if necessary */
	if (adjust && (preadjust ? in->space : out->space) != icSigRgbData) {
		warning("Color adjustments implemented only in RGB; clearing all adjustment flags");
		adjust = auto_adjust = 0;
	}
	
	if (auto_adjust) analyze_image(in);
	if (verbose) dump_adjustment();
	
	
	/* Do color correction */
	if (!icc_pipe && !adjust)
		copy_loop(in, out);
	else if (icc_pipe && !adjust)
		calibration_loop(in, out);
	else
		adjustment_loop(in, out);
	
	
	/* Clean up after ourselves */
	IMG_Close(in); IMG_Close(out);
	icc_pipe = free_profiles(icc_pipe);
	
	return 0;
}
