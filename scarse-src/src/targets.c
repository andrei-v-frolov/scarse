/* $Id: targets.c,v 1.5 2001/02/27 04:00:53 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * IT8.7 calibration target data handling routines.
 * 
 * Copyright (C) 2000 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#define _GNU_SOURCE

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fnmatch.h>
#include <math.h>

#include <icc.h>
#include "scarse.h"


/****************** Target file indexing ******************************/

/* List target types and batches found in the target index */
void list_targets(char *type)
{
	char *t, *b, *f;
	size_t bsize = 128;
	char *buffer = (char *)xmalloc(bsize);
	FILE *fp = xfetch("targets", "TARGETS", "r"), *out = popen("sort | uniq", "w");
	
	if (!out) error("Could not run sort on target list");
	
	while (getline(&buffer, &bsize, fp) != -1) {
		/* skip over comments and empty lines */
		if (*buffer == '#' || *buffer == '\n' || *buffer == '\r') continue;
		
		/* index file format: type batch file */
		if (sscanf(buffer, "%as %as %as", &t, &b, &f) != 3)
			error("Error parsing target index file");
		
		if (type) {
			/* list all target batches for a given type */
			if (!strcmp(t, type) && !strchr(b, '*')) fprintf(out, "%s\n", b);
		} else {
			/* list all available target types */
			if (!strchr(t, '*') && !strchr(b, '*')) fprintf(out, "%s\n", t);
		}
		
		free(t);
		free(b);
		free(f);
	}
	
	free(buffer);
	pclose(out);
	fclose(fp);
}

/* Find target data and layout files in the target index */
void find_target(char *type, char *batch, char **data, char **layout)
{
	char *t, *b, *f;
	size_t bsize = 128;
	char *buffer = (char *)xmalloc(bsize);
	FILE *fp = xfetch("targets", "TARGETS", "r");
	
	*data = *layout = NULL;
	
	while (getline(&buffer, &bsize, fp) != -1) {
		/* skip over comments and empty lines */
		if (*buffer == '#' || *buffer == '\n' || *buffer == '\r') continue;
		
		/* index file format: type batch file */
		if (sscanf(buffer, "%as %as %as", &t, &b, &f) != 3)
			error("Error parsing target index file");
		
		/* the first match is it */
		if (!(*data) && !strcmp(t, type) && !strcmp(b, batch)) *data = xstrdup(f);
		if (!(*layout) && !fnmatch(t, type, 0) && !fnmatch(b, batch, 0)) *layout = xstrdup(f);
		
		free(t);
		free(b);
		free(f);
	}
	
	/* if not found in index, assume type and batch are layout and data files respectively */
	if (!(*data)) *data = xstrdup(batch);
	if (!(*layout)) *layout = xstrdup(type);
	
	free(buffer);
	fclose(fp);
}



/****************** IT8.7 data parser *********************************/

/* IT8.7 data parser is implemented as a finite state machine */
typedef enum {UNDEF, OK, FORMAT, DATA, LAYOUT} ParserState;

static char *state_id[] = /* state id strings for error reporting */
	{"undef", "header", "data format", "data", "target layout"};

/* Parser state */
static ParserState state = UNDEF;

/* Input stream */
static char token[65];
static char *it87file = "";
static FILE *source = NULL;

/* IT8.7 data */
static char        *orig = NULL;
static char        *desc = NULL;
static char     *created = NULL;
static char    *manufact = NULL;
static char        *date = NULL;
static char      *serial = NULL;
static char    *material = NULL;

static int           pts = 0;
static int         count = 0;
static int        fields = 0;
static char     **format = NULL;

static struct _data { char *label; double *values; } *data = NULL;

/* Layout data */
static int          rows = 0;
static int          cols = 0;
static int      *rowspan = NULL;
static int      *colspan = NULL;
static char    ***layout = NULL;

static int       graypts = 0;
static char  **grayscale = NULL;


/* Reset parser to initial state, free allocated memory */
static void reset_parser()
{
	int i, j;
	
	#define INITPTR(p) {if (p) free(p); p = NULL;}
	
	state = UNDEF;
	
	INITPTR(orig);
	INITPTR(desc);
	INITPTR(created);
	INITPTR(manufact);
	INITPTR(date);
	INITPTR(serial);
	INITPTR(material);
	
	if (format) for (i = 1; i < fields; i++)
		INITPTR(format[i]);
	INITPTR(format);
	
	if (data) for (i = 0; i < pts; i++) {
		INITPTR(data[i].label);
		INITPTR(data[i].values);
	}
	INITPTR(data);
	
	pts = count = fields = 0;
	
	INITPTR(rowspan);
	INITPTR(colspan);
	
	if (layout) for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) INITPTR(layout[i][j]);
		INITPTR(layout[i]);
	}
	INITPTR(layout);
	
	if (grayscale) for (i = 1; i < graypts; i++)
		INITPTR(grayscale[i]);
	INITPTR(grayscale);
	
	rows = cols = graypts = 0;
	
	#undef INITPTR
}



/****************** Low-level parser functions: ***********************/

/* Report unexpected parser input */
static void unexpected(char *expected)
{
	int l;
	char buffer[33];
	
	fgets(buffer, 33, source); l = strlen(buffer);
	if (l > 0 && buffer[l-1] == '\n') buffer[--l] = 0;
	if (l == 32) buffer[31] = buffer[30] = buffer[29] = '.';
	
	error("%s: error parsing %s tag: expected %s, got '%s'\n",
		it87file, token, expected, buffer);
}


/* Parse integer */
static void parse_int(void *i)
{
	int t;
	
	if (fscanf(source, "%i", i ? (int *)i : &t) != 1) unexpected("integer");
}

/* Parse double */
static void parse_double(void *d)
{
	double t;
	
	if (fscanf(source, "%lf", d ? (double *)d : &t) != 1) unexpected("double");
}

/* Parse word */
static void parse_word(void *w)
{
	char *t = NULL;
	
	if (fscanf(source, "%as", w ? (char **)w : &t) != 1) unexpected("bare word");
	
	if (t) free(t);
}

/* Parse quote-delimetered string */
static void parse_string(void *s)
{
	char *t = NULL;
	
	if (fscanf(source, " \"%a[^\"\n]\"", s ? (char **)s : &t) != 1)
		unexpected("quote-delimitered string");
	
	if (t) free(t);
}

/* Parse format specification */
static void parse_format(void *ptr)
{
	int i;
	char **t = (char **)xmalloc(fields*sizeof(char *));
	
	if (!fields) error("%s: no data fields declared - missing NUMBER_OF_FIELDS tag?", it87file);
	
	t[0] = "ID";
	
	for (i = 1; i < fields; i++)
		parse_word(&(t[i]));
	
	*((char ***)(ptr)) = t;
}

/* Parse calibration data */
static void parse_data(void *ptr)
{
	int i;
	struct _data **dp = (struct _data **)(ptr), *d = *dp;
	
	if (!d) {
		if (!fields) error("%s: no data fields declared - missing NUMBER_OF_FIELDS tag?", it87file);
		if (!pts) error("%s: no data points declared - missing NUMBER_OF_SETS tag?", it87file);
		count = 0; d = *dp = (struct _data *)xmalloc(pts*sizeof(struct _data));
	}
	
	if (count >= pts) error("%s: excess data encountered (%i pts declared)", it87file, pts);
	
	d[count].label = xstrdup(token);
	d[count].values = vector(fields);
	
	for (i = 1; i < fields; i++)
		parse_double(&(d[count].values[i]));
	
	count++;
}

/* Parse row or column span for target layout */
static void parse_span(void *ptr)
{
	int i;
	int n = (ptr == &rowspan) ? rows : cols;
	int *t = (int *)xmalloc(n*sizeof(int));
	
	if (!n) error("%s: no target %s declared - missing NUMBER_OF_%s tag?", it87file,
		(ptr == &rowspan) ? "rows" : "columns", (ptr == &rowspan) ? "ROWS" : "COLUMNS");
	
	for (i = 0; i < n; i++)
		parse_int(&(t[i]));
	
	*((int **)(ptr)) = t;
}

/* Parse target layout */
static void parse_layout(void *ptr)
{
	int i, j;
	char ***t = (char ***)xmalloc(rows*sizeof(char **));
	
	if (!rows) error("%s: no target rows declared - missing NUMBER_OF_ROWS tag?", it87file);
	if (!cols) error("%s: no target columns declared - missing NUMBER_OF_COLUMNS tag?", it87file);
	
	for (i = 0; i < rows; i++) {
		t[i] = (char **)xmalloc(cols*sizeof(char *));
		
		for (j = 0; j < cols; j++)
			parse_word(&t[i][j]);
	}
	
	*((char ****)(ptr)) = t;
}

/* Parse grayscale for target layout */
static void parse_scale(void *ptr)
{
	int i;
	char **t = (char **)xmalloc(graypts*sizeof(char *));
	
	if (!graypts) error("%s: no grayscale points declared - missing GRAYSCALE_STEPS tag?", it87file);
	
	for (i = 0; i < graypts; i++)
		parse_word(&(t[i]));
	
	*((char ***)(ptr)) = t;
}



/****************** Parser transitions: *******************************/

typedef struct {
	char *tag; ParserState from, to;
	void (*action)(void *arg); void *arg;
} ParserTransitions;

/* Transition table */
static ParserTransitions trans_table[] = {
	{"IT8.7/*",		UNDEF, OK,	NULL, NULL},
	{"ORIGINATOR",		OK, OK,		parse_string, &orig},
	{"DESCRIPTOR",		OK, OK,		parse_string, &desc},
	{"CREATED",		OK, OK,		parse_string, &created},
	{"MANUFACTURER",	OK, OK,		parse_string, &manufact},
	{"PROD_DATE",		OK, OK,		parse_string, &date},
	{"SERIAL",		OK, OK,		parse_string, &serial},
	{"MATERIAL",		OK, OK,		parse_string, &material},
	{"KEYWORD",		OK, OK,		parse_string, NULL},
	{"NUMBER_OF_SETS",	OK, OK,		parse_int, &pts},
	{"NUMBER_OF_FIELDS",	OK, OK,		parse_int, &fields},
	{"BEGIN_DATA_FORMAT",	OK, FORMAT,	NULL, NULL},
	{"END_DATA_FORMAT",	FORMAT, OK,	NULL, NULL},
	{"SAMPLE_ID",		FORMAT, FORMAT,	parse_format, &format},
	{"BEGIN_DATA",		OK, DATA,	NULL, NULL},
	{"END_DATA",		DATA, OK,	NULL, NULL},
	{"*",			DATA, DATA,	parse_data, &data},
	/* Extended tags for target layout specification */
	{"BEGIN_LAYOUT",	OK, LAYOUT,	NULL, NULL},
	{"END_LAYOUT",		LAYOUT, OK,	NULL, NULL},
	{"NUMBER_OF_ROWS",	LAYOUT, LAYOUT,	parse_int, &rows},
	{"NUMBER_OF_COLUMNS",	LAYOUT, LAYOUT,	parse_int, &cols},
	{"ROW_SPAN",		LAYOUT, LAYOUT,	parse_span, &rowspan},
	{"COLUMN_SPAN",		LAYOUT, LAYOUT,	parse_span, &colspan},
	{"LAYOUT_TABLE",	LAYOUT, LAYOUT,	parse_layout, &layout},
	{"GRAYSCALE_STEPS",	LAYOUT, LAYOUT,	parse_int, &graypts},
	{"GRAYSCALE",		LAYOUT, LAYOUT,	parse_scale, &grayscale},
	NULL
};


/* Parse a token from input stream */
static void parse_token()
{
	int i;
	
	if (fscanf(source, "%64s", token) != 1) return;
	if (*token == '#') { fscanf(source, "%*[^\n]\n"); return; }
	
	for (i = 0; trans_table[i].tag; i++) {
		if (trans_table[i].from != state) continue;
		if (fnmatch(trans_table[i].tag, token, FNM_CASEFOLD)) continue;
		
		/* Found a match */
		state = trans_table[i].to;
		if (trans_table[i].action)
			(*trans_table[i].action)(trans_table[i].arg);
		
		return;
	}
	
	/* No valid transition */
	if (state == UNDEF)
		error("%s: missing magic 'IT8.7' tag (got '%s' instead) - wrong file type?\n", it87file, token);
	else
		error("%s: invalid tag '%s' encountered while parsing %s", it87file, token, state_id[state]);
}



/****************** High-level parser routines ************************/

/* Rescale XYZ from 100.0 (IT8.7) to 1.0 (us) */
static void scale100to1(double in[], double out[])
{
	out[0] = in[0]/100.0;
	out[1] = in[1]/100.0;
	out[2] = in[2]/100.0;
}

/* Convert Yxy to XYZ and scale from 100.0 (IT8.7) to 1.0 (us) */
static void Yxy2XYZ_100to1(double in[], double out[])
{
	double XYZ[3];
	
	Yxy2XYZ(in, XYZ);
	scale100to1(XYZ, out);
}

/* Convert Gray to XYZ and scale from 100.0 (IT8.7) to 1.0 (us) */
static void Gray2XYZ_100to1(double in[], double out[])
{
	double XYZ[3];
	
	Gray2XYZ(in, XYZ);
	scale100to1(XYZ, out);
}

/* find a field in format specification */
static int field_idx(const char *id, int strict)
{
	int i;
	
	if (strict) {
		for (i = 1; i < fields; i++)
			if (!strcmp(format[i], id)) return i;
	} else {
		for (i = 1; i < fields; i++)
			if (!strcasecmp(format[i], id)) return i;
	}
	
	return 0;
}

/* figure out how to convert supplied data to XYZ */
static transform fields2XYZ(int idx[])
{
	/* Lab is less susceptible to roundoff errors, so try it first */
	if (	(idx[0] = field_idx("Lab_L", 0)) &&
		(idx[1] = field_idx("Lab_a", 0)) &&
		(idx[2] = field_idx("Lab_b", 0)) ) return Lab2XYZ;
	
	if (	(idx[0] = field_idx("Luv_L", 0)) &&
		(idx[1] = field_idx("Luv_u", 0)) &&
		(idx[2] = field_idx("Luv_v", 0)) ) return Luv2XYZ;
	
	/* if not Lab, there must be XYZ data around somewhere */
	/* remember, IT8.7 has XYZ scaled to 100.0, not 1.0 like us! */
	if (	(idx[0] = field_idx("XYZ_X", 0)) &&
		(idx[1] = field_idx("XYZ_Y", 0)) &&
		(idx[2] = field_idx("XYZ_Z", 0)) ) return scale100to1;
	
	if (	(idx[0] = field_idx("Yxy_Y", 1)) &&
		(idx[1] = field_idx("Yxy_x", 1)) &&
		(idx[2] = field_idx("Yxy_y", 1)) ) return Yxy2XYZ_100to1;
	
	/* OK, OK, so maybe it's a grayscale... */
	if ( (idx[0] = idx[1] = idx[2] = field_idx("GRAY", 0)) ) return Gray2XYZ_100to1;
	
	/* Nah, data format is too weird! */
	error("%s: Don't know how to convert supplied data to XYZ\n", it87file);
	
	return NULL;
}

/* find a label in calibration data */
static int label_idx(const char *id)
{
	int i;
	
	if (!strcmp("*", id)) return -1;
	
	for (i = 0; i < pts; i++)
		if (!strcmp(data[i].label, id)) return i;
	
	error("%s: label %s from layout specification not found in calibration data", it87file, id);
	
	return -1;
}


/* Parse IT8.7 data file and store the data in the target structure */
void parse_IT87_target(target *tg, char *data_file, char *layout_file)
{
	/* parse IT8.7 data file */
	it87file = data_file;
	source = xfetch("targets", data_file, "r");
	while (!feof(source) && !ferror(source)) parse_token();
	if (state != OK) warning("%s: IT8.7 parser left in intermediate state - missing END_... tag?", it87file);
	
	fclose(source);
	state = UNDEF;
	
	/* store parsed calibration data into target structure */
	{
		double in[3];
		int i, j, idx[3];
		transform data2XYZ;
		
		/* sanity checks */
		if (!data) error("%s: no calibration data found", it87file);
		if (count != pts) error("%s: missing data (%i pts declared, %i parsed)", it87file, pts, count);
		
		if (format) data2XYZ = fields2XYZ(idx);
		else if (fields > 3) {
			warning("%s: no calibration data format found - assuming straight XYZ values", it87file);
			
			data2XYZ = scale100to1;
			idx[0] = 1; idx[1] = 2; idx[2] = 3;
		} else error("%s: no calibration data format found", it87file);
		
		tg->pts = pts;
		tg->data = (struct _target_data *)xmalloc(pts*sizeof(struct _target_data));
		
		for (i = 0; i < pts; i++) {
			tg->data[i].label = xstrdup(data[i].label);
			
			for (j = 0; j < 3; j++) in[j] = data[i].values[idx[j]];
			(*data2XYZ)(in, tg->data[i].XYZ);
		}
	}
	
	/* parse IT8.7 layout file */
	it87file = layout_file;
	source = xfetch("targets", layout_file, "r");
	while (!feof(source) && !ferror(source)) parse_token();
	if (state != OK) warning("%s: IT8.7 parser left in intermediate state - missing END_... tag?", it87file);
	
	fclose(source);
	state = UNDEF;
	
	/* store parsed layout data into target structure */
	{
		int i, j, k;
		
		/* sanity checks */
		if (!layout) error("%s: no target layout found", it87file);
		if (!grayscale) error("%s: no grayscale definition found", it87file);
		
		/* init layout table */
		tg->rows = rows;
		tg->cols = cols;
		
		tg->layout = (int **)xmalloc(rows*sizeof(int *));
		for (i = 0; i < rows; i++) {
			tg->layout[i] = (int *)xmalloc(cols*sizeof(int));
			
			for (j = 0; j < cols; j++)
				tg->layout[i][j] = label_idx(layout[i][j]);
		}
		
		/* init subrow & subcolumn accelerator tables */
		for (tg->subrows = 0, i = 0; i < rows; i++)
			tg->subrows += rowspan ? rowspan[i] : 1;
		tg->subridx = (int *)xmalloc(tg->subrows*sizeof(int));
		
		for (i = 0, j = 0; i < rows; i++)
			for (k = 0; k < (rowspan ? rowspan[i] : 1); k++)
				tg->subridx[j++] = i;
		
		for (tg->subcols = 0, i = 0; i < cols; i++)
			tg->subcols += colspan ? colspan[i] : 1;
		tg->subcidx = (int *)xmalloc(tg->subcols*sizeof(int));
		
		for (i = 0, j = 0; i < cols; i++)
			for (k = 0; k < (colspan ? colspan[i] : 1); k++)
				tg->subcidx[j++] = i;
		
		/* init grayscale */
		tg->graypts = graypts;
		tg->grayscale = (int *)xmalloc(graypts*sizeof(int));
		for (i = 0; i < graypts; i++) {
			tg->grayscale[i] = label_idx(grayscale[i]);
			
			if (tg->grayscale[i] < 0)
				error("%s: void data labels are not allowed in grayscale specification", it87file);
		}
	}
	
	/* Clean up after ourselves */
	reset_parser();
}



/****************** Reading and writing targets ***********************/

/* Parse geometry specification of the data grid */
static void parse_geometry(char *geometry, double *left, double *top, double *right, double *bottom)
{
	int confused = 0;
	char *p = geometry, *q;
	double w = 1.0, h = 1.0, x = 0.0, y = 0.0;
	
	if (!geometry) {
		*left = *top = 0.0;
		*right = *bottom = 1.0;
		return;
	}
	
	/* Parse geometry string */
	w = strtod(p, &q)/100.0;   if (q == p || *q != 'x') confused = 1;
	h = strtod(++q, &p)/100.0; if (q == p) confused = 1;
	
	if (*p) { /* displacement is optional */
		x = strtod(p, &q)/100.0; if (q == p) confused = 1;
		y = strtod(q, &p)/100.0; if (q == p) confused = 1;
	}
	
	if (confused)
		error("Invalid geometry specification: '%s'", geometry);
	
	if (x < 0.0) { *left = 1.0 + x - w; *right = 1.0 + x; }
	else { *left = x; *right = x + w; }
	
	if (y < 0.0) { *top = 1.0 + y - h; *bottom = 1.0 + y; }
	else { *top = y; *bottom = y + h; }
}


/* Render IT8.7 target data into raster image */
void render_IT87_target(target *tg, char *file, char *geometry)
{
	int r, c;
	unsigned long x, y, w, h;
	image *img = IMG_Open(file, "w");
	double **Lab = matrix(tg->pts+1, 3) + 1;
	
	/* Translate target data to Lab colorspace */
	{
		int i;
		double t[3], *Dmin = tg->data[tg->grayscale[0]].XYZ;
		
		/* background pixel - make it neutral gray */
		Lab[-1][0] = 50.0; Lab[-1][1] = Lab[-1][2] = 0.0;
		
		for (i = 0; i < tg->pts; i++) {
			/* map Dmin to white point as per ICC specs */
			t[0] = tg->data[i].XYZ[0]/Dmin[0]*XYZ_WPT[0];
			t[1] = tg->data[i].XYZ[1]/Dmin[1]*XYZ_WPT[1];
			t[2] = tg->data[i].XYZ[2]/Dmin[2]*XYZ_WPT[2];
			
			XYZ2Lab(t, Lab[i]);
		}
	}
	
	/* Setup image size */
	w = tg->subcols << 5; h = tg->subrows << 5;
	
	if (geometry) {
		double l, t, r, b;
		
		parse_geometry(geometry, &l, &t, &r, &b);
		w = 100.0*(r-l); h = 100.0*(b-t);
	}
	
	IMG_SetField(img, TIFFTAG_COLORSPACE, icSigLabData);
	IMG_SetField(img, TIFFTAG_IMAGEWIDTH, w);
	IMG_SetField(img, TIFFTAG_IMAGELENGTH, h);
	
	IMG_Alloc(img);
	
	
	/* Render target */
	for (y = 0; y < h; y++) {
		r = tg->subridx[y*tg->subrows/h];
		
		for (x = 0; x < w; x++) {
			c = tg->subcidx[x*tg->subcols/w];
			
			IMG_WritePixel(img, x, Lab[tg->layout[r][c]]);
		}
		
		IMG_WriteRow(img, y);
	}
	
	free_matrix(Lab-1);
	IMG_Close(img);
}


/* Read IT8.7 target data from raster image */
void read_IT87_target(target *tg, char *file, char *geometry)
{
	int i, r, c;
	int *xidx, *yidx;
	unsigned long x, y;
	double p[3], deviation;
	image *img = IMG_Open(file, "r");
	struct { unsigned long l, r, t, b, w, h; } crop;
	double ***data = (double ***)xmalloc((tg->pts)*sizeof(double **));
	unsigned long *pxls = uvector(tg->pts), *size = uvector(tg->pts);
	
	
	IMG_Alloc(img);
	
	/* Sanity check */
	if (img->space != icSigRgbData)
		error("%s: RGB color space expected", img->file);
	
	/* Initialize data structures */
	for (i = 0; i < tg->pts; i++) {
		pxls[i] = 0;
		size[i] = 1 << 12;
		data[i] = matrix(3, size[i]);
	}
	
	/* Image cropping */
	{
		double l, t, r, b;
		
		parse_geometry(geometry, &l, &t, &r, &b);
		
		crop.l = l * img->w;
		crop.r = r * img->w;
		crop.t = t * img->h;
		crop.b = b * img->h;
		
		crop.w = crop.r - crop.l;
		crop.h = crop.b - crop.t;
	}
	
	/* Column lookup table */
	xidx = ivector(img->w); for (x = 0; x < img->w; x++) {
		if (x < crop.l || x >= crop.r) { xidx[x] = -1; continue; }
		
		c = (x-crop.l)*tg->subcols/crop.w; xidx[x] = tg->subcidx[c];
		
		if ( (c == 0 || tg->subcidx[c] != tg->subcidx[c-1]) &&
			(x < crop.l + (4*c+1) * crop.w/tg->subcols/4) ) xidx[x] = -1;
		if ( (c == tg->subcols || tg->subcidx[c] != tg->subcidx[c+1]) &&
			(x > crop.l + (4*c+3) * crop.w/tg->subcols/4) ) xidx[x] = -1;
	}
	
	/* Row lookup table */
	yidx = ivector(img->h); for (y = 0; y < img->h; y++) {
		if (y < crop.t || y >= crop.b) { yidx[y] = -1; continue; }
		
		r = (y-crop.t)*tg->subrows/crop.h; yidx[y] = tg->subridx[r];
		
		if ( (r == 0 || tg->subridx[r] != tg->subridx[r-1]) &&
			(y < crop.t + (4*r+1) * crop.h/tg->subrows/4) ) yidx[y] = -1;
		if ( (r == tg->subrows || tg->subridx[r] != tg->subridx[r+1]) &&
			(y > crop.t + (4*r+3) * crop.h/tg->subrows/4) ) yidx[y] = -1;
	}
	
	
	/* Read target data */
	for (y = 0; y < img->h; y++) {
		IMG_ReadRow(img, y); /* serialize access */
		r = yidx[y]; if (r < 0) continue;
		
		for (x = 0; x < img->w; x++) {
			c = xidx[x]; if (c < 0) continue;
			i = tg->layout[r][c]; if (i < 0) continue;
			
			IMG_ReadPixel(img, x, p);
			
			if (pxls[i] >= size[i])
				data[i] = grow_matrix(data[i], 3, size[i] <<= 1);
			
			data[i][0][pxls[i]] = p[0];
			data[i][1][pxls[i]] = p[1];
			data[i][2][pxls[i]] = p[2];
			
			pxls[i]++;
		}
	}
	
	free_vector(xidx);
	free_vector(yidx);
	
	
	/* Process target data */
	for (i = 0; i < tg->pts; i++) {
		if (!pxls[i]) error("%s: No data at point %s", img->file, tg->data[i].label);
		
		/* Average values */
		tg->data[i].RGB[0] = avg(pxls[i], data[i][0]);
		tg->data[i].RGB[1] = avg(pxls[i], data[i][1]);
		tg->data[i].RGB[2] = avg(pxls[i], data[i][2]);
		
		/* Deviations */
		tg->data[i].RGB[3] = stddev(pxls[i], data[i][0], tg->data[i].RGB[0]);
		tg->data[i].RGB[4] = stddev(pxls[i], data[i][1], tg->data[i].RGB[1]);
		tg->data[i].RGB[5] = stddev(pxls[i], data[i][2], tg->data[i].RGB[2]);
		
		deviation = sqrt( (
			tg->data[i].RGB[3] * tg->data[i].RGB[3] +
			tg->data[i].RGB[4] * tg->data[i].RGB[4] +
			tg->data[i].RGB[5] * tg->data[i].RGB[5])/3.0
		);
		
		if (deviation > 0.05)
			warning("%s: Fluctuations are big in box %s (SDev = %4.3g)",
				img->file, tg->data[i].label, deviation);
		
		free_matrix(data[i]);
	}
	
	free_vector(pxls);
	free_vector(size);
	free(data);
	
	IMG_Close(img);
}
