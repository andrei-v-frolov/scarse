/* $Id: scarse.h,v 1.4 2001/06/27 03:58:57 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Function & external variables declarations.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include "util.h"
#include "spaces.h"
#include "imageio/imageio.h"


/**********************************************************************/

#ifndef __SCARSE_H__
#define __SCARSE_H__

#define DEFAULT_GAMMA 2.2


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


/* Calibration data (ipb.c) */

typedef struct { /* curves data and fit */
	int n;			/* data points */
	double **data;		/* [y,x,dx][n] */
	double *fit;		/* curve fit */
} curve;

typedef struct { /* LUT data and fit */
	int flag;			/* outlier? */
	char *label;			/* patch label */
	/* Measured values */
	double  in[MAXCHANNELS];	/* input */
	double out[MAXCHANNELS];	/* output */
	double var[MAXCHANNELS];	/* variance */
	/* Curve pullback values and variances */
	double DEV[MAXCHANNELS*2];	/* (linearized) device */
	double XYZ[6];			/* PCS */
} datapt;



#endif /* __SCARSE_H__ */
