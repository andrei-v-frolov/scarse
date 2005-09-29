/* $Id: calibrate.c,v 1.7 2005/09/29 06:31:01 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Collect and prepare data for device calibration.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

/* CREDITS:
 *
 */

/* TODO:
 *
 */

#define SELF "calibrate"

#define _GNU_SOURCE
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include "scarse.h"



/******************* Options and defaults *****************************/

/* Usage */
char *program_name = SELF;
char *usage_msg[] = {
	"Input/output device calibrator, Version " VERSION,
	"Author: Andrei Frolov <andrei@phys.ualberta.ca>",
	"",
	"Usage: " SELF " -d [scanner|display|printer] [...] profile.icm",
	"For now, only scanner calibration is fully operational.",
	"",
	" General options:",
	"  -h		print this message",
	"  -v		verbosity (cumulative)",
	"  --longopt	expand option macro 'longopt' as defined in etc/" SELF ".options",
	"",
	" ICC profile generation:",
	"  -D		dump measured data to stdout, instead of running ipb",
	"  -O options	pass options to ICC profile builder, see ipb help",
	"",
	" Calibration target options:",
	"  -l		list all known target types (or batches if type was given by -t)",
	"  -t type	specify target type (or IT8.7 target layout file to parse)",
	"  -b batch	specify target batch (or IT8.7 target data file to parse)",
	"  -g geometry	size and position of target grid (geometry is like in X)",
	"  -i file	import calibration target data from raster image",
	"  -r file	render calibration target data into raster image and exit",
	NULL
};


/* Device to be calibrated */
typedef enum {SCANNER = 0, DISPLAY = 1, PRINTER = 2} DeviceType;

/* Options */
static int debug = 0;
static int verbose = 0;

static DeviceType   device = SCANNER;
static char       *profile = NULL;
static char   *ipb_options = "";

/* Calibration targets */
static char         *ttype = NULL;
static char         *batch = NULL;
static char      *geometry = NULL;
static target         *ctg = NULL;



/**********************************************************************/

/* Output scanner color correction data */
static void scanner_correction(FILE *fp, target *tg)
{
	int i, j; double WPT_CAT[3][3], t[3];
	double *Dmin = tg->data[tg->grayscale[0]].XYZ;
	double *Dmax = tg->data[tg->grayscale[tg->graypts-1]].XYZ;
	
	static char *channel[] = {"red", "green", "blue"};
	
	
	/* Profile info */
	fprintf(fp, "PROFILE: input RGB -> XYZ;\n\n");
	
	/* Media white & black points */
	fprintf(fp, "WHITEPOINT: %12.10g %12.10g %12.10g;\n", Dmin[0], Dmin[1], Dmin[2]);
	fprintf(fp, "BLACKPOINT: %12.10g %12.10g %12.10g;\n\n", Dmax[0], Dmax[1], Dmax[2]);
	
	/* map media white point (Dmin) to standard illuminant as per ICC specs */
	XYZ_CAT(Dmin, XYZ_ILLUM, WPT_CAT);
	
	/* Curves */
	for (j = 0; j < 3; j++) {
		fprintf(fp, "CURVE IN%i (%s):\n", j, channel[j]);
		
		for (i = 0; i < tg->graypts; i++) {
			apply33(WPT_CAT, tg->data[tg->grayscale[i]].XYZ, t);
			
			fprintf(fp, "      %12.10g %12.10g \t%12.10g\n",
				tg->data[tg->grayscale[i]].RGB[j], t[1],
				tg->data[tg->grayscale[i]].RGB[3+j]);
		}
		
		fprintf(fp, ";\n\n");
	}
	
	/* LUT */
	if (tg->pts > tg->graypts) {
		fprintf(fp, "LUT:\n");
		
		for (i = 0; i < tg->pts; i++)
			apply33(WPT_CAT, tg->data[i].XYZ, t);
			
			fprintf(fp, " %4s %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g \t%12.10g %12.10g %12.10g\n",
				tg->data[i].label,
				/* RGB values */
				tg->data[i].RGB[0],
				tg->data[i].RGB[1],
				tg->data[i].RGB[2],
				/* XYZ values */
				t[0], t[1], t[2],
				/* measured deviation */
				tg->data[i].RGB[3],
				tg->data[i].RGB[4],
				tg->data[i].RGB[5]
			);
		
		fprintf(fp, ";\n\n");
	}
}



/**********************************************************************/

/* Open pipe to ICC profile builder backend */
static FILE *pipe2ipb()
{
	FILE *fp;
	char *buffer, *verb;
	char *class[3] = {"i -U", "d -M", "o"};
	
	if (debug) return stdout;
	
	verb = xstrdup("-vvvvv");
	if (verbose) { if (verbose < 5) verb[verbose+1] = 0; } else *verb = 0;
	
	asprintf(&buffer, "ipb %s -c%s -C- %s %s",
		verb, class[device], ipb_options, profile);
	
	if (!(fp = popen(buffer, "w")))
		error("Could not start profile builder, command line:\n\t%s", buffer);
	
	free(buffer); free(verb);
	
	return fp;
}


/* Main routine */
int main(int argc, char *argv[])
{
	char c;
	
	/* RGB setup */
	{
		double primaries[4][2];
		
		LookupPrimaries("Adobe", primaries, NULL);
		SetPrimaries(primaries);
	}
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv-:d:DO:lt:b:g:i:r:")) != -1)
	switch (c) {
	/* General options */
		case 'h':				/* Help message */
			usage(); break;
		case 'v':				/* Verbosity */
			verbose++; break;
		case '-':				/* Long options */
			readopt(&argc, &argv, SELF ".options", optarg);
			break;
	
	/* Device to be calibrated */
		case 'd':
			switch (*optarg) {
				case 's':
					device = SCANNER;
					break;
				case 'd':
				case 'm':
					device = DISPLAY;
					break;
				case 'p':
					device = PRINTER;
					break;
				default:
					usage();
			}
			break;
	
	/* ICC profile builder options */
		case 'D':				/* Data dump */
			debug = 1;
			break;
		case 'O':				/* Pass options to ipb */
			ipb_options = optarg;
			break;
	
	/* Calibration targets */
		case 'l':				/* List targets */
			if (verbose) {
				fprintf(stderr, "%s\n\n", usage_msg[0]);
				if (ttype) fprintf(stderr, "Known %s target batches:\n", ttype);
				else fprintf(stderr, "Supported calibration target types:\n");
				
				fflush(stderr);
			}
			
			list_targets(ttype); exit(0);
			break;
		case 't':				/* Target type */
			ttype = optarg;
			break;
		case 'b':				/* Target batch */
			batch = optarg;
			break;
		case 'g':				/* Grid geometry */
			geometry = optarg;
			break;
		case 'i':				/* Import image */
			{
				char *data, *layout;
				target *tg = xmalloc(sizeof(target));
				
				if (!ttype) ttype = "Q60E3";
				if (!batch) batch = "generic";
				
				if (verbose) {
					fprintf(stderr, "%s\n", usage_msg[0]);
					fprintf(stderr, "Reading %s(%s) calibration target data from '%s'...\n", ttype, batch, optarg);
				}
				
				find_target(ttype, batch, &data, &layout);
				parse_IT87_target(tg, data, layout);
				read_IT87_target(tg, optarg, geometry);
				ctg = tg;
			}
			break;
		case 'r':				/* Render image */
			{
				char *data, *layout;
				target *tg = xmalloc(sizeof(target));
				
				if (!ttype) ttype = "Q60E3";
				if (!batch) batch = "generic";
				
				if (verbose) {
					fprintf(stderr, "%s\n", usage_msg[0]);
					fprintf(stderr, "Rendering %s(%s) calibration target into '%s'...\n", ttype, batch, optarg);
				}
				
				find_target(ttype, batch, &data, &layout);
				parse_IT87_target(tg, data, layout);
				render_IT87_target(tg, optarg, geometry);
				exit(0);
			}
			break;
	
	/* Default response */
		default:
			usage();
	}
	
	if (argc != optind+1) usage();
	profile = argv[optind++];
	
	
	/* Output correction data */
	switch (device) {
		case SCANNER:
			if (!ctg) error("No calibration target given, nothing to do...");
			if (verbose) fprintf(stderr, "Generating scanner profile '%s' from collected data\n", profile);
			{ FILE *ipb = pipe2ipb(); scanner_correction(ipb, ctg); pclose(ipb); }
			break;
		case DISPLAY:
			usage();
			break;
		case PRINTER:
			usage();
			break;
		default:
			usage();
	}
	
	return 0;
}
