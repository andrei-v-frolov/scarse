/* $Id: ppm.c,v 1.1.1.1 2001/01/26 22:45:32 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * PPM image format wrappers (read-only).
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * PPM wrapper by Peter Kirchgessner <peter@kirchgessner.net>
 *
 */

/* CREDITS:
 *
 */

/* TODO:
 *   Write access would be nice for completeness...
 */

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <icc.h>

#include "imageio.h"
#include "../util.h"


/* Private data for PPM image wrappers */
typedef struct { short bits; unsigned long imgstart, r; } ppm_private;



/******************* Reading/writing tags *****************************/

/* Get TIFF image tag */
static int PPM_IMG_GetField(image *img, ttag_t tag, void *arg)
{
	int r = 1;
	
	switch (tag) {
		case TIFFTAG_COLORSPACE:
			*((unsigned int *)(arg)) = img->space; break;
		case TIFFTAG_IMAGEWIDTH:
			*((uint32 *)(arg)) = img->w; break;
		case TIFFTAG_IMAGELENGTH:
			*((uint32 *)(arg)) = img->h; break;
		case TIFFTAG_BITSPERSAMPLE:
			*((uint16 *)(arg)) = ((ppm_private *)(img->private))->bits; break;
		case TIFFTAG_PHOTOMETRIC:
			*((uint16 *)(arg)) = PHOTOMETRIC_RGB; break;
		case TIFFTAG_SAMPLESPERPIXEL:
			*((uint16 *)(arg)) = 3; break;
		case TIFFTAG_PLANARCONFIG:
			*((uint16 *)(arg)) = PLANARCONFIG_CONTIG; break;
		case TIFFTAG_ORIENTATION:
			*((uint16 *)(arg)) = ORIENTATION_TOPLEFT; break;
		default:
			r = 0; break;
	}
	
	return r;
}



/******************* Row access ***************************************/

/* Read row into buffer from binary PPM file */
static void PPM_IMG_GetRow_binary(image *img, unsigned long r)
{
	unsigned long current_row = ((ppm_private *)(img->private))->r;
	unsigned long imgstart = ((ppm_private *)(img->private))->imgstart;
	unsigned long l = (img->w)*3*(((ppm_private *)(img->private))->bits)/8;
	
	if (r == current_row) return;
	
	/* Check for random access and reposition if possible */
	if (r != current_row+1)
		if (fseek((FILE *)(img->fd), imgstart + r*l, 0))
			error("%s: Random access is not supported on this stream", img->file);
	
	/* Read row data */
	if (fread(img->buffer, 1, l, (FILE *)(img->fd)) != l)
		error("%s: Error reading row %ul", img->file, r);
	
	((ppm_private *)(img->private))->r = r;
}

/* Read row into buffer from ASCII PPM file */
static void PPM_IMG_GetRow_ascii(image *img, unsigned long r)
{
	int v;
	unsigned long i, l = (img->w)*3;
	unsigned long current_row = ((ppm_private *)(img->private))->r;
	unsigned long imgstart = ((ppm_private *)(img->private))->imgstart;
	
	if (r == current_row) return;
	
	/* Check for random access and reposition if possible */
	if (r != current_row+1) {
		if (fseek((FILE *)(img->fd), imgstart, 0))
			error("%s: Random access is not supported on this stream", img->file);
		
		for (i = 0; i < r*l; i++) if (fscanf((FILE *)(img->fd), "%i", &v) != 1)
			error("%s: Error reading ASCII PPM image data", img->file);
	}
	
	/* Read row data */
	for (i = 0; i < l; i++) {
		if (fscanf((FILE *)(img->fd), "%i", &v) != 1)
			error("%s: Error reading ASCII PPM image data", img->file);
		
		if (((ppm_private *)(img->private))->bits == 16)
			((unsigned short *)(img->buffer))[i] = (unsigned short)v;
		else
			((unsigned char *)(img->buffer))[i] = (unsigned char)v;
	}
	
	((ppm_private *)(img->private))->r = r;
}

/* Free row buffer */
static void PPM_IMG_FreeRow(image *img)
{
	free(img->buffer);
}



/******************* Image stats **************************************/

/* Print PPM image stats */
static void PPM_IMG_PrintStats(image *img, FILE *fp)
{
	fprintf(fp, "\n%s: ", img->file);
	fprintf(fp, "PPM data at offset 0x%lx\n", ((ppm_private *)(img->private))->imgstart);
	fprintf(fp, "  Image Width: %lu Image Length: %lu\n", img->w, img->h);
	fprintf(fp, "  Bits/Sample: %i\n", ((ppm_private *)(img->private))->bits);
}



/******************* PPM image input/output ***************************/

/* Allocate buffers */
static void PPM_IMG_Alloc(image *img)
{
	unsigned long l = (img->w)*3*(((ppm_private *)(img->private))->bits)/8;
	
	img->buffer = xmalloc((size_t)(l));
}

/* Close PPM image */
static void PPM_IMG_Close(image *img)
{
	IMG_Clean(img);
	
	free(img->file);
	fclose((FILE *)(img->fd));
	free(img->private);
	
	free(img);
}

/* Open PPM image */
image *PPM_IMG_Open(const char *name, const char *mode)
{
	FILE *fp;
	size_t bsize = 128;
	char *buffer = (char *)xmalloc(bsize);
	image *img = (image *)xmalloc(sizeof(image));
	int ascii, binary, max, bits;
	
	
	if (*mode != 'r') error("%s: Writing PPM is not supported", name);
	
	img->file = xstrdup(name);
	
	/* Open file and check PPM signature */
	fp = fopen(name, mode); img->fd = (void *)fp;
	if (!fp) error("Can't open file '%s'", name);
	
	if (getline(&buffer, &bsize, fp) == -1)
		error("%s: Can't read PPM signature", name);
	ascii = (buffer[0] == 'P') && (buffer[1] == '3');
	binary = (buffer[0] == 'P') && (buffer[1] == '6');
	if ((!binary) && (!ascii)) error("%s: Unsupported PPM format", name);
	
	
	/* PPM is always in RGB space */
	img->space = icSigRgbData;
	
	
	/* Read image width and height */
	while (getline(&buffer, &bsize, fp) != -1)
		if (*buffer != '#') break;
	if (sscanf(buffer, "%lu%lu", &(img->w), &(img->h)) != 2)
		error("%s: Error reading width/height of PPM file", name);
	
	img->buffer = NULL;
	img->private = xmalloc(sizeof(ppm_private));
	
	
	/* Read bits/sample value */
	while (getline(&buffer, &bsize, fp) != -1)
		if (*buffer != '#') break;
	if (sscanf(buffer, "%i", &max) != 1)
		error("%s: Error reading maxval of PPM file", name);
	switch (max) {
		case 255:
			bits = 8; break;
		case 65535:
			bits = 16; break;
		default:
			error("%s: Only 8 and 16 bit/sample supported; maxval was %i", name, max);
	}
	
	
	/* Initialize private data */
	((ppm_private *)(img->private))->bits = bits;
	((ppm_private *)(img->private))->r = -2;
	((ppm_private *)(img->private))->imgstart = ftell((FILE *)(img->fd));
	
	
	/* Image methods */
	img->alloc = PPM_IMG_Alloc;
	img->clean = PPM_IMG_FreeRow;
	img->close = PPM_IMG_Close;
	img->get_field = (field_method)PPM_IMG_GetField;
	img->set_field = (field_method)uninitialized;
	img->get_row = ascii ? PPM_IMG_GetRow_ascii : PPM_IMG_GetRow_binary;
	img->set_row = (row_method)uninitialized;
	img->get_pixel = unpack(img->space, bits, PLANARCONFIG_CONTIG);
	img->set_pixel = (set_pixel_method)uninitialized;
	img->print_stats = PPM_IMG_PrintStats;
	
	free(buffer);
	
	return img;
}
