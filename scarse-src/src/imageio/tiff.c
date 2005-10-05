/* $Id: tiff.c,v 1.2 2005/10/05 06:29:32 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * TIFF image format wrappers.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

/* CREDITS:
 *   TIFF code based on libtiff & tools by Sam Leffler
 */

/* TODO:
 *
 */

#include "../scarse.h"


/* Private data for TIFF image wrappers */
typedef struct { short samples; } tiff_private;



/******************* Reading/writing tags *****************************/

/* Get TIFF image color space */
static icColorSpaceSignature TIFF_IMG_GetColorSpace(image *img)
{
	short photo, samples;
	TIFF *t = (TIFF *)(img->fd);
	
	if (img->space) return img->space;
	
	TIFFGetField(t, TIFFTAG_PHOTOMETRIC, &photo);
	TIFFGetField(t, TIFFTAG_SAMPLESPERPIXEL, &samples);
	
	switch(photo) {
		case PHOTOMETRIC_MINISWHITE:
		case PHOTOMETRIC_MINISBLACK:
			return icSigGrayData;
		case PHOTOMETRIC_RGB:
			return icSigRgbData;
		case PHOTOMETRIC_SEPARATED:
			if (samples == 3) return icSigCmyData;
			if (samples == 4) return icSigCmykData;
		case PHOTOMETRIC_YCBCR:
			return icSigYCbCrData;
		case PHOTOMETRIC_CIELAB:
			return icSigLabData;
		default:
			error("PhotometricIntent [%i] is not supported by ICC", photo);
	}
	
	return 0;
}

/* Set TIFF image color space */
static void TIFF_IMG_SetColorSpace(image *img, icColorSpaceSignature space)
{
	short photo, samples = 3;
	TIFF *t = (TIFF *)(img->fd);
	
	switch(space) {
		case icSigGrayData:
			photo = PHOTOMETRIC_MINISBLACK;
			samples = 1;
			break;
		case icSigRgbData:
			photo = PHOTOMETRIC_RGB;
			break;
		case icSigCmyData:
			photo = PHOTOMETRIC_SEPARATED;
			break;
		case icSigCmykData:
			photo = PHOTOMETRIC_SEPARATED;
			samples = 4;
			break;
		case icSigYCbCrData:
			photo = PHOTOMETRIC_YCBCR;
			break;
		case icSigLabData:
			photo = PHOTOMETRIC_CIELAB;
			break;
		default:
			error("Color space [%s] is not supported by TIFF", ColorSpaceSignature2str(space));
	}
	
	img->space = space;
	
	TIFFSetField(t, TIFFTAG_PHOTOMETRIC, photo);
	TIFFSetField(t, TIFFTAG_SAMPLESPERPIXEL, samples);
}


/* Get TIFF image tag */
static int TIFF_IMG_GetField(image *img, ttag_t tag, ...)
{
	int r = 0;
	va_list args;
	
	va_start(args, tag);
	
	switch (tag) {
		case TIFFTAG_COLORSPACE:
			{
				unsigned int *sp = va_arg(args, unsigned int *);
				
				*sp = TIFF_IMG_GetColorSpace(img); r = 1;
			}
			break;
		default:			/* Pass args to libtiff */
			r = TIFFVGetField((TIFF *)(img->fd), tag, args);
			break;
	}
	
	va_end(args);
	return r;
}

/* Set TIFF image tag */
static int TIFF_IMG_SetField(image *img, ttag_t tag, ...)
{
	int r = 0;
	va_list args;
	
	va_start(args, tag);
	
	switch (tag) {
		case TIFFTAG_COLORSPACE:
			{
				unsigned int sp = va_arg(args, unsigned int);
				
				TIFF_IMG_SetColorSpace(img, sp); r = 1;
			}
			break;
		default:			/* Pass args to libtiff */
			r = TIFFVSetField((TIFF *)(img->fd), tag, args);
			break;
	}
	
	va_end(args);
	return r;
}



/******************* Row access ***************************************/

/* Read contigious row into buffer */
static void TIFF_IMG_GetRow_contigious(image *img, unsigned long r)
{
	if (TIFFReadScanline((TIFF *)(img->fd), img->buffer, r, 0) < 0)
		error("%s: Error reading row %lu", img->file, r);
}

/* Write contigious row from buffer */
static void TIFF_IMG_SetRow_contigious(image *img, unsigned long r)
{
	if (TIFFWriteScanline((TIFF *)(img->fd), img->buffer, r, 0) < 0)
		error("%s: Error writing row %lu", img->file, r);
}

/* Free contigious row buffer */
static void TIFF_IMG_FreeRow_contigious(image *img)
{
	free(img->buffer);
}


/* Read interlaced row into buffer */
static void TIFF_IMG_GetRow_interlaced(image *img, unsigned long r)
{
	int i;
	
	for (i = 0; i < ((tiff_private *)(img->private))->samples; i++)
		if (TIFFReadScanline((TIFF *)(img->fd), ((void **)(img->buffer))[i], r, i) < 0)
			error("%s: Error reading row %lu", img->file, r);
}

/* Write interlaced row from buffer */
static void TIFF_IMG_SetRow_interlaced(image *img, unsigned long r)
{
	int i;
	
	for (i = 0; i < ((tiff_private *)(img->private))->samples; i++)
		if (TIFFWriteScanline((TIFF *)(img->fd), ((void **)(img->buffer))[i], r, i) < 0)
			error("%s: Error writing row %lu", img->file, r);
}

/* Free interlaced row buffer */
static void TIFF_IMG_FreeRow_interlaced(image *img)
{
	int i;
	
	for (i = 0; i < ((tiff_private *)(img->private))->samples; i++)
		free(((void **)(img->buffer))[i]);
	
	free(img->buffer);
}



/******************* Image stats **************************************/

/* Print TIFF image stats */
static void TIFF_IMG_PrintStats(image *img, FILE *fp)
{
	fprintf(fp, "\n%s: ", img->file);
	TIFFPrintDirectory((TIFF *)(img->fd), fp, TIFFPRINT_NONE);
}



/******************* TIFF image input/output **************************/

/* Allocate buffers and initialize access methods */
static void TIFF_IMG_Alloc(image *img)
{
	int i;
	uint32 w, h, strip;
	unsigned int space;
	short samples, bits, planarcfg;
	TIFF *t = (TIFF *)(img->fd);
	
	/* Read required tags */
	IMG_GetField(img, TIFFTAG_IMAGEWIDTH, &w);
	IMG_GetField(img, TIFFTAG_IMAGELENGTH, &h);
	IMG_GetField(img, TIFFTAG_COLORSPACE, &space);
	IMG_GetField(img, TIFFTAG_SAMPLESPERPIXEL, &samples);
	
	/* Read optional tags, set to default values if undefined */
	if (!IMG_GetField(img, TIFFTAG_BITSPERSAMPLE, &bits))
		IMG_SetField(img, TIFFTAG_BITSPERSAMPLE, bits=8);
	if (!IMG_GetField(img, TIFFTAG_PLANARCONFIG, &planarcfg))
		IMG_SetField(img, TIFFTAG_PLANARCONFIG, planarcfg=PLANARCONFIG_CONTIG);
	
	if (!IMG_GetField(img, TIFFTAG_ROWSPERSTRIP, &strip))
		IMG_SetField(img, TIFFTAG_ROWSPERSTRIP,
			planarcfg == PLANARCONFIG_CONTIG ? (2<<19)/(w*samples*bits)+1 : 1);
	else if (strip != 1 && planarcfg != PLANARCONFIG_CONTIG) {
		warning("%s: interlaced configuration requested, forcing rows per strip to 1", img->file);
		IMG_SetField(img, TIFFTAG_ROWSPERSTRIP, 1);
	}
	
	/* Initialize the rest of image structure fields */
	img->space = space;
	img->w = w; img->h = h;
	((tiff_private *)(img->private))->samples = samples;
	
	switch (planarcfg) {
		case PLANARCONFIG_CONTIG:
			img->buffer = xmalloc(TIFFScanlineSize(t));
			
			img->clean = TIFF_IMG_FreeRow_contigious;
			img->get_row = TIFF_IMG_GetRow_contigious;
			img->set_row = TIFF_IMG_SetRow_contigious;
			break;
		case PLANARCONFIG_SEPARATE:
			img->buffer = xmalloc(samples*sizeof(void *));
			for (i = 0; i < samples; i++)
				((void **)(img->buffer))[i] = xmalloc(TIFFScanlineSize(t));
			
			img->clean = TIFF_IMG_FreeRow_interlaced;
			img->get_row = TIFF_IMG_GetRow_interlaced;
			img->set_row = TIFF_IMG_SetRow_interlaced;
			break;
		default:
			error("%s: Unsupported PlanarConfig (%i)", img->file, planarcfg);
	}
	
	img->get_pixel = unpack(space, bits, planarcfg);
	img->set_pixel = pack(space, bits, planarcfg);
}

/* Close TIFF image */
static void TIFF_IMG_Close(image *img)
{
	IMG_Clean(img);
	
	free(img->file);
	TIFFClose((TIFF *)(img->fd));
	free(img->private);
	
	free(img);
}

/* Open TIFF image */
image *TIFF_IMG_Open(const char *name, const char *mode)
{
	image *img = (image *)xmalloc(sizeof(image));
	
	img->file = xstrdup(name);
	
	img->fd = (void *)TIFFOpen(name, mode);
	if (!img->fd) error("Can't open file '%s'", name);
	
	img->space = 0;
	img->w = img->h = 0;
	img->buffer = NULL;
	img->private = xmalloc(sizeof(tiff_private));
	
	img->alloc = TIFF_IMG_Alloc;
	img->clean = (image_method)uninitialized;
	img->close = TIFF_IMG_Close;
	img->get_field = TIFF_IMG_GetField;
	img->set_field = TIFF_IMG_SetField;
	img->get_row = img->set_row = (row_method)uninitialized;
	img->get_pixel = (get_pixel_method)uninitialized;
	img->set_pixel = (set_pixel_method)uninitialized;
	img->print_stats = TIFF_IMG_PrintStats;
	
	return img;
}
