/* $Id: imageio.c,v 1.1.1.1 2001/01/26 22:45:32 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * High-level image IO routines.
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

#define _GNU_SOURCE

#include <fnmatch.h>
#include "imageio.h"
#include "../util.h"


/**********************************************************************/

/* Open image */
image *IMG_Open(const char *name, const char *mode)
{
	if (!fnmatch("*.tif", name, FNM_CASEFOLD)) return TIFF_IMG_Open(name, mode);
	if (!fnmatch("*.p[np]m", name, FNM_CASEFOLD)) return PPM_IMG_Open(name, mode);
	
	warning("%s: Confused by file extension, will open as TIFF", name);
	return TIFF_IMG_Open(name, mode);
}


/* Copy image header */
void IMG_CpHeader(image *in, image *out)
{
	unsigned int space;
	short shortv, shortv2;
	short *shortav;
	uint32 longv;
	float floatv;
	char *stringv;
	
	#define	CopyField(tag, v) \
		if ((*in->get_field)(in, tag, &v)) (*out->set_field)(out, tag, v)
	#define	CopyField2(tag, v1, v2) \
		if ((*in->get_field)(in, tag, &v1, &v2)) (*out->set_field)(out, tag, v1, v2)
	
	CopyField(TIFFTAG_COLORSPACE, space);
	
	CopyField(TIFFTAG_SUBFILETYPE, longv);
	CopyField(TIFFTAG_IMAGEWIDTH, longv);
	CopyField(TIFFTAG_IMAGELENGTH, longv);
	CopyField(TIFFTAG_BITSPERSAMPLE, shortv);
	CopyField(TIFFTAG_COMPRESSION, shortv);
	CopyField(TIFFTAG_PREDICTOR, shortv);
	CopyField(TIFFTAG_PHOTOMETRIC, shortv);
	CopyField(TIFFTAG_ORIENTATION, shortv);
	CopyField(TIFFTAG_SAMPLESPERPIXEL, shortv);
	CopyField(TIFFTAG_XRESOLUTION, floatv);
	CopyField(TIFFTAG_YRESOLUTION, floatv);
	CopyField(TIFFTAG_RESOLUTIONUNIT, shortv);
	CopyField(TIFFTAG_PLANARCONFIG, shortv);
	CopyField(TIFFTAG_ROWSPERSTRIP, longv);
	CopyField(TIFFTAG_XPOSITION, floatv);
	CopyField(TIFFTAG_YPOSITION, floatv);
	CopyField2(TIFFTAG_EXTRASAMPLES, shortv, shortav);
	CopyField2(TIFFTAG_PAGENUMBER, shortv, shortv2);
	CopyField(TIFFTAG_ARTIST, stringv);
	CopyField(TIFFTAG_IMAGEDESCRIPTION, stringv);
	CopyField(TIFFTAG_MAKE, stringv);
	CopyField(TIFFTAG_MODEL, stringv);
	CopyField(TIFFTAG_SOFTWARE, stringv);
	CopyField(TIFFTAG_DATETIME, stringv);
	CopyField(TIFFTAG_HOSTCOMPUTER, stringv);
	CopyField(TIFFTAG_PAGENAME, stringv);
	CopyField(TIFFTAG_DOCUMENTNAME, stringv);
	
	#undef CopyField
	#undef CopyField2
}

