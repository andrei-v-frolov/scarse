/* $Id: imageio.h,v 1.1.1.1 2001/01/26 22:45:32 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Public declarations for image IO wrapper library.
 * 
 * Copyright (C) 1999 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#include <stdio.h>
#include <tiffio.h>


/**********************************************************************/

#ifndef __IMAGEIO_H__
#define __IMAGEIO_H__


/* Image structure */
typedef struct _image {
	char *file; void *fd;		/* Image file name & descriptor */
	unsigned int space;		/* Image color space */
	unsigned long w, h;		/* Width and height */
	void *buffer;			/* Buffer for row data */
	void *private;			/* Pointer to private data */
	
	/* object methods */
	void (*alloc)(struct _image *img);
	void (*clean)(struct _image *img);
	void (*close)(struct _image *img);
	 int (*get_field)(struct _image *img, ttag_t tag, ...);
	 int (*set_field)(struct _image *img, ttag_t tag, ...);
	void (*get_row)(struct _image *img, unsigned long r);
	void (*set_row)(struct _image *img, unsigned long r);
	void (*get_pixel)(struct _image *img, unsigned long x, double v[]);
	 int (*set_pixel)(struct _image *img, unsigned long x, double v[]);
	void (*print_stats)(struct _image *img, FILE *fp);
} image;

/* Method typedefs */
typedef void (*image_method)(image *img);
typedef  int (*field_method)(image *img, ttag_t tag, ...);
typedef void (*row_method)(image *img, unsigned long r);
typedef void (*get_pixel_method)(image *img, unsigned long x, double v[]);
typedef  int (*set_pixel_method)(image *img, unsigned long x, double v[]);
typedef void (*stats_method)(image *img, FILE *fp);

/* Convenience macros for calling object methods */
#define IMG_Alloc(img) (*(img->alloc))(img)
#define IMG_Clean(img) (*(img->clean))(img)
#define IMG_Close(img) (*(img->close))(img)
#define IMG_GetField(img, tag, val) (*(img->get_field))(img, tag, val)
#define IMG_SetField(img, tag, val) (*(img->set_field))(img, tag, val)
#define IMG_ReadRow(img, r) (*(img->get_row))(img, r)
#define IMG_WriteRow(img, r) (*(img->set_row))(img, r)
#define IMG_ReadPixel(img, x, v) (*(img->get_pixel))(img, x, v)
#define IMG_WritePixel(img, x, v) (*(img->set_pixel))(img, x, v)
#define IMG_PrintStats(img, fp) (*(img->print_stats))(img, fp)

/* Pseudo-tag for setting/reading colorspace */
#define TIFFTAG_COLORSPACE 65566

/* Flag for low-level packing primitives */
extern int noisify;

/* Dummy image access method */
void uninitialized(image *img, ...);

/* Pixel access methods for various configurations */
set_pixel_method pack(unsigned int space, short bits, short planarcfg);
get_pixel_method unpack(unsigned int space, short bits, short planarcfg);

/* Open image */
image *IMG_Open(const char *name, const char *mode);
image *PPM_IMG_Open(const char *name, const char *mode);
image *TIFF_IMG_Open(const char *name, const char *mode);

/* Copy image header */
void IMG_CpHeader(image *in, image *out);


#endif /* __IMAGEIO_H__ */
