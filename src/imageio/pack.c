/* $Id: pack.c,v 1.4 2005/10/20 06:15:04 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Primitives for packing/unpacking FP pixel values.
 * 
 * Copyright (C) 1999-2005 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <frolov@cita.utoronto.ca>
 * 
 */

#include "../scarse.h"


/******************* Clipping and quantizing **************************/

/* Use system random numbers - not like we are doing cryptography here! */
static double ran()
{
	register double x = rand()/(RAND_MAX+1.0);
	
	return x;
}

/* Quantization sub-bit noise flag */
int noisify = 0;

/* Quantize value to specified number of bits */
static int quantize(int bits, double v)
{
	register double max = (double)((1 << bits) - 1);
	register double sv = max*v + (noisify ? ran() : 0.5);
	
	return (int)floor(sv);
}

/* Clip value to specified range */
static double clip(double v, double min, double max, int *clipped)
{
	if (v < min) { *clipped |= 0x01; return min; }
	if (v > max) { *clipped |= 0x01; return max; }
	
	return v;
}



/******************* Unitialized method placeholders ******************/

/* Dummy image access method */
void uninitialized(image *img, ...)
{
	error("Internal error: uninitialized method called while accessing %s", img->file);
}



/******************* Grayscale images *********************************/

/* Pack gray pixel into 8-bit buffer */
static int pack1_8(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned char *p = &(((unsigned char *)(img->buffer))[x]);
	
	*p = quantize(8, clip(*v, 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack gray pixel from 8-bit buffer */
static void unpack1_8(image *img, unsigned long x, double v[])
{
	unsigned char *p = &(((unsigned char *)(img->buffer))[x]);
	
	*v = (double)(*p)/255.0;
}


/* Pack gray pixel into 16-bit buffer */
static int pack1_16(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned short *p = &(((unsigned short *)(img->buffer))[x]);
	
	*p = quantize(16, clip(*v, 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack gray pixel from 16-bit buffer */
static void unpack1_16(image *img, unsigned long x, double v[])
{
	unsigned short *p = &(((unsigned short *)(img->buffer))[x]);
	
	*v = (double)(*p)/65535.0;
}



/******************* RGB, XYZ and related spaces **********************/

/* Pack RGB pixel into 8-bit contigious buffer */
static int pack3_8c(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned char *p = &(((unsigned char *)(img->buffer))[3*x]);
	
	p[0] = quantize(8, clip(v[0], 0.0, 1.0, &rv));
	p[1] = quantize(8, clip(v[1], 0.0, 1.0, &rv));
	p[2] = quantize(8, clip(v[2], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack RGB pixel from 8-bit contigious buffer */
static void unpack3_8c(image *img, unsigned long x, double v[])
{
	unsigned char *p = &(((unsigned char *)(img->buffer))[3*x]);
	
	v[0] = (double)(p[0])/255.0;
	v[1] = (double)(p[1])/255.0;
	v[2] = (double)(p[2])/255.0;
}


/* Pack RGB pixel into 8-bit interlaced buffer */
static int pack3_8i(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned char **p = (unsigned char **)(img->buffer);
	
	p[0][x] = quantize(8, clip(v[0], 0.0, 1.0, &rv));
	p[1][x] = quantize(8, clip(v[1], 0.0, 1.0, &rv));
	p[2][x] = quantize(8, clip(v[2], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack RGB pixel from 8-bit interlaced buffer */
static void unpack3_8i(image *img, unsigned long x, double v[])
{
	unsigned char **p = (unsigned char **)(img->buffer);
	
	v[0] = (double)(p[0][x])/255.0;
	v[1] = (double)(p[1][x])/255.0;
	v[2] = (double)(p[2][x])/255.0;
}


/* Pack RGB pixel into 16-bit contigious buffer */
static int pack3_16c(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned short *p = &(((unsigned short *)(img->buffer))[3*x]);
	
	p[0] = quantize(16, clip(v[0], 0.0, 1.0, &rv));
	p[1] = quantize(16, clip(v[1], 0.0, 1.0, &rv));
	p[2] = quantize(16, clip(v[2], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack RGB pixel from 16-bit contigious buffer */
static void unpack3_16c(image *img, unsigned long x, double v[])
{
	unsigned short *p = &(((unsigned short *)(img->buffer))[3*x]);
	
	v[0] = (double)(p[0])/65535.0;
	v[1] = (double)(p[1])/65535.0;
	v[2] = (double)(p[2])/65535.0;
}


/* Pack RGB pixel into 16-bit interlaced buffer */
static int pack3_16i(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned short **p = (unsigned short **)(img->buffer);
	
	p[0][x] = quantize(16, clip(v[0], 0.0, 1.0, &rv));
	p[1][x] = quantize(16, clip(v[1], 0.0, 1.0, &rv));
	p[2][x] = quantize(16, clip(v[2], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack RGB pixel from 16-bit interlaced buffer */
static void unpack3_16i(image *img, unsigned long x, double v[])
{
	unsigned short **p = (unsigned short **)(img->buffer);
	
	v[0] = (double)(p[0][x])/65535.0;
	v[1] = (double)(p[1][x])/65535.0;
	v[2] = (double)(p[2][x])/65535.0;
}



/******************* CMYK and related spaces **************************/

/* Pack CMYK pixel into 8-bit contigious buffer */
static int pack4_8c(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned char *p = &(((unsigned char *)(img->buffer))[4*x]);
	
	p[0] = quantize(8, clip(v[0], 0.0, 1.0, &rv));
	p[1] = quantize(8, clip(v[1], 0.0, 1.0, &rv));
	p[2] = quantize(8, clip(v[2], 0.0, 1.0, &rv));
	p[3] = quantize(8, clip(v[3], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack CMYK pixel from 8-bit contigious buffer */
static void unpack4_8c(image *img, unsigned long x, double v[])
{
	unsigned char *p = &(((unsigned char *)(img->buffer))[4*x]);
	
	v[0] = (double)(p[0])/255.0;
	v[1] = (double)(p[1])/255.0;
	v[2] = (double)(p[2])/255.0;
	v[3] = (double)(p[3])/255.0;
}


/* Pack CMYK pixel into 8-bit interlaced buffer */
static int pack4_8i(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned char **p = (unsigned char **)(img->buffer);
	
	p[0][x] = quantize(8, clip(v[0], 0.0, 1.0, &rv));
	p[1][x] = quantize(8, clip(v[1], 0.0, 1.0, &rv));
	p[2][x] = quantize(8, clip(v[2], 0.0, 1.0, &rv));
	p[3][x] = quantize(8, clip(v[3], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack CMYK pixel from 8-bit interlaced buffer */
static void unpack4_8i(image *img, unsigned long x, double v[])
{
	unsigned char **p = (unsigned char **)(img->buffer);
	
	v[0] = (double)(p[0][x])/255.0;
	v[1] = (double)(p[1][x])/255.0;
	v[2] = (double)(p[2][x])/255.0;
	v[3] = (double)(p[3][x])/255.0;
}


/* Pack CMYK pixel into 16-bit contigious buffer */
static int pack4_16c(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned short *p = &(((unsigned short *)(img->buffer))[4*x]);
	
	p[0] = quantize(16, clip(v[0], 0.0, 1.0, &rv));
	p[1] = quantize(16, clip(v[1], 0.0, 1.0, &rv));
	p[2] = quantize(16, clip(v[2], 0.0, 1.0, &rv));
	p[3] = quantize(16, clip(v[3], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack CMYK pixel from 16-bit contigious buffer */
static void unpack4_16c(image *img, unsigned long x, double v[])
{
	unsigned short *p = &(((unsigned short *)(img->buffer))[4*x]);
	
	v[0] = (double)(p[0])/65535.0;
	v[1] = (double)(p[1])/65535.0;
	v[2] = (double)(p[2])/65535.0;
	v[3] = (double)(p[3])/65535.0;
}


/* Pack CMYK pixel into 16-bit interlaced buffer */
static int pack4_16i(image *img, unsigned long x, double v[])
{
	int rv = 0;
	unsigned short **p = (unsigned short **)(img->buffer);
	
	p[0][x] = quantize(16, clip(v[0], 0.0, 1.0, &rv));
	p[1][x] = quantize(16, clip(v[1], 0.0, 1.0, &rv));
	p[2][x] = quantize(16, clip(v[2], 0.0, 1.0, &rv));
	p[3][x] = quantize(16, clip(v[3], 0.0, 1.0, &rv));
	
	return rv;
}

/* Unpack CMYK pixel from 16-bit interlaced buffer */
static void unpack4_16i(image *img, unsigned long x, double v[])
{
	unsigned short **p = (unsigned short **)(img->buffer);
	
	v[0] = (double)(p[0][x])/65535.0;
	v[1] = (double)(p[1][x])/65535.0;
	v[2] = (double)(p[2][x])/65535.0;
	v[3] = (double)(p[3][x])/65535.0;
}



/******************* CIE Lab packing **********************************/

/* Pack Lab pixel into 8-bit contigious buffer */
static int packL_8c(image *img, unsigned long x, double v[])
{
	int rv = 0;
	signed char *p = &(((signed char *)(img->buffer))[3*x]);
	
	*((unsigned char *)p) = quantize(8, clip(v[0]/100.0, 0.0, 1.0, &rv));
	p[1] = quantize(8, clip(v[1]/256.0, -0.5, 0.5, &rv));
	p[2] = quantize(8, clip(v[2]/256.0, -0.5, 0.5, &rv));
	
	return rv;
}

/* Unpack Lab pixel from 8-bit contigious buffer */
static void unpackL_8c(image *img, unsigned long x, double v[])
{
	signed char *p = &(((signed char *)(img->buffer))[3*x]);
	
	v[0] = (double)(*((unsigned char *)p))/2.55;
	v[1] = (double)(p[1]);
	v[2] = (double)(p[2]);
}


/* Pack Lab pixel into 8-bit interlaced buffer */
static int packL_8i(image *img, unsigned long x, double v[])
{
	int rv = 0;
	signed char **p = (signed char **)(img->buffer);
	
	((unsigned char *)p[0])[x] = quantize(8, clip(v[0]/100.0, 0.0, 1.0, &rv));
	p[1][x] = quantize(8, clip(v[1]/256.0, -0.5, 0.5, &rv));
	p[2][x] = quantize(8, clip(v[2]/256.0, -0.5, 0.5, &rv));
	
	return rv;
}

/* Unpack Lab pixel from 8-bit interlaced buffer */
static void unpackL_8i(image *img, unsigned long x, double v[])
{
	signed char **p = (signed char **)(img->buffer);
	
	v[0] = (double)(((unsigned char *)p[0])[x])/2.55;
	v[1] = (double)(p[1][x]);
	v[2] = (double)(p[2][x]);
}


/* Pack Lab pixel into 16-bit contigious buffer */
static int packL_16c(image *img, unsigned long x, double v[])
{
	int rv = 0;
	signed short *p = &(((signed short *)(img->buffer))[3*x]);
	
	*((unsigned short *)p) = quantize(16, clip(v[0]/100.0, 0.0, 1.0, &rv));
	p[1] = quantize(16, clip(v[1]/256.0, -0.5, 0.5, &rv));
	p[2] = quantize(16, clip(v[2]/256.0, -0.5, 0.5, &rv));
	
	return rv;
}

/* Unpack Lab pixel from 16-bit contigious buffer */
static void unpackL_16c(image *img, unsigned long x, double v[])
{
	signed short *p = &(((signed short *)(img->buffer))[3*x]);
	
	v[0] = (double)(*((unsigned short *)p))/655.35;
	v[1] = (double)(p[1])/256.0;
	v[2] = (double)(p[2])/256.0;
}


/* Pack Lab pixel into 16-bit interlaced buffer */
static int packL_16i(image *img, unsigned long x, double v[])
{
	int rv = 0;
	signed short **p = (signed short **)(img->buffer);
	
	((unsigned short *)p[0])[x] = quantize(16, clip(v[0]/100.0, 0.0, 1.0, &rv));
	p[1][x] = quantize(16, clip(v[1]/256.0, -0.5, 0.5, &rv));
	p[2][x] = quantize(16, clip(v[2]/256.0, -0.5, 0.5, &rv));
	
	return rv;
}

/* Unpack Lab pixel from 16-bit interlaced buffer */
static void unpackL_16i(image *img, unsigned long x, double v[])
{
	signed short **p = (signed short **)(img->buffer);
	
	v[0] = (double)(((unsigned short *)p[0])[x])/655.35;
	v[1] = (double)(p[1][x])/256.0;
	v[2] = (double)(p[2][x])/256.0;
}



/******************* Select packing methods ***************************/

/* Return set_pixel method */
set_pixel_method pack(unsigned int space, short bits, short planarcfg)
{
	if (bits != 8 && bits != 16)
		error("Unsupported packing configuration: %i bits/channel", bits);
	
	if (planarcfg != PLANARCONFIG_CONTIG && planarcfg != PLANARCONFIG_SEPARATE)
		error("Unsupported packing configuration: %i PlanarConfig", planarcfg);
	
	switch (space) {
		case icSigGrayData:
			return (bits == 8) ? pack1_8 : pack1_16;
		case icSigRgbData:
		case icSigCmyData:
		case icSigXYZData:
		case icSigYxyData:
			return (planarcfg == PLANARCONFIG_CONTIG) ?
				((bits == 8) ? pack3_8c : pack3_16c) :
				((bits == 8) ? pack3_8i : pack3_16i);
		case icSigCmykData:
			return (planarcfg == PLANARCONFIG_CONTIG) ?
				((bits == 8) ? pack4_8c : pack4_16c) :
				((bits == 8) ? pack4_8i : pack4_16i);
		case icSigLabData:
			return (planarcfg == PLANARCONFIG_CONTIG) ?
				((bits == 8) ? packL_8c : packL_16c) :
				((bits == 8) ? packL_8i : packL_16i);
		default:
			error("[%s]: Color space not supported",
				ColorSpaceSignature2str(space));
	}
	
	return (set_pixel_method)uninitialized;
}

/* Return get_pixel method */
get_pixel_method unpack(unsigned int space, short bits, short planarcfg)
{
	if (bits != 8 && bits != 16)
		error("Unsupported packing configuration: %i bits/channel", bits);
	
	if (planarcfg != PLANARCONFIG_CONTIG && planarcfg != PLANARCONFIG_SEPARATE)
		error("Unsupported packing configuration: %i PlanarConfig", planarcfg);
	
	switch (space) {
		case icSigGrayData:
			return (bits == 8) ? unpack1_8 : unpack1_16;
		case icSigRgbData:
		case icSigCmyData:
		case icSigXYZData:
		case icSigYxyData:
			return (planarcfg == PLANARCONFIG_CONTIG) ?
				((bits == 8) ? unpack3_8c : unpack3_16c) :
				((bits == 8) ? unpack3_8i : unpack3_16i);
		case icSigCmykData:
			return (planarcfg == PLANARCONFIG_CONTIG) ?
				((bits == 8) ? unpack4_8c : unpack4_16c) :
				((bits == 8) ? unpack4_8i : unpack4_16i);
		case icSigLabData:
			return (planarcfg == PLANARCONFIG_CONTIG) ?
				((bits == 8) ? unpackL_8c : unpackL_16c) :
				((bits == 8) ? unpackL_8i : unpackL_16i);
		default:
			error("[%s]: Color space not supported",
				ColorSpaceSignature2str(space));
	}
	
	return (get_pixel_method)uninitialized;
}
