/* $Id: spaces.c,v 1.5 2001/06/27 03:58:57 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Generic color space conversions.
 * 
 * Copyright (C) 1999-2001 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

/* CREDITS:
 *   Based on Color Spaces FAQ by David Bourgin
 *   CIE Lab code borrowed from icclib & examples by Graeme W. Gill
 *   HSV code adopted from ???
 */

/* TODO:
 *
 */

#include <string.h>
#include <math.h>
#include <icc.h>

#include "spaces.h"
#include "util.h"



/******************* CIE XYZ based spaces *****************************/

/* IMPORTANT NOTE:
 * ===============
 *   XYZ values used here are based on RELATIVE colorimetry, meaning they
 *   are scaled so that physical white point is always mapped to the same
 *   value - the so-called PCS illuminant (specified to be D50 by ICC specs).
 *   In addition, the max Y value we use is 1.0, not 100.0 as it's often set!
 */

/* PCS illuminant and white point (in absolute XYZ representation) */
double XYZ_ILLUM[3] = {0.9642, 1.0000, 0.8249};	/* Always CIE illuminant D50 */
double   XYZ_WPT[3] = {0.9505, 1.0000, 1.0891};	/* Adobe RGB (D65) by default */


/* CIE (relative) XYZ to CIE (absolute) Yxy */
void XYZ2Yxy(double in[], double out[])
{
	double *I = XYZ_ILLUM, *W = XYZ_WPT;
	double X = (in[0]/I[0])*W[0];
	double Y = (in[1]/I[1])*W[1];
	double Z = (in[2]/I[2])*W[2];
	double S = X + Y + Z;
	
	out[0] = Y;
	out[1] = X/S;
	out[2] = Y/S;
}

/* CIE (absolute) Yxy to CIE (relative) XYZ */
void Yxy2XYZ(double in[], double out[])
{
	double *I = XYZ_ILLUM, *W = XYZ_WPT;
	double Y = in[0], x = in[1], y = in[2];
	double z = 1.0 - (x + y), X = (x/y)*Y, Z = (z/y)*Y;
	
	out[0] = (X/W[0])*I[0];
	out[1] = (Y/W[1])*I[1];
	out[2] = (Z/W[2])*I[2];
}


/* CIE (relative) XYZ to perceptual Lab */
void XYZ2Lab(double in[], double out[])
{
	double *I = XYZ_ILLUM, fx, fy, fz, L;
	double x = in[0]/I[0], y = in[1]/I[1], z = in[2]/I[2];
	
	if (x > 0.008856451586)
		fx = pow(x, 1.0/3.0);
	else
		fx = 7.787036979 * x + 16.0/116.0;
	
	if (y > 0.008856451586) {
		fy = pow(y, 1.0/3.0);
		L = 116.0 * fy - 16.0;
	} else {
		fy = 7.787036979 * y + 16.0/116.0;
		L = 903.2963058 * y;
	}
	
	if (z > 0.008856451586)
		fz = pow(z, 1.0/3.0);
	else
		fz = 7.787036979 * z + 16.0/116.0;
	
	out[0] = L;
	out[1] = 500.0 * (fx - fy);
	out[2] = 200.0 * (fy - fz);
}

/* Perceptual Lab to CIE (relative) XYZ */
void Lab2XYZ(double in[], double out[])
{
	double L = in[0], a = in[1], b = in[2];
	double *I = XYZ_ILLUM, x, y, z, fx, fy, fz;
	
	if (L > 8.0) {
		fy = (L + 16.0)/116.0;
		y = pow(fy, 3.0);
	} else {
		y = L/903.2963058;
		fy = 7.787036979 * y + 16.0/116.0;
	}
	
	fx = a/500.0 + fy;
	if (fx > 24.0/116.0)
		x = pow(fx, 3.0);
	else
		x = (fx - 16.0/116.0)/7.787036979;
	
	fz = fy - b/200.0;
	if (fz > 24.0/116.0)
		z = pow(fz, 3.0);
	else
		z = (fz - 16.0/116.0)/7.787036979;
	
	out[0] = x * I[0];
	out[1] = y * I[1];
	out[2] = z * I[2];
}


/* CIE (relative) XYZ to Luv */
void XYZ2Luv(double in[], double out[])
{
	double X = in[0], Y = in[1], Z = in[2];
	double *I = XYZ_ILLUM, y = Y/I[1], L, u, v, un, vn;
	
	if (y > 0.008856451586)
		L = 116.0 * pow(y, 1.0/3.0) - 16.0;
	else
		L = 903.2963058 * y;
	
	u = 4.0 * X/(X + 15.0*Y + 3.0*Z);
	v = 9.0 * Y/(X + 15.0*Y + 3.0*Z);
	
	un = 4.0 * I[0]/(I[0] + 15.0*I[1] + 3.0*I[2]);
	vn = 9.0 * I[1]/(I[0] + 15.0*I[1] + 3.0*I[2]);
	
	out[0] = L;
	out[1] = 13.0 * L * (u - un);
	out[2] = 13.0 * L * (v - vn);
}

/* Luv to CIE (relative) XYZ */
void Luv2XYZ(double in[], double out[])
{
	double L = in[0], u = in[1], v = in[2];
	double *I = XYZ_ILLUM, y, un, vn, X, Y, Z;
	
	if (L > 8.0)
		y = pow((L + 16.0)/116.0, 3.0);
	else
		y = L/903.2963058;
	
	un = 4.0 * I[0]/(I[0] + 15.0*I[1] + 3.0*I[2]);
	vn = 9.0 * I[1]/(I[0] + 15.0*I[1] + 3.0*I[2]);
	
	u = u / 13.0 / L + un;
	v = v / 13.0 / L + vn;
	
	Y = y * I[1];
	X = 9.0/4.0 * Y * u/v;
	Z = 3.0*Y/v - 5.0*Y - X/3.0;
	
	out[0] = X;
	out[1] = Y;
	out[2] = Z;
}


/* Gray to CIE (relative) XYZ */
void Gray2XYZ(double *in, double out[])
{
	double g = *in, *I = XYZ_ILLUM;
	
	out[0] = g * I[0];
	out[1] = g * I[1];
	out[2] = g * I[2];
}

/* CIE (relative) XYZ to gray */
void XYZ2Gray(double in[], double *out)
{
	double *I = XYZ_ILLUM;
	
	*out = in[1]/I[1];
}



/***************** Standard RGB primaries and spaces ******************/

/* Standard illuminants (in xy representation) */
static const double D50[2] = {0.3457, 0.3585};
static const double D55[2] = {0.3324, 0.3474};
static const double D65[2] = {0.3127, 0.3290};
static const double D75[2] = {0.2990, 0.3149};
static const double D93[2] = {0.2848, 0.2932};
static const double IlA[2] = {0.4476, 0.4074};
static const double IlB[2] = {0.3484, 0.3516};
static const double IlC[2] = {0.3101, 0.3162};
static const double IlE[2] = {0.3333, 0.3333};

/* Standard RGB primaries (in xy representation) */
static const double WideG[6] = {
	0.7347, 0.2653,
	0.1152, 0.8264,
	0.1566, 0.0177
};
static const double Adobe[6] = {
	0.6400, 0.3300,
	0.2100, 0.7100,
	0.1500, 0.0600
};
static const double Bruce[6] = {
	0.6400, 0.3300,
	0.2800, 0.6500,
	0.1500, 0.0600
};
static const double CIE[6] = {
	0.7350, 0.2650,
	0.2740, 0.7170,
	0.1670, 0.0090
};
static const double EBU[6] = {
	0.6400, 0.3300,
	0.2900, 0.6000,
	0.1500, 0.0600
};
static const double HDTV[6] = {
	0.6400, 0.3300,
	0.3000, 0.6000,
	0.1500, 0.0600
};
static const double NTSC[6] = {
	0.6700, 0.3300,
	0.2100, 0.7100,
	0.1400, 0.0800
};
static const double P22[6] = {
	0.6300, 0.3400,
	0.2950, 0.6050,
	0.1550, 0.0770
};
static const double SMPTE[6] = {
	0.6300, 0.3400,
	0.3100, 0.5950,
	0.1550, 0.0700
};
static const double Trini[6] = {
	0.6250, 0.3400,
	0.2800, 0.5950,
	0.1550, 0.0700
};

/* Index of standard illuminants, primaries and color spaces */
static struct { char *label; const double *wpt, *rgb, g; } primaries_idx[] = {
	/* Standard illuminants */
	{"D50",			D50,	NULL,	0.0},
	{"D55",			D55,	NULL,	0.0},
	{"D65",			D65,	NULL,	0.0},
	{"D75",			D75,	NULL,	0.0},
	{"D93",			D93,	NULL,	0.0},
	{"A",			IlA,	NULL,	0.0},
	{"B",			IlB,	NULL,	0.0},
	{"C",			IlC,	NULL,	0.0},
	{"E",			IlE,	NULL,	0.0},
	
	/* Standard RGB primaries */
	{"700/525/450nm",	NULL,	WideG,	0.0},
	{"EBU",			NULL,	EBU,	0.0},
	{"HDTV",		NULL,	HDTV,	0.0},
	{"P22",			NULL,	P22,	0.0},
	{"Trinitron",		NULL,	Trini,	0.0},
	
	/* Standard RGB spaces */
	{"Adobe",		D65,	Adobe,	2.2},
	{"Apple",		D65,	Trini,	1.8},
	{"Bruce",		D65,	Bruce,	2.2},
	{"ColorMatch",		D50,	P22,	1.8},
	{"CIE",			IlE,	CIE,	2.2},
	{"NTSC",		IlC,	NTSC,	2.2},
	{"PAL/SECAM",		D65,	EBU,	2.2},
	{"sRGB",		D65,	HDTV,	2.2},
	{"SMPTE-C",		D65,	SMPTE,	2.2},
	{"WideGamut",		D50,	WideG,	2.2},
	NULL
};


/* Copy white point and RGB primaries from array */
int LookupPrimaries(char *p, double dest[4][2], double *gamma)
{
	int i;
	
	for (i = 0; primaries_idx[i].label; i++) {
		if (strcasecmp(p, primaries_idx[i].label)) continue;
		
		if (primaries_idx[i].wpt) {	/* copy white point */
			dest[0][0] = primaries_idx[i].wpt[0];
			dest[0][1] = primaries_idx[i].wpt[1];
		}
		
		if (primaries_idx[i].rgb) {	/* copy RGB primaries */
			dest[1][0] = primaries_idx[i].rgb[0];
			dest[1][1] = primaries_idx[i].rgb[1];
			dest[2][0] = primaries_idx[i].rgb[2];
			dest[2][1] = primaries_idx[i].rgb[3];
			dest[3][0] = primaries_idx[i].rgb[4];
			dest[3][1] = primaries_idx[i].rgb[5];
		}
		
		if (primaries_idx[i].g > 0.0) {	/* default RGB gamma */
			if (gamma) *gamma = primaries_idx[i].g;
		}
		
		return 1;
	}
	
	/* No match found */
	return 0;
}


/***************** RGB based device dependant spaces ******************/

/* RGB <=> XYZ conversion matrices */
double M_RGB2XYZ[3][3] = {		/* Adobe RGB by default */
	{0.5850, 0.1882, 0.1910},
	{0.2973, 0.6274, 0.0753},
	{0.0205, 0.0535, 0.7509}
};

double M_XYZ2RGB[3][3] = {		/* Adobe RGB by default */
	{ 2.0125, -0.5650, -0.4551},
	{-0.9554,  1.8760,  0.0549},
	{ 0.0133, -0.1184,  1.3403}
};


/* Set white point and RGB primaries */
void SetPrimaries(double xy[4][2])
{
	double E1[3][3], Y[3];
	
	/* white point (in absolute XYZ representation) */
	double W[3] = {
		        xy[0][0]/xy[0][1],
		               1.0,
		 (1.0-xy[0][0]-xy[0][1])/xy[0][1]
	};
	
	/* RGB primaries matrix (in absolute XYZ representation) */
	double E[3][3] = {
		{       xy[1][0]/xy[1][1],                xy[2][0]/xy[2][1],                xy[3][0]/xy[3][1]        },
		{              1.0,                              1.0,                              1.0               },
		{(1.0-xy[1][0]-xy[1][1])/xy[1][1], (1.0-xy[2][0]-xy[2][1])/xy[2][1], (1.0-xy[3][0]-xy[3][1])/xy[3][1]}
	};
	
	/* relative XYZ scaling factor */
	double *I = XYZ_ILLUM, S[3] = { I[0]/W[0], I[1]/W[1], I[2]/W[2] };
	
	
	inv33(E, E1); apply33(E1, W, Y);
	biscale33(S, E, Y, M_RGB2XYZ);
	inv33(M_RGB2XYZ, M_XYZ2RGB);
	
	XYZ_WPT[0] = W[0];
	XYZ_WPT[1] = W[1];
	XYZ_WPT[2] = W[2];
}


/* RGB with specified primaries to CIE XYZ */
void RGB2XYZ(double in[], double out[])
{
	double t[3] = {in[0], in[1], in[2]};
	
	apply33(M_RGB2XYZ, t, out);
}

/* CIE XYZ to RGB with specified primaries */
void XYZ2RGB(double in[], double out[])
{
	double t[3] = {in[0], in[1], in[2]};
	
	apply33(M_XYZ2RGB, t, out);
}


/* Gray to RGB */
void Gray2RGB(double *in, double out[])
{
	out[0] = out[1] = out[2] = *in;
}


/* RGB to HSV */
/* h = [0,360], s = [0,1], v = [0,1] */
void RGB2HSV(double in[], double out[])
{
	double r = in[0], g = in[1], b = in[2];
	double min, max, delta, h, s, v;
	
	min = r;  if (g < min) min = g;  if (b < min) min = b;
	max = r;  if (g > max) max = g;  if (b > max) max = b;
	
	delta = max - min;
	
	if (max == 0.0) {
		out[0] = out[1] = out[2] = 0.0; return;
	}
	
	v = max;
	s = delta/max;
		
	if (r == max)
		h = (g - b)/delta;		// between yellow & magenta
	else if (g == max)
		h = 2.0 + (b - r)/delta;	// between cyan & yellow
	else
		h = 4.0 + (r - g)/delta;	// between magenta & cyan
		
	h *= 60.0; if (h < 0.0) h += 360.0;
	
	out[0] = h;
	out[1] = s;
	out[2] = v;
}

/* HSV to RGB */
void HSV2RGB(double in[], double out[])
{
	double h = in[0], s = in[1], v = in[2];
	double f, p, q, t, r, g, b;
	int i;
	
	if (s == 0) {
		out[0] = out[1] = out[2] = v; return;
	}
	
	h /= 60.0; i = floor(h); f = h - i;
	
	p = v * (1.0 - s);
	q = v * (1.0 - s * f);
	t = v * (1.0 - s * (1.0 - f));
	
	switch (i) {
		case 0:
			r = v; g = t; b = p; break;
		case 1:
			r = q; g = v; b = p; break;
		case 2:
			r = p; g = v; b = t; break;
		case 3:
			r = p; g = q; b = v; break;
		case 4:
			r = t; g = p; b = v; break;
		case 5: default:
			r = v; g = p; b = q; break;
	}
	
	out[0] = r;
	out[1] = g;
	out[2] = b;
}



/***************** CMY based device dependant spaces ******************/

/* CMY to RGB */
void CMY2RGB(double in[], double out[])
{
	double c = in[0], m = in[1], y = in[2];
	
	out[0] = 1.0 - c;
	out[1] = 1.0 - m;
	out[2] = 1.0 - y;
}

/* RGB to CMY */
void RGB2CMY(double in[], double out[])
{
	double r = in[0], g = in[1], b = in[2];
	
	out[0] = 1.0 - r;
	out[1] = 1.0 - g;
	out[2] = 1.0 - b;
}


/* CMY to CMYK */
void CMY2CMYK(double in[], double out[])
{
	double c = in[0], m = in[1], y = in[2];
	double b = c;  if (m < b) b = m;  if (y < b) b = y;
	
	out[0] = (c - b)/(1.0 - b);
	out[1] = (m - b)/(1.0 - b);
	out[2] = (y - b)/(1.0 - b);
	out[3] = b;
}

/* CMYK to CMY */
void CMYK2CMY(double in[], double out[])
{
	double c = in[0], m = in[1], y = in[2], b = in[3];
	
	out[0] = c * (1.0 - b) + b;  if (out[0] > 1.0) out[0] = 1.0;
	out[1] = m * (1.0 - b) + b;  if (out[1] > 1.0) out[1] = 1.0;
	out[2] = y * (1.0 - b) + b;  if (out[2] > 1.0) out[2] = 1.0;
}



/***************** Arbitrary color space conversions ******************/

/* Return number of channels in color space */
int channels(icColorSpaceSignature any)
{
	switch (any) {
		case icSigGrayData:
			return 1;
		case icSigXYZData:
		case icSigLabData:
		case icSigLuvData:
		case icSigYxyData:
		case icSigRgbData:
		case icSigHsvData:
		case icSigCmyData:
			return 3;
		case icSigCmykData:
			return 4;
		default:
			error("[%s]: Color space not supported",
				ColorSpaceSignature2str(any));
	}
	
	return 0;
}

/* Convert string into ICC color space signature */
icColorSpaceSignature str2ColorSpaceSignature(char *s)
{
	if (!strcasecmp(s, "XYZ")) return icSigXYZData;
	if (!strcasecmp(s, "Lab")) return icSigLabData;
	if (!strcasecmp(s, "Luv")) return icSigLuvData;
	if (    !strcmp(s, "Yxy")) return icSigYxyData;
	if (!strcasecmp(s, "RGB")) return icSigRgbData;
	if (!strcasecmp(s, "GRAY")) return icSigGrayData;
	if (!strcasecmp(s, "HSV")) return icSigHsvData;
	if (!strcasecmp(s, "CMY")) return icSigCmyData;
	if (!strcasecmp(s, "CMYK")) return icSigCmykData;
	
	error("[%s]: Color space not supported", s); return 0;
}


/* Identity transformation of 3-channel data */
void identity3(double in[], double out[])
{
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}

/* Identity transformation of multi-channel data */
void identity(double in[], double out[], int channels)
{
	int i;
	
	for (i = 0; i < channels; i++) out[i] = in[i];
}

/* Gamma transformation of multi-channel data */
void gamma_scale(double in[], double out[], int channels, double gamma)
{
	int i;
	
	for (i = 0; i < channels; i++) out[i] = ppow(in[i], gamma);
}


/* Composite of two transformations */
#define composite(a, b, ab) void ab(double in[], double out[])\
		{ double tmp[MAXCHANNELS]; a(in, tmp); b(tmp, out); }

composite (HSV2RGB, RGB2XYZ, HSV2XYZ)
composite (XYZ2RGB, RGB2HSV, XYZ2HSV)

composite (CMY2RGB, RGB2XYZ, CMY2XYZ)
composite (XYZ2RGB, RGB2CMY, XYZ2CMY)

composite (CMYK2CMY, CMY2XYZ, CMYK2XYZ)
composite (XYZ2CMY, CMY2CMYK, XYZ2CMYK)


/* Return function transforming from any color space to CIE XYZ */
transform toXYZ(icColorSpaceSignature any)
{
	switch (any) {
		case icSigXYZData:
			return identity3;
		case icSigLabData:
			return Lab2XYZ;
		case icSigLuvData:
			return Luv2XYZ;
		case icSigYxyData:
			return Yxy2XYZ;
		case icSigRgbData:
			return RGB2XYZ;
		case icSigGrayData:
			return Gray2XYZ;
		case icSigHsvData:
			return HSV2XYZ;
		case icSigCmyData:
			return CMY2XYZ;
		case icSigCmykData:
			return CMYK2XYZ;
		default:
			error("[%s]: Color space not supported",
				ColorSpaceSignature2str(any));
	}
	
	return NULL;
}

/* Return function transforming from CIE XYZ to any color space */
transform fromXYZ(icColorSpaceSignature any)
{
	switch (any) {
		case icSigXYZData:
			return identity3;
		case icSigLabData:
			return XYZ2Lab;
		case icSigLuvData:
			return XYZ2Luv;
		case icSigYxyData:
			return XYZ2Yxy;
		case icSigRgbData:
			return XYZ2RGB;
		case icSigGrayData:
			return XYZ2Gray;
		case icSigHsvData:
			return XYZ2HSV;
		case icSigCmyData:
			return XYZ2CMY;
		case icSigCmykData:
			return XYZ2CMYK;
		default:
			error("[%s]: Color space not supported",
				ColorSpaceSignature2str(any));
	}
	
	return NULL;
}



/***************** Error metric on various spaces *********************/

/* CIE Lab */
double Lab_dE(double Lab1[], double Lab2[])
{
	double d[3] = {Lab2[0]-Lab1[0], Lab2[1]-Lab1[1], Lab2[2]-Lab1[2]};
	
	return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

/* CIE XYZ */
double XYZ_dE(double XYZ1[], double XYZ2[])
{
	double Lab1[3], Lab2[3];
	
	XYZ2Lab(XYZ1, Lab1);
	XYZ2Lab(XYZ2, Lab2);
	
	return Lab_dE(Lab1, Lab2);
}

/* Jacobian of CIE XYZ -> perceptual Lab transformation */
void dLab_dXYZ(double XYZ[], double dLab[3][3])
{
	double *I = XYZ_ILLUM, dx, dy, dz, dL;
	double x = XYZ[0]/I[0], y = XYZ[1]/I[1], z = XYZ[2]/I[2];
	
	if (x > 0.008856451586)
		dx = pow(x, -2.0/3.0)/3.0/I[0];
	else
		dx = 7.787036979/I[0];
	
	if (y > 0.008856451586) {
		dy = pow(y, -2.0/3.0)/3.0/I[1];
		dL = 116.0 * dy;
	} else {
		dy = 7.787036979/I[1];
		dL = 903.2963058 * dy;
	}
	
	if (z > 0.008856451586)
		dz = pow(z, -2.0/3.0)/3.0/I[2];
	else
		dz = 7.787036979/I[2];
	
	
	dLab[0][0] = 0.0;		/* dL/dX */
	dLab[0][1] = dL;		/* dL/dY */
	dLab[0][2] = 0.0;		/* dL/dZ */
	dLab[1][0] = 500.0*dx;		/* da/dX */
	dLab[1][1] = -500.0*dy;		/* da/dY */
	dLab[1][2] = 0.0;		/* da/dZ */
	dLab[2][0] = 0.0;		/* db/dX */
	dLab[2][1] = 200.0*dy;		/* db/dY */
	dLab[2][2] = -200.0*dz;		/* db/dZ */
}

/* Error metric on CIE XYZ space */
void XYZ_metric(double XYZ[], double g[3][3])
{
	int i, j;
	double J[3][3];
	
	dLab_dXYZ(XYZ, J);
	
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
		g[i][j] = J[0][i]*J[0][j] + J[1][i]*J[1][j] + J[2][i]*J[2][j];
}
