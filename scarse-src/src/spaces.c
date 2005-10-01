/* $Id: spaces.c,v 1.10 2005/10/01 01:18:40 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Generic color space conversions.
 * 
 * Copyright (C) 1999-2005 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <afrolov@stanford.edu>
 * 
 */

/* CREDITS:
 *   Originally based on Color Spaces FAQ by David Bourgin.
 *   CIE Lab code based on icclib & examples by Graeme W. Gill,
 *   with coefficients replaced by rationals to avoid discontinuity.
 *   Futher color space info from http://www.brucelindbloom.com/
 *   and http://www.aim-dtp.net/aim/technology/cie_xyz/cie_xyz.htm
 *   Chromatic adaptation code based on research paper by
 *   S. Susstrunk, J. Holm, and G. D. Finlayson, SPIE 4300 (2001).
 */

/* IMPORTANT NOTE:
 * ===============
 *   XYZ values used here are based on RELATIVE colorimetry, meaning they
 *   are transformed so that physical white point is always mapped to the
 *   same value - the so-called PCS illuminant (specified to be D50).
 *   Currently, this is done by Bradford chromatic adaptation transform.
 *   In addition, the max Y value we use is 1.0, not 100.0 as it's often set!
 */

/* TODO:
 *
 */

#include <string.h>
#include <math.h>
#include <icc.h>

#include "spaces.h"
#include "util.h"



/***************** Chromatic adaptation transforms ********************/

/* PCS illuminant and white point (in absolute XYZ representation) */
double XYZ_ILLUM[3] = {0.964203, 1.000000, 0.824905}; /* Always CIE illuminant D50 */
double   XYZ_WPT[3] = {0.950455, 1.000000, 1.089050}; /* Adobe RGB (D65) by default */

/* CAT matrices; default is Bradford (as used by Photoshop) */
static double XYZscale[3][3] = {
	{ 1.0000, 0.0000, 0.0000 },
	{ 0.0000, 1.0000, 0.0000 },
	{ 0.0000, 0.0000, 1.0000 }
};
static double vonKries[3][3] = {
	{ 0.3897, 0.6890,-0.0787 },
	{-0.2298, 1.1834, 0.0464 },
	{ 0.0000, 0.0000, 1.0000 }
};
static double Bradford[3][3] = {
	{ 0.8951, 0.2664,-0.1614 },
	{-0.7502, 1.7135, 0.0367 },
	{ 0.0389,-0.0685, 1.0296 }
};
static double Sharp[3][3] = {
	{ 1.2694,-0.0988,-0.1706 },
	{-0.8364, 1.8006, 0.0357 },
	{ 0.0297,-0.0315, 1.0018 }
};
static double CMCCAT[3][3] = {
	{ 0.7982, 0.3389,-0.1371 },
	{-0.5918, 1.5512, 0.0406 },
	{ 0.0008, 0.0239, 0.9753 }
};
static double M709[3][3] = {
	{ 3.0803,-1.5373,-0.5430 },
	{-0.9211, 1.8758, 0.0453 },
	{ 0.0528,-0.2040, 1.1511 }
};
static double ROMM[3][3] = {
	{ 1.2977,-0.2556,-0.0422 },
	{-0.5251, 1.5082, 0.0169 },
	{ 0.0000, 0.0000, 1.0000 }
};
static double Prime[3][3] = {
	{ 2.0016,-0.5576,-0.4440 },
	{-0.7997, 1.6627, 0.1371 },
	{ 0.0089,-0.0190, 1.0100 }
};


/* Calculate matrix M that maps XYZ space from illuminant 1 to illuminant 2 */
void XYZ_CAT(double IL1[3], double IL2[3], double M[3][3])
{
	double WPT1[3], WPT2[3], S1[3][3], S2[3][3], T[3][3];
	
	apply33(M_CAT, IL1, WPT1); diag33(WPT1, T); mult33(T, M_CAT, S1);
	apply33(M_CAT, IL2, WPT2); diag33(WPT2, T); mult33(T, M_CAT, S2);
	
	inv33(S1, T); mult33(T, S2, M);
}



/******************* CIE XYZ based spaces *****************************/

/* CIE XYZ to CIE Yxy */
void XYZ2Yxy(double in[], double out[])
{
	double X = in[0], Y = in[1], Z = in[2], S = X + Y + Z;
	
	out[0] = Y; out[1] = X/S; out[2] = Y/S;
}

/* CIE Yxy to CIE XYZ */
void Yxy2XYZ(double in[], double out[])
{
	double Y = in[0], x = in[1], y = in[2], z = 1.0 - (x + y);
	
	out[0] = (x/y)*Y; out[1] = Y; out[2] = (z/y)*Y;
}


/* CIE XYZ to perceptual Lab */
void XYZ2Lab(double in[], double out[])
{
	double *I = XYZ_ILLUM, fx, fy, fz;
	double x = in[0]/I[0], y = in[1]/I[1], z = in[2]/I[2];
	
	fx = (x > EPSILON) ? pow(x, GAMMA1) : (KAPPA * x + BETA)/ALPHA;
	fy = (y > EPSILON) ? pow(y, GAMMA1) : (KAPPA * y + BETA)/ALPHA;
	fz = (z > EPSILON) ? pow(z, GAMMA1) : (KAPPA * z + BETA)/ALPHA;
	
	out[0] = ALPHA * fy - BETA;
	out[1] = 500.0 * (fx - fy);
	out[2] = 200.0 * (fy - fz);
}

/* Perceptual Lab to CIE XYZ */
void Lab2XYZ(double in[], double out[])
{
	double L = in[0], a = in[1], b = in[2];
	double *I = XYZ_ILLUM, x, y, z, fx, fy, fz;
	
	fy = (L + BETA)/ALPHA;
	fx = a/500.0 + fy;
	fz = fy - b/200.0;
	
	y = (L > 8.0) ? pow(fy, GAMMA) : L/KAPPA;
	x = (fx > (8.0+BETA)/ALPHA) ? pow(fx, GAMMA) : (ALPHA*fx - BETA)/KAPPA;
	z = (fz > (8.0+BETA)/ALPHA) ? pow(fz, GAMMA) : (ALPHA*fz - BETA)/KAPPA;
	
	out[0] = x * I[0];
	out[1] = y * I[1];
	out[2] = z * I[2];
}


/* CIE XYZ to Luv */
void XYZ2Luv(double in[], double out[])
{
	double X = in[0], Y = in[1], Z = in[2];
	double *I = XYZ_ILLUM, y = Y/I[1], L, u, v, un, vn;
	
	L = (y > EPSILON) ? ALPHA*pow(y, GAMMA1) - BETA : KAPPA * y;
	
	u = 4.0 * X/(X + 15.0*Y + 3.0*Z);
	v = 9.0 * Y/(X + 15.0*Y + 3.0*Z);
	
	un = 4.0 * I[0]/(I[0] + 15.0*I[1] + 3.0*I[2]);
	vn = 9.0 * I[1]/(I[0] + 15.0*I[1] + 3.0*I[2]);
	
	out[0] = L;
	out[1] = 13.0 * L * (u - un);
	out[2] = 13.0 * L * (v - vn);
}

/* Luv to CIE XYZ */
void Luv2XYZ(double in[], double out[])
{
	double L = in[0], u = in[1], v = in[2];
	double *I = XYZ_ILLUM, y, un, vn, X, Y, Z;
	
	y = (L > 8.0) ? pow((L + BETA)/ALPHA, GAMMA) : L/KAPPA;
	
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


/* Gray to CIE XYZ */
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
static const double Adobe[6] = {
	0.6400, 0.3300,
	0.2100, 0.7100,
	0.1500, 0.0600
};
static const double Best[6] = {
	0.7347, 0.2653,
	0.2150, 0.7750,
	0.1300, 0.0350
};
static const double Beta[6] = {
	0.6888, 0.3112,
	0.1986, 0.7551,
	0.1265, 0.0352
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
static const double Don4[6] = {
	0.6960, 0.3000,
	0.2150, 0.7650,
	0.1300, 0.0350
};
static const double EBU[6] = {
	0.6400, 0.3300,
	0.2900, 0.6000,
	0.1500, 0.0600
};
static const double Ekta[6] = {
	0.6950, 0.3050,
	0.2600, 0.7000,
	0.1100, 0.0050
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
static const double ProPh[6] = {
	0.7347, 0.2653,
	0.1596, 0.8404,
	0.0366, 0.0001
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
static const double WideG[6] = {
	0.7347, 0.2653,
	0.1152, 0.8264,
	0.1566, 0.0177
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
	{"Best",		D50,	Best,	2.2},
	{"Beta",		D50,	Beta,	2.2},
	{"Bruce",		D65,	Bruce,	2.2},
	{"CIE",			IlE,	CIE,	2.2},
	{"ColorMatch",		D50,	P22,	1.8},
	{"Don4",		D50,	Don4,	2.2},
	{"ECI",			D50,	NTSC,	1.8},
	{"EktaSpace",		D50,	Ekta,	2.2},
	{"NTSC",		IlC,	NTSC,	2.2},
	{"PAL/SECAM",		D65,	EBU,	2.2},
	{"ProPhoto",		D50,	ProPh,	1.8},
	{"SMPTE-C",		D65,	SMPTE,	2.2},
	{"sRGB",		D65,	HDTV,	2.2},
	{"WideGamut",		D50,	WideG,	2.2},
	{NULL,			NULL,	NULL,	0.0}
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
	{ 0.609741, 0.205276, 0.149185 },
	{ 0.311111, 0.625671, 0.063217 },
	{ 0.019470, 0.060867, 0.744568 }
};

double M_XYZ2RGB[3][3] = {		/* Adobe RGB by default */
	{ 1.962529,-0.610675,-0.341372 },
	{-0.978754, 1.916152, 0.033418 },
	{ 0.028692,-0.140673, 1.349255 }
};


/* Set white point and RGB primaries */
void SetPrimaries(double xy[4][2])
{
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
	
	double E1[3][3], T[3][3], Q[3][3], Y[3];
	
	/* Calculate RGB -> absolute XYZ transfer matrix Q */
	inv33(E, E1); apply33(E1, W, Y); diag33(Y, T); mult33(E, T, Q);
	
	/* Map XYZ white point to standard illuminant using CAT */
	XYZ_CAT(W, XYZ_ILLUM, T); mult33(T, Q, M_RGB2XYZ);
	
	/* Calculate XYZ -> RGB transfer matrix */
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


/* RGB to HSV (h = [0,360], s = [0,1], v = [0,1]) */
void RGB2HSV(double in[], double out[])
{
	double min, max, delta, h, s, v;
	double r = in[0], g = in[1], b = in[2];
	
	min = r;  if (g < min) min = g;  if (b < min) min = b;
	max = r;  if (g > max) max = g;  if (b > max) max = b;
	
	if (max == 0.0) { out[0] = out[1] = out[2] = 0.0; return; }
	
	delta = max - min; v = max; s = delta/max;
		
	if (r == max)      h = (g - b)/delta;		// between yellow & magenta
	else if (g == max) h = 2.0 + (b - r)/delta;	// between cyan & yellow
	else               h = 4.0 + (r - g)/delta;	// between magenta & cyan
		
	h *= 60.0; if (h < 0.0) h += 360.0;
	
	out[0] = h;
	out[1] = s;
	out[2] = v;
}

/* HSV to RGB */
void HSV2RGB(double in[], double out[])
{
	int i; double f, p, q, t, r, g, b;
	double h = in[0], s = in[1], v = in[2];
	
	if (s == 0) { out[0] = out[1] = out[2] = v; return; }
	
	h /= 60.0; i = floor(h); f = h - i;
	
	p = v * (1.0 - s);
	q = v * (1.0 - s * f);
	t = v * (1.0 - s * (1.0 - f));
	
	switch (i%6) {
		case 0: r = v; g = t; b = p; break;
		case 1: r = q; g = v; b = p; break;
		case 2: r = p; g = v; b = t; break;
		case 3: r = p; g = q; b = v; break;
		case 4: r = t; g = p; b = v; break;
		case 5: r = v; g = p; b = q; break;
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
		case icSigXYZData:  return vcopy3;
		case icSigLabData:  return Lab2XYZ;
		case icSigLuvData:  return Luv2XYZ;
		case icSigYxyData:  return Yxy2XYZ;
		case icSigRgbData:  return RGB2XYZ;
		case icSigGrayData: return Gray2XYZ;
		case icSigHsvData:  return HSV2XYZ;
		case icSigCmyData:  return CMY2XYZ;
		case icSigCmykData: return CMYK2XYZ;
		default: error("[%s]: Color space not supported",
				ColorSpaceSignature2str(any));
	}
	
	return NULL;
}

/* Return function transforming from CIE XYZ to any color space */
transform fromXYZ(icColorSpaceSignature any)
{
	switch (any) {
		case icSigXYZData:  return vcopy3;
		case icSigLabData:  return XYZ2Lab;
		case icSigLuvData:  return XYZ2Luv;
		case icSigYxyData:  return XYZ2Yxy;
		case icSigRgbData:  return XYZ2RGB;
		case icSigGrayData: return XYZ2Gray;
		case icSigHsvData:  return XYZ2HSV;
		case icSigCmyData:  return XYZ2CMY;
		case icSigCmykData: return XYZ2CMYK;
		default: error("[%s]: Color space not supported",
				ColorSpaceSignature2str(any));
	}
	
	return NULL;
}


/* Convert string into ICC color space signature */
icColorSpaceSignature str2ColorSpaceSignature(char *s)
{
	if (!strcasecmp(s, "XYZ"))  return icSigXYZData;
	if (!strcasecmp(s, "Lab"))  return icSigLabData;
	if (!strcasecmp(s, "Luv"))  return icSigLuvData;
	if (    !strcmp(s, "Yxy"))  return icSigYxyData;
	if (!strcasecmp(s, "RGB"))  return icSigRgbData;
	if (!strcasecmp(s, "GRAY")) return icSigGrayData;
	if (!strcasecmp(s, "HSV"))  return icSigHsvData;
	if (!strcasecmp(s, "CMY"))  return icSigCmyData;
	if (!strcasecmp(s, "CMYK")) return icSigCmykData;
	
	error("[%s]: Color space not supported", s); return 0;
}

/* Return number of channels in color space */
int channels(icColorSpaceSignature any)
{
	switch (any) {
		case icSigGrayData: return 1;
		case icSigXYZData:
		case icSigLabData:
		case icSigLuvData:
		case icSigYxyData:
		case icSigRgbData:
		case icSigHsvData:
		case icSigCmyData:  return 3;
		case icSigCmykData: return 4;
		default: error("[%s]: Color space not supported",
				ColorSpaceSignature2str(any));
	}
	
	return 0;
}



/***************** Error metric on various spaces *********************/

/* CIE Lab */
double dE_Lab(double Lab1[], double Lab2[])
{
	double d[3] = {Lab2[0]-Lab1[0], Lab2[1]-Lab1[1], Lab2[2]-Lab1[2]};
	
	return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

/* CIE XYZ */
double dE_XYZ(double XYZ1[], double XYZ2[])
{
	double Lab1[3], Lab2[3];
	
	XYZ2Lab(XYZ1, Lab1);
	XYZ2Lab(XYZ2, Lab2);
	
	return dE_Lab(Lab1, Lab2);
}

/* Jacobian of CIE XYZ -> perceptual Lab transformation */
static void dLab_dXYZ(double XYZ[], double dLab[3][3])
{
	double *I = XYZ_ILLUM, dx, dy, dz;
	double x = XYZ[0]/I[0], y = XYZ[1]/I[1], z = XYZ[2]/I[2];
	
	dx = (x > EPSILON) ? GAMMA1*pow(x, GAMMA1-1.0) : KAPPA/ALPHA; dx /= I[0];
	dy = (y > EPSILON) ? GAMMA1*pow(y, GAMMA1-1.0) : KAPPA/ALPHA; dy /= I[1];
	dz = (z > EPSILON) ? GAMMA1*pow(z, GAMMA1-1.0) : KAPPA/ALPHA; dz /= I[2];
	
	dLab[0][0] = 0.0;		/* dL/dX */
	dLab[0][1] = ALPHA*dy;		/* dL/dY */
	dLab[0][2] = 0.0;		/* dL/dZ */
	dLab[1][0] = 500.0*dx;		/* da/dX */
	dLab[1][1] = -500.0*dy;		/* da/dY */
	dLab[1][2] = 0.0;		/* da/dZ */
	dLab[2][0] = 0.0;		/* db/dX */
	dLab[2][1] = 200.0*dy;		/* db/dY */
	dLab[2][2] = -200.0*dz;		/* db/dZ */
}

/* Error metric on CIE XYZ space */
void gXYZ(double XYZ[], double g[3][3])
{
	int i, j;
	double J[3][3];
	
	dLab_dXYZ(XYZ, J);
	
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
		g[i][j] = J[0][i]*J[0][j] + J[1][i]*J[1][j] + J[2][i]*J[2][j];
}
