/* $Id: spaces.h,v 1.8 2005/09/30 20:38:49 afrolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Generic color spaces conversions - declarations.
 * 
 * Copyright (C) 1999-2005 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */


/**********************************************************************/

#ifndef __SPACES_H__
#define __SPACES_H__

#define MAXCHANNELS 4		/* Compile-time default; change to suit */
#define BLEND_GAMMA 2.2		/* Compile-time default; change to suit */


/* Lab conversion constants */
#define ALPHA 116.0
#define BETA 16.0
#define GAMMA 3.0
#define GAMMA1 (1.0/3.0)
#define EPSILON (216.0/24389.0)
#define KAPPA (24389.0/27.0)

/* Chromatic adaptation transform */
#define M_CAT Bradford

/* Illuminant, white point, and primaries */
extern double XYZ_ILLUM[3], XYZ_WPT[3];
extern double M_RGB2XYZ[3][3], M_XYZ2RGB[3][3];

void SetPrimaries(double xy[4][2]);
int LookupPrimaries(char *p, double dest[4][2], double *gamma);
void XYZ_CAT(double IL1[3], double IL2[3], double M[3][3]);

/* Color space conversions */
void XYZ2Yxy(double in[], double out[]);
void Yxy2XYZ(double in[], double out[]);
void XYZ2Lab(double in[], double out[]);
void Lab2XYZ(double in[], double out[]);
void XYZ2Luv(double in[], double out[]);
void Luv2XYZ(double in[], double out[]);
void Gray2XYZ(double *in, double out[]);
void XYZ2Gray(double in[], double *out);

void RGB2XYZ(double in[], double out[]);
void XYZ2RGB(double in[], double out[]);
void Gray2RGB(double *in, double out[]);
void RGB2HSV(double in[], double out[]);
void HSV2RGB(double in[], double out[]);

void CMY2RGB(double in[], double out[]);
void RGB2CMY(double in[], double out[]);
void CMY2CMYK(double in[], double out[]);
void CMYK2CMY(double in[], double out[]);

typedef void (*transform)(double in[], double out[]);
transform toXYZ(int any); transform fromXYZ(int any);

int str2ColorSpaceSignature(char *s);
int channels(int any);

/* Color space error metrics */
double dE_Lab(double Lab1[], double Lab2[]);
double dE_XYZ(double XYZ1[], double XYZ2[]);
void gXYZ(double XYZ[], double g[3][3]);


#endif /* __SPACES_H__ */
