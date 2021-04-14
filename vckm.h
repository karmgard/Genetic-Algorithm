#ifndef __VCKM_H
#define __VCKM_H

#include "global.h"

#include <parameters.h>
#include <utilities.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>

/*-- Define the quark masses in GeV (as of 14 May 13)--*/
#define Mu 0.0023f
#define Md 0.0049f
#define Ms 0.095f
#define Mc 1.29f
#define Mb 4.18f
#define Mt 172.9f

/* V_CKM ranges as of 14 May 13 */
float V_Min[] = { 0.97413, 0.2246,  0.00335,
		  0.2245,  0.97329, 0.0403,
		  0.00842, 0.0396,  0.999107 };

float V_Max[] = { 0.97443, 0.2260,  0.00363,
		  0.2259,  0.97360, 0.0421,
		  0.00888, 0.0414,  0.999182 };

static char  VCKM_MATRIX;
static float VCKM_MAX_FITNESS;
static float VCKM_EXIT_LIMIT;

void ckm_fitness( void * );
void  dumpMatrix( void * );

#ifndef Mu
void Vckm(char which, complex *m, 
	  float C1u, float C1d, float C2u, float C2d,
	  float d1, float d2,
	  float mu, float md, float ms, float mc, float mb, float mt);
#else
void Vckm(char which, complex *m, 
	  float C1u, float C1d, float C2u, float C2d,
	  float d1, float d2);
#endif

#endif
