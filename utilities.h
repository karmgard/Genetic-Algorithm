#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <global.h>

/*-- Define a couple of new (personal) variable types --*/
typedef struct { float r, i; } complex;

double fround( double number, int digits );
float Cabs( complex z );

#endif
