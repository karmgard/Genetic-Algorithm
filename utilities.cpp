#include "utilities.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int UTILITY_ACCURACY;

/*-- Simple, locally available, function to round a floating point number --*/
double fround( double number, int digits ) {
  double rounded_number =  0.0f;

  if ( !UTILITY_ACCURACY )
    UTILITY_ACCURACY = params->getInt("ACCURACY");

  char   rounded_number_string[UTILITY_ACCURACY];
  char   format[8];
  
  sprintf( format, "%%.%if", digits );
  sprintf( rounded_number_string, format, number );

  rounded_number = strtod(rounded_number_string, NULL);

  return rounded_number;
}

/*-- Simple local function to return |z| rapidly --*/
float Cabs(complex z) { 
  float x, y, ans, temp; 

  x = fabs(z.r); 
  y = fabs(z.i); 

  if ( x == 0.0 ) 
    ans = y; 
  else if ( y == 0.0 ) 
    ans = x; 
  else if ( x > y ) { 
    temp = y/x; 
    ans = x*sqrt( 1.0 + temp*temp ); 
  } else { 
    temp = x/y; 
    ans = y*sqrt( 1.0 + temp*temp ); 
  }

  return ans; 
} 

