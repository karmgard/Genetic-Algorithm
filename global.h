#ifndef __GLOBAL_H
#define __GLOBAL_H

#define GLIB_VERSION_MAX_ALLOWED  GLIB_VERSION_2_32
#define GLIB_VERSION_MIN_REQUIRED GLIB_VERSION_2_26

#include <stdlib.h>
#include <parameters.h>
#include <limits>

#define MAX_UL_INT 0xffffffff
#define MAX_INT    0x7fffffff
const float ONE_OVER_RAND_MAX = 1/(RAND_MAX + 1.0f);
const float INFINITY = numeric_limits<float>::infinity();

extern parameters *params;

inline float randf() {
  return ONE_OVER_RAND_MAX*random();
}

#endif
