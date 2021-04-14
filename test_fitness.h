#ifndef __TEST_FITNESS_H
#define __TEST_FITNESS_H

#include "global.h"

#include <utilities.h>
#include <parameters.h>

#include <math.h>

static bool  TEST_INITIALIZED;

void initialize_test( void );
void test_fitness( void * );
void dumpTest( individual * );

#endif
