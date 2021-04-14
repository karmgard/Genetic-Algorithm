#include "test_fitness.h"

/* Test the program by using the distance from a 
 * random point in an N-dimensional space as the 
 * measure of fitness.
 */

static double *testPoint;

void initialize_test() {
  int i = 0;

  int number_of_genes = params->getInt("NUMBER_OF_GENES");
  float *pLo = params->getFloatArray("LOWER_LIMITS");
  float *pHi = params->getFloatArray("UPPER_LIMITS");

  if ( !pLo || !pHi ) {
    cerr << "Limit arrays do not exist!\n";

    float lower = params->getFloat("LOWER_LIMIT_ALL");
    float upper = params->getFloat("UPPER_LIMIT_ALL");

    if ( lower == 0.0 && upper == 0.0 ) {
      perror("Invalid limits on parameters!\n");
      exit (2);
    }

    pLo = (float *)malloc( number_of_genes * sizeof(float) );
    pHi = (float *)malloc( number_of_genes * sizeof(float) );

    for ( int i=0; i<number_of_genes; i++ ) {
      pLo[i] = lower;
      pHi[i] = upper;
    }

  }

  int verbose = params->getInt("VERBOSE");
  int accuracy = params->getInt("ACCURACY");  

  if ( !testPoint ) {
    testPoint = new double [number_of_genes];

    for (i=0; i<number_of_genes; i++)
      testPoint[i] = pLo[i] + randf()*(pHi[i] - pLo[i]);

    if ( verbose ) {
      printf("Seeking ");
      for ( int i=0; i<accuracy+4; i++ )
	printf(" ");
      printf("( ");
      for ( int i=0; i<number_of_genes; i++ )
	printf("%+0.*f, ", accuracy, testPoint[i]);
      printf("\b\b )\n");
    }
  }

  TEST_INITIALIZED = true;
  return;
}

void test_fitness( void *p ) {
  int i = 0;
  float fitness = 0.0f;

  individual *person = (individual *)p;

  if ( !TEST_INITIALIZED )
    initialize_test();

  for ( i=0; i<person->nGenes; i++ ) {
      fitness += 
	(testPoint[i] - person->gene[i])*(testPoint[i] - person->gene[i]);
  }

  if ( fitness >= 0.0f )
    fitness = sqrt( fitness );
  else
    fitness = person->max_fitness;

  if ( fitness < 0.0f || fitness > person->max_fitness )
    fitness = person->max_fitness;

  fitness = fround(fitness, person->accuracy);

  person->fitness = fitness;

  return;
}

void dumpTest( void *p ) {
  //individual *person = (individual *)p;
  return;
}
