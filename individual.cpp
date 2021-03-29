#include "individual.h"
#include <fitness.h>

extern char state[256];

/*-- Default instantiator, mostly for creating temporary individuals --*/
individual::individual( void ) {

  this->nGenes = params->getInt("NUMBER_OF_GENES");

  this->gene = new float [this->nGenes];
  if ( this->gene == NULL )
    exit(2);

  for (int i=0; i<this->nGenes; i++)
    this->gene[i] = 0.0f;

  pLo = params->getFloatArray("LOWER_LIMITS");
  pHi = params->getFloatArray("UPPER_LIMITS");

  if ( !pLo || !pHi ) {
    cerr << "Limit arrays do not exist!\n";

    float lower = params->getFloat("LOWER_LIMIT_ALL");
    float upper = params->getFloat("UPPER_LIMIT_ALL");

    if ( lower == upper ) {
      perror("Invalid limits on parameters!\n");
      exit (2);
    }

    float *lo = (float *)malloc( nGenes * sizeof(float) );
    float *hi = (float *)malloc( nGenes * sizeof(float) );

    for ( int i=0; i<nGenes; i++ ) {
      lo[i] = lower;
      hi[i] = upper;
    }
    params->setFloatArray("LOWER_LIMITS", lo, nGenes);
    params->setFloatArray("UPPER_LIMITS", hi, nGenes);

    pLo = params->getFloatArray("LOWER_LIMITS");
    pHi = params->getFloatArray("UPPER_LIMITS");
  }

  mutation_rate = params->getFloat("MUTATION_RATE");
  simple_mutation = params->getBool("MUTATE_SIMPLE");

  num_threads = params->getInt("NUM_THREADS");

  this->fitness = this->max_fitness = params->getFloat("MAX_FITNESS");
  this->accuracy = params->getInt("ACCURACY");
  this->count = -1;
  this->generation = 0;
  this->previous = NULL;
  this->next     = NULL;

  return;
};

/*-- Basic instantiator, for creating individuals --*/
individual::individual( individual **new_person, bool initialize , bool FIRST ) {

  this->nGenes = params->getInt("NUMBER_OF_GENES");

  this->gene = new float [this->nGenes];
  if ( this->gene == NULL )
    exit(3);

  this->max_fitness = params->getFloat("MAX_FITNESS");
  this->accuracy    = params->getInt("ACCURACY");

  mutation_rate = params->getFloat("MUTATION_RATE");
  simple_mutation = params->getBool("MUTATE_SIMPLE");

  pLo = params->getFloatArray("LOWER_LIMITS");
  pHi = params->getFloatArray("UPPER_LIMITS");

  if ( !pLo || !pHi ) {
    cerr << "Limit arrays do not exist!\n";

    float lower = params->getFloat("LOWER_LIMIT_ALL");
    float upper = params->getFloat("UPPER_LIMIT_ALL");

    if ( lower == upper ) {
      perror("Invalid limits on parameters!\n");
      exit (2);
    }

    float *lo = (float *)malloc( nGenes * sizeof(float) );
    float *hi = (float *)malloc( nGenes * sizeof(float) );

    for ( int i=0; i<nGenes; i++ ) {
      lo[i] = lower;
      hi[i] = upper;
    }
    params->setFloatArray("LOWER_LIMITS", lo, nGenes);
    params->setFloatArray("UPPER_LIMITS", hi, nGenes);

    pLo = params->getFloatArray("LOWER_LIMITS");
    pHi = params->getFloatArray("UPPER_LIMITS");
  }

  if ( initialize ) {
    for (int i=0; i<this->nGenes; i++)
      this->gene[i] = pLo[i] + randf()*(pHi[i] - pLo[i]);
    this->testFitness();
  }
  else {
    for (int i=0; i<this->nGenes; i++)
      this->gene[i] = 0.0f;
    this->fitness = max_fitness;
  }

  this->previous    = 0x0;
  this->next        = 0x0;
  this->count       = 0;
  this->generation  = 0;

  if ( FIRST ) {
    this->count = 1;
  } else if ( *new_person ) {
    (*new_person)->next = this;
    this->previous = *new_person;
    this->count = (*new_person)->count + 1;
  }

  *new_person = this;
  return;
}

/*-- Default destructor --*/
individual::~individual( void ) {
  delete [] gene;
  this->next = 0x0;
  this->previous = 0x0;
  this->count = 0;
  this->generation = 0;
  this->fitness = this->max_fitness;
  this->pLo = this->pHi = 0x0;

  return;
}

/* Create a new set of genes for this individual */
void individual::set_genes( void ) {
  for (int i=0; i<this->nGenes; i++)
    this->gene[i] = pLo[i] + randf()*(pHi[i] - pLo[i]);
  this->testFitness();
  return;
}

bool individual::isClone( individual * person ) {
  if ( this->fitness != person->fitness )
    return false;

  bool clone = true;
  for ( int i=0; i<this->nGenes; i++ )
    clone = clone && (this->gene[i] == person->gene[i]);
  return clone;
}

individual *individual::get_mate( int size, individual *population[] ) {
  int number = (int)(randf()*((float)size-1.0f));

  int counter = 0;

  bool condition = 
    population[number]->count != this->count &&
    number < size && !this->isClone(population[number]);  

  while ( !condition ) {
    number = (int)(randf()*((float)size-1.0f));

    condition = population[number]->count != this->count &&
      number < size && !this->isClone(population[number]);

    if ( counter++ > 5*abs(size) ) {
      fprintf(stderr, "\nUnable to find any mates in a population of %i\n", size);
      return NULL;
    }

  }

  return population[number];
}

individual *individual::make_baby( individual *mommy ) {
  static individual *baby = new individual();

  for (int i=0; i<this->nGenes; i+=2) {

    if (randf() > 0.50f) {
      baby->gene[i] = this->gene[i];

      if ( i+1 < this->nGenes )
	baby->gene[i+1] = mommy->gene[i+1];

    } else {
      baby->gene[i] = mommy->gene[i];

      if ( i+1 < this->nGenes )
	baby->gene[i+1] = this->gene[i+1];

    }
  }

  baby->generation = 0;
  
  baby->mutate();
  if ( !num_threads )  
    baby->testFitness();

  return baby;
}

/*-- Test the fitness of this individual --*/
void individual::testFitness( void ) {
  getFitness((void *)this);
  return;
}

void individual::mutate_simple( void ) {

  /* 
   * Toss a random number in [0,1].... as long as this number
   * is less than params->MUTATION_RATE, randomly select
   * one gene to mutate by replacement. Slightly better than
   * brute force testing each gene for mutation.
   */
  int counter = 1;
  while ( randf() < mutation_rate ) {   // Are we going to mutate?
    // Yuppers.... select an integer in [0,this->nGenes-1]
    unsigned int which = rand() %  this->nGenes;
    this->gene[which] = pLo[which] + randf()*(pHi[which] - pLo[which]);
    counter++;
  }

  this->testFitness();
  return;
}

/*-- Mutate this individuals DNA by flipping random bits --*/
void individual::mutate( void ) {

  /* If the mutation rate is <= 0.0, just bounce.
   * Nothing is going to happen anyhow, so don't
   * waste the CPU cycles for no result.
   */
  if ( mutation_rate <= 0.0f )
    return;

  /*-- If a simple mutation scheme was called for, call it here and bail --*/
  if ( simple_mutation ) {
    mutate_simple();
    return;
  }

  /* We need to iterate over the bits in the parameters structure
   * this will be a pointer to those bits since floats don't go
   * into bitwise operators nicely. unsignend longs are 32 bits,
   * same as floats, so this matches each of the parameters in
   * struct parameters
   */
  unsigned int *fltPointer;

  /*-- Save the original parameter, in case the mutation is really bad --*/
  float saveParam;
  bool stillborn = false;

  /* this is the number of bits that we'll look at in the parameter
   * structure. (sizeof returns bytes, not bits, hence the 8 bits/byte 
   * conversion factor)
   */ 
  const unsigned long int number_of_bits = 8*sizeof( typeof(this->gene[0]) );

  /* This yields a cummulative probability that 
   * 1 or more bits will get flipped in the routine
   */
  static float probability_per_bit = mutation_rate/((float)number_of_bits);
  static float probability_per_gene = mutation_rate/((float)this->nGenes);
  unsigned short int i;

  /* Randomly flip bits with some cummulative probability.
   * According to IEEE Standard 754,  a single precision 
   * floating point number is represented as ... 
   * bit 31 = sign, bits 30 - 23 = exponent, bits 22 - 0 = 
   * fraction with an exponent bias of 127 (Little Endian).
   */
  errno = 0;
  for (i=0; i<this->nGenes; i++) {

    // One test to see if something is getting flipped
    if ( randf() < probability_per_gene ) {

      // So at least one is getting changed....
      unsigned short int which = (unsigned short int)(randf()*number_of_bits);

      /*-- Cast the memory location of the first variable to an unsigned long --*/
      fltPointer = (unsigned int *)&this->gene[i];
      saveParam = this->gene[i];
      
      // Flip the bit
      *fltPointer ^= (1<<which); 

      // Take a shot at flipping the others
      while ( randf() < probability_per_bit ) {
	which = (unsigned short int)randf()*number_of_bits;
	*fltPointer ^= (1<<which); 
      }

      // Check for errors induced in the bit flipping
      stillborn = this->gene[i] != this->gene[i] ||
	fabs(this->gene[i]) == INFINITY ||
	this->gene[i] < pLo[i] || this->gene[i] > pHi[i];

      if ( stillborn ) {
	this->gene[i] = saveParam;
	i--;
      }

    } // End of if ( randf() < params->MUTATION_RATE )

  } // End of for (i=0; i<this->nGenes; i++)

  if ( errno ) {
    printf("Error %i in ", errno);
    perror("mutate_genes(): ");
    exit(errno);
  }
  
  return;
}

/*-- Make a complete copy of person, excluding the links -*/
void individual::copy( individual *person, bool deep ) {

  *this->gene = {0.0f};
  memcpy( this->gene, person->gene, this->nGenes*sizeof(float) );

  this->fitness = person->fitness;
  this->generation = person->generation;

    if ( deep ) {
      if ( person->next )
	this->next = person->next;
      else
	this->next = 0x0;

      if ( person->previous )
	this->previous = person->previous;
      else
	this->previous = 0x0;
    }

  return;
}

void individual::output(bool newline) {

  printf("%i %0.*f ( ", this->count, this->accuracy, this->fitness);
  for ( int i=0; i<this->nGenes; i++ )
    printf("%+0.*f, ", this->accuracy, this->gene[i]);
  printf("\b\b )\n");

  if ( this->previous )
    printf(" following %i",this->previous->count);
  else
    printf(" first ");
  if ( this->next )
    printf(" leading %i", this->next->count);
  else
    printf(" last");

  if ( newline )
    printf("\n");


  return;
}
