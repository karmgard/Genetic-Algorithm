#include "population.h"

/*-- We'll need a temporary population --*/
static population *newPopulation;

/*-- Base constructor. Creates a new population with randomly filled individuals --*/
population::population( void ) {

  this->first = new individual( &this->first, true, true );
  person = this->first;

  for (int i=1; i<params->INITIAL_POPULATION; i++)
    person = new individual( &person );

  this->last = person;
  this->recount();
  
  this->get_fittest();

  this->allocation = params->INITIAL_POPULATION;
  this->fitness_array = new double [this->allocation];

  this->generation = 0;
  this->clones = 0;
  this->average = 0.0f;
  this->stdev = 0.0f;
  this->variation = 0.0f;
  this->mating_in_progress = false;

  return;
}

/*-- Destructor for class population --*/
population::~population( void ) {

  if ( this != newPopulation ) {
    if (newPopulation)
      delete newPopulation;
  }
  this->flush();

  if ( allocation && fitness_array )
    delete [] fitness_array;

  this->person = this->first->next;
  delete first;

  while ( person ) {
    individual * temp = person;
    person = person->next;
    delete temp;
  }

  return;
}

void population::copy_elites(void) {

  /*-- Keep the most fit individual(s) for a few generations, by request --*/
  if (params->ELITISM_GENERATIONS > 0) {

    /*-- If we're carrying some % of the population over... pack 'em in --*/
    if ( params->PERCENT_ELITES_KEPT ) {
      int number_to_keep = (int)fround((params->PERCENT_ELITES_KEPT*this->last->count), 0);

      /*-- Start from the first person in the old population --*/
      person = this->first;

      while ( person->count < number_to_keep ) {

	if ( person->generation++ < params->ELITISM_GENERATIONS ) {
	  newPopulation->unshift(person);

	  if ( params->KEEP_STABLE_POPULATION )
	    newPopulation->pop();

	} else
	  person->generation = 0;

	person = person->next;

      }

    } else {              // We're only keeping the most fit individual
      if ( this->mostfit->generation++ < params->ELITISM_GENERATIONS ) {

	newPopulation->unshift(this->mostfit);

	if ( params->KEEP_STABLE_POPULATION )
	  newPopulation->pop();

      } else
	this->mostfit->generation = 0;
    }
  }

  return;
}

void population::roulette_fill( void ) {

  /*
   * Roullette fill: The number of children an individual can have
   * is how far above the average fitness of the population they are
   * (survival of the fittest) + a chance at 1 more.
   */
  individual *person = this->first;
  float average_fitness = params->MAX_FITNESS - get_avg_fitness();

  if ( average_fitness == 0.0f ) {

    while (person) {
      person->progeny = 1;
      person = person->next;
    }
    if ( params->VERBOSE == 3 )
      printf("No spread... everyone gets a baby\n");
    return;
  }

  float population_control = 5*params->INITIAL_POPULATION/(this->last->count);

  while ( person ) {
    float fitness = params->MAX_FITNESS - person->fitness;
    unsigned int number_of_copies = (unsigned int)(fitness/average_fitness);
    float chance = (fitness/average_fitness) - (int)(fitness/average_fitness);

    person->progeny = number_of_copies;
    if ( params->KEEP_STABLE_POPULATION ) {
      if ( randf() < chance )
	person->progeny++;
    } else {
      if ( randf() < chance*population_control )
	person->progeny++;
    }
    if ( params->VERBOSE == 3 )
      printf("Individual %i has %i progeny\n", person->count, person->progeny);

    person = person->next;
  }
  return;
}

void population::mate( void ) {
  /*-- Create some temporary individuals for mating --*/
  individual *daddy;
  individual *mommy;
  individual *baby;

  // Raise the mating flag.... Ahoy maties :P
  this->mating_in_progress = true;

  /*-- Declare the internal variables we'll need to do our job --*/
  int newCount = 0;

  if ( params->VERBOSE == 3 ) {
    printf("/================ mating population ===================/\n");
    this->dump();
    printf("/======================================================/\n");
  }

  if ( !newPopulation )
    newPopulation = new population();

  newPopulation->mating_in_progress = true;

  // Figure out how many kids each individual can have
  this->roulette_fill();

  /*-- New individuals are poked onto the new population --*/
  baby = newPopulation->first;
  daddy = this->first;

  // Build an array of pointers to the population
  mommy = this->first;
  individual *pop_array [this->last->count];
  while ( mommy ) {
    pop_array[newCount++] = mommy;
    mommy = mommy->next;
  }
  newCount = 1;

  while ( daddy ) {

    if ( daddy->progeny <= 0 ) {
      daddy = daddy->next;
      continue;
    }

    mommy = daddy->get_mate(this->last->count, pop_array);

    if ( daddy == NULL || mommy == NULL )
      break;
    else if ( daddy->count == mommy->count ) {
      fprintf(stderr,"\nSelf replication!\n");
      fprintf(stderr,"Mating %i with %i\n", daddy->count, mommy->count);

    } else if ( daddy->fitness == mommy->fitness ) {
      // Could be cloning...
      bool cloning = true;
      for ( int i=0; i<params->NUMBER_OF_GENES; i++ )
	cloning = cloning && (daddy->gene[i] == mommy->gene[i]);

      if ( cloning )
	fprintf(stderr,"\nCloning!\nMating %i,%f with %i,%f\n",
	     daddy->count, daddy->fitness, mommy->count, mommy->fitness);
    }

    if ( params->VERBOSE == 3 ) {
      daddy->output(true);
      mommy->output(true);
    }

    if ( baby == NULL ) {
      if ( params->KEEP_STABLE_POPULATION )
	break;

      if ( params->VERBOSE == 3 )
	printf("Adding brand new baby\n");

      baby = newPopulation->push(daddy->make_baby( mommy ));

    } else
      baby->copy( daddy->make_baby( mommy ) );

    if ( params->VERBOSE == 3 )
      baby->output(true);

    newCount++;
    baby = baby->next;
    daddy->progeny--;

    if ( daddy->progeny <= 0 )
      daddy = daddy->next;

    if ( newCount >= params->INITIAL_POPULATION && params->KEEP_STABLE_POPULATION )
      break;
  }

  // Get the fitness of each member of the population
  if ( params->NUM_THREADS ) {
    lock();
    baby = newPopulation->first;
    while ( baby ) {
      getFitness((void *)baby);
      baby = baby->next;
    }
    unlock();
    wait_for_threads();
  }

  if ( newCount < newPopulation->last->count ) {
    newPopulation->trim( newCount );
    newPopulation->recount();
  }

  if ( params->VERBOSE == 3 ) {
    printf("/================ new population ======================/\n");
    newPopulation->dump();
    printf("/======================================================/\n");
  }


  // Keep the elitist bastards around for a while
  if (params->ELITISM_GENERATIONS > 0)
    this->copy_elites();

  // Either way we go, we'll need a sorted population
  newPopulation->sort();

  //newPopulation->get_fittest();

  /*-- Now, copy the new population into the old one --*/
  this->copy(newPopulation);

  // Spring is over.... enter the summer of our life
  this->mating_in_progress = false;
  newPopulation->mating_in_progress = false;

  if ( params->VERBOSE == 3 ) {
    individual *newP = newPopulation->first;
    individual *oldP = this->first;
    printf("\nGeneration %i complete\n", this->generation);
    while ( newP && oldP ) {
      printf("newPop %02i = %f\tPop %02i = %f\n", newP->count, newP->fitness, oldP->count, oldP->fitness);
      newP = newP->next;
      oldP = oldP->next;
    }
    printf("/======================================================/\n\n");
  }

  return;
}

void population::copy( population *newPop ) {

  individual *oldP, *newP;

  oldP = this->first;
  newP = newPop->first;

  while ( newP ) {

    if ( oldP == NULL ) {
      if ( params->KEEP_STABLE_POPULATION ) {
	break;
      }
      oldP = this->push( newP );
    } else
      oldP->copy( newP );

    oldP = oldP->next;
    newP = newP->next;

  }

  if ( newPop->last->count < this->last->count && !params->KEEP_STABLE_POPULATION )
    this->trim(newPop->last->count);

  this->count = 0;
  this->clones = 0;
  this->average = 0.0;
  this->stdev = 0.0;
  this->variation = 0.0;

  /*-- Keep a pointer to the most fit individual for efficiency --*/
  this->get_fittest();
  this->generation++;
  this->count = this->last->count;
  this->check_for_clones();
  this->get_statistics();

  if ( this->generation && !(this->generation % 50 ) )
    this->mutation_gain();

  return;
}

double *population::get_population_fitness() {

  if ( this->allocation < this->count ) {
    this->allocation = this->count;
    if ( fitness_array )
      delete [] this->fitness_array;
    this->fitness_array = new double [this->allocation];
  }
  memset(this->fitness_array, 0, this->allocation*sizeof(double));

  this->person = this->first;
  int counter = 0;
  while (this->person) {
    fitness_array[counter++] = this->person->fitness;
    this->person = this->person->next;
  }

  return fitness_array;
}


void population::check_for_clones() {

  individual *person1 = this->first;
  individual *person2 = this->first->next;

  int number_of_clones = 0;

  while ( person2 ) {
    number_of_clones += (person1->isClone(person2)) ? 1 : 0;
    person1 = person2;
    person2 = person2->next;
  }

  this->clones = number_of_clones;
  return;
}

void population::print(void) {

  if ( params->VERBOSE > 0 && params->VERBOSE <= 2 ) {
    printf("Most fit %0.*f ", 
	   params->ACCURACY, this->mostfit->fitness);

    printf("Generation %i ", this->generation);

    if ( params->VERBOSE == 1 )
      printf("     \r");
    else if ( params->VERBOSE == 2 ) {
      printf(" Pop. %i (%i clones) Stats: Avg = %0.1f StDev = %0.1f Var = %0.1f persist = %i rate = %0.2f ",
	     this->count, this->clones, this->average, this->stdev, this->variation, this->mostfit->generation,
	     params->MUTATION_RATE);
    }

    fflush(stdout);
  }
  return;
}

void population::dump( int max ) {
  individual *temp = this->first;
  int counter = 0;

  if ( !max )
    max = this->last->count;

  while ( temp && counter++ < max ) {
    printf("%03i (", temp->count);
    for ( int i=0; i<params->NUMBER_OF_GENES; i++ )
	    printf(" %+.*f,", params->ACCURACY, temp->gene[i]);
	  printf("\b )\t=> fitness = %.*f\n", params->ACCURACY, temp->fitness);
    temp = temp->next;
  }
  return;
}

/*-- Zero all of the individual parameters in the population --*/
void population::flush( void ) {

  this->person = this->first;
  while ( this->person != NULL ) {
    *this->person->gene = {0.0f};
    this->person->fitness = params->MAX_FITNESS;

    this->person = this->person->next;
  }
  return;
}

individual * population::get_individual( int who ) {

  if ( who >= this->last->count )
    return this->last;
  else if ( who == 1 )
    return this->first;

  int counter = 0;
  individual *whichPerson = this->first;

  while ( whichPerson != NULL && ++counter < who )
    whichPerson = whichPerson->next;

  return whichPerson;
}

/*-- Adjust the size of a population --*/
void population::trim( int newSize ) {

  this->last = this->get_individual(newSize);
  this->person = this->last->next;

  individual *temp;
  while ( person ) {
    temp = this->person;
    this->person = this->person->next;
    delete temp;
  }
  this->last->next = NULL;

  return;
}

/*-- Remove the head of the list --*/
individual *population::shift(void) {
  person = this->first;
  this->first = this->first->next;
  this->first->previous = NULL;

  if ( !this->mating_in_progress )
    this->recount();

  delete person;

  return NULL;
}

/*-- Add a new (empty) member to the head of the population --*/
individual *population::unshift(void) {

  person = new individual();
  this->first->previous = person;
  person->next = this->first;
  this->first = person;

  if ( !this->mating_in_progress )
    this->recount();
  return person;
}

/*-- Add a new member to the head of the population --*/
individual *population::unshift(individual *head) {

  person = new individual();
  person->copy(head);
  this->first->previous = person;
  person->next = this->first;
  this->first = person;

  if ( !this->mating_in_progress )
    this->recount();

  return person;
}

/*-- Pop the last individual off the end of the list --*/
individual *population::pop( void ) {
  person = this->last;
  this->last = this->last->previous;
  this->last->next = NULL;

  person->previous = NULL;
  person->next = NULL;
  person->fitness = params->MAX_FITNESS;

  delete person;

  return NULL;
}

/*-- Add a new member to the end of the population --*/
individual *population::push( individual *tail ) {

  person = this->last;
  person = new individual();

  person->copy(tail);

  this->last->previous = person;
  this->last->next     = person;

  person->previous = this->last;
  person->next = NULL;
  person->count = person->previous->count+1;

  this->last = person;
  
  return person;
}

/*-- Add a new (empty) member to the end of the population --*/
individual *population::push( void ) {

  person = this->last;
  person = new individual( &person );

  this->last->previous = person;
  this->last->next     = person;

  person->previous = this->last;
  person->next = NULL;
  person->count = person->previous->count+1;

  this->last = person;

  return person;
}

int population::get_count( void ) {
  return this->last->count;
}

/*-- Get some population statistics about the fitness --*/
void population::get_statistics( void ) {

  float avg   = this->get_avg_fitness();
  float stdev = this->get_stdev_fitness();

  float var;
  if ( avg )
    var = stdev/avg;
  else
    var = MAX_INT;

  // Unbiased estimator of the coefficient of variation for normal populations
  this->variation = (1 + 1/(4*this->count))*var;

  return;
}

float population::get_avg_fitness( void ) {

  this->average = 0.0f;

  person = this->first;
  while ( person != NULL ) {
    this->average += person->fitness;
    person = person->next;
  }

  this->average /= this->last->count;

  return this->average;
}

float population::get_stdev_fitness( void ) {

  this->stdev = 0.0f;
  person = this->first;
  while ( person != NULL ) {
    this->stdev += (person->fitness - this->average)*(person->fitness - this->average);
    person = person->next;
  }

  if ( this->stdev >= 0.0f )
    this->stdev = sqrt(this->stdev/this->last->count);
  else
    this->stdev = -1.0f;

  return this->stdev;
}

void population::get_fittest( void ) {
  this->mostfit = this->first;
  return;

  float most_fit = params->MAX_FITNESS;
  this->mostfit = this->first;

  person = this->first;
  while ( person != NULL ) {
    if ( person->fitness < most_fit ) {
      most_fit = person->fitness;
      this->mostfit = person;
    }
    person = person->next;
  }
  if ( most_fit > params->MAX_FITNESS )
    most_fit = params->MAX_FITNESS;

  if ( params->VERBOSE == 3 && this->mostfit != this->first )
    cout << "\n######################## MOST FIT NOT FIRST! ########################\n";

  return;
}

void population::recount( void ) {
  int counter = 1;

  person = this->first;
  while ( person != NULL ) {
    person->count = counter++;
    person = person->next;
  }

  return;
}

/*-- Modify the rate of mutation based on population distribution --*/
void population::mutation_gain( void ) {

  if ( params->MUTATION_RATE >= 0.5 ) {
    params->MUTATION_RATE = 0.5f;
    return;
  }

  if ( this->mostfit->generation > 500 ) {
    if ( params->VERBOSE > 1 )
      cout << "\nincreasing rate for persistance of " << this->mostfit->generation << "\n";
    params->MUTATION_RATE += params->MUTATION_GAIN;
  }

  if ( this->variation < 0.5f ) {
    if ( params->VERBOSE > 1 )
      cout << "\nincreasing rate for variation of " << this->variation << "\n";
    params->MUTATION_RATE += params->MUTATION_GAIN;
  }

  if  ( this->clones > 0.025*this->last->count ) {
    if ( params->VERBOSE > 1 )
      cout << "\nincreasing rate for cloning rate of " << (float)this->clones/(float)this->last->count << "\n";
    params->MUTATION_RATE += params->MUTATION_GAIN;
  }

  if ( this->variation > 5.0f ) {
    if ( params->VERBOSE > 1 )
      cout << "\ndecreasing rate for variation of " << this->variation << "\n";
    params->MUTATION_RATE -= params->MUTATION_GAIN;
  }

  if ( params->MUTATION_RATE > 0.5 )
    params->MUTATION_RATE = 0.5f;

  return;
}

// Sort ascending by fitness
void population::sort() {

  // Make sure we've got a good count
  this->recount();

  individual *p = this->first;
  unsigned int n = this->last->count;
  unsigned long i = 0;
  individual * popArray [n+1];

  while ( p ) {
    popArray[ ++i ] = p;
    p = p->next;
  }

  if ( params->SORT_TYPE == "QUICK" )
    this->quick_sort((void **)popArray);
  else if ( params->SORT_TYPE == "HEAP" )
    this->heap_sort((void **)popArray);
  else {
    fprintf(stderr, "\nUnknown sort routine\nDefaulting to heap sort\n");
    params->SORT_TYPE = "HEAP";
    this->heap_sort((void **)popArray);
  }

  for ( i=1; i<=n; i++) {
    popArray[i]->count = i;
    if ( i==1 ) {
      popArray[1]->previous = NULL;
      popArray[1]->next = popArray[2];
      this->first = this->mostfit = popArray[1];

    } else if ( i == n ) {
      popArray[n]->next = NULL;
      popArray[n]->previous = popArray[n-1];
      this->last = popArray[n];

    } else {
      popArray[i]->previous = popArray[i-1];
      popArray[i]->next = popArray[i+1];

    }
  }

  return;
}

void population::heap_sort(void **p) {

  unsigned int n = this->last->count;
  unsigned long i = 0, ir = 0,j = 0,l = 0;
  individual *rra = NULL;

  individual **ra = (individual **)p;

  if ( n < 2 )
    return;

  l = (n >> 1) + 1;
  ir = n;

  for (;;) {
    if ( l > 1 ) {
      rra = ra[--l];
    } else {
      rra = ra[ir];           // clear a space at the end of the array
      ra[ir] = ra[1];         // retire the top of the heap into it
      if ( --ir == 1 ) {
	ra[1] = rra;          // bottom of the heap
	break;
      }
    }

    i = l;
    j = l + l;
    while ( j <= ir ) {
      if ( j < ir && ra[j]->fitness < ra[j+1]->fitness )
	j++;
      if ( rra->fitness < ra[j]->fitness ) {
	ra[i] = ra[j];        // If jth person is fitter than ith person, swap 'em
	i = j;
	j <<= 1;
      } else
	break;
    }
    ra[i] = rra;
  }

  return;
}

#define SWAP(a,b) {individual *temp=a;a=b;b=temp;}

void population::quick_sort(void **p) {

  individual **arr = (individual **)p;

  unsigned int n = this->last->count;
  const unsigned short NSTACK = (n+1);
  unsigned long i=0,ir=n,j=0,k=0,l=1,istack[NSTACK];
  int jstack=0;
  individual *a;

  const unsigned short M = 8;

  for (;;) {
    if ( ir - l < M ) {        // Use straight insertion for small subarays
      for ( j=l+1; j<=ir; j++ ) {
	a = arr[j];
	for ( i=j-1; i>=l; i-- ) {
	  if ( arr[i]->fitness <= a->fitness )
	    break;
	  arr[i+1] = arr[i];
	}
	arr[i+1] = a;
      }
      if ( jstack == 0 )
	break;

      ir = istack[jstack--];   // Pop the stack and start a new round of partitioning
      l  = istack[jstack--];

    } else {
      k = (l+ir) >> 1;         // Pivot = median of (left,center,right)
      SWAP(arr[k],arr[l+1]);   // Move the pivot out of the way

      // swap so a[l]<=a[l+1]<=a[ir]
      if ( arr[l]->fitness > arr[ir]->fitness )
	SWAP(arr[l],arr[ir]);

      if ( arr[l+1]->fitness > arr[ir]->fitness )
	SWAP(arr[l+1],arr[ir]);

      if ( arr[l]->fitness > arr[l+1]->fitness )
	SWAP(arr[l],arr[l+1]);

      i=l+1;                   // Initialize the pointers for partitioning
      j=ir;
      a=arr[l+1];              // Pivot

      for (;;) {               // Innermost loop

	do i++; while (arr[i]->fitness < a->fitness);  // Scan up   to find first element > a
	do j--; while (arr[j]->fitness > a->fitness);  // Scan down to find first element < a

	if ( j < i )           // Pointers crossed each other, all done
	  break;
	SWAP(arr[i],arr[j]);   // Exchange the elements
      }                        // End of inner loop

      arr[l+1] = arr[j];       // Insert the pivot into place

      arr[j] = a;
      jstack += 2;

      // Take the pointers to the larger of the subarrays
      // and push them onto the stack. Sort the smaller 
      // partition immediately
      if ( jstack > NSTACK ) {
	fprintf(stderr, "Stack size %i in QuickSort is too small: need %i\n",
		NSTACK, jstack);
      }

      if ( ir-i+1 >= j-1 ) {
	istack[jstack]   = ir;
	istack[jstack-1] = i;
	ir = j - 1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }

  return;
}

