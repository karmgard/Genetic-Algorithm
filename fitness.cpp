#include <fitness.h>

#include "test_fitness.cpp" 
#include "vckm.cpp"

static void (*outFunc)( void * );
static void (*fitFunc)( void * );
static threadpool *pool;
static unsigned short Nthreads;
void initialize_fitness_library( void ) {

  string FITNESS_FUNCTION = params->getString("FITNESS_FUNCTION");
  Nthreads = params->getUInt("NUM_THREADS");

  if ( FITNESS_FUNCTION == "TEST" ) {
    initialize_test();               // Call the test init function here in case we're running 
                                     // with multiple threads, otherwise it'll likely initialize
                                     // with different solutions for each thread. Which can be 
                                     // quite painful.
    fitFunc = test_fitness;
    outFunc = dumpTest;
  } else if ( FITNESS_FUNCTION.compare(0, 3, "CKM", 3) == 0 ) {
    fitFunc = ckm_fitness;
    outFunc = dumpMatrix;
  } else {
    cout << "Unable to find fitness function " << FITNESS_FUNCTION << "\n";
    exit (2);
  }

  // Are we going to run the fitness calculations in parallel?
  if ( Nthreads )
    pool = new threadpool( fitFunc, Nthreads );

  return;
}

void getFitness( void *person ) {

  if ( Nthreads )
    pool->enqueue(person);
  else
    (*fitFunc)(person);
  return;

}

void lock(void) {
  pool->queue_lock();
  return;
}

void unlock(void) {
  pool->queue_unlock();
  return;
}

void wait_for_threads( void ) {
  pool->wait_until_empty();
  return;
}

void outputIndividual( void *person ) {
  (*outFunc)(person);
  return;
}
