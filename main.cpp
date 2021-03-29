#include "global.h"
#include "population.h"
#include "individual.h"
//#include "gnuplot.h"
#include <fitness.h>

#include <signal.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <values.h>
#include <fcntl.h>

#include <glibtop.h>
#include <glibtop/cpu.h>
#include <glibtop/mem.h>
#include <glibtop/proctime.h>

/*-- Prototypes --*/
static int naptime( int );
//void randomize();

/*-- Global statements --*/
char state[512];
parameters *params;

static bool STOPNOW = false;
static int  renice_priority = 15;

/*-- Install a handler for SIGSTOP so we can exit cleanly on <CTRL>-c --*/
static void sig_stop(int sig_number) {
  STOPNOW = true;
}

int main( int argc, char** argv ) {

  struct timeval tv1, tv2;
  unsigned int elapsed_time = 0;

  //GnuPlot *gplot = new GnuPlot();
  //double * fitness_array
  //char range[128];

  const unsigned int nbins = 100;
  double * ordinate = new double[nbins];

  /*-- Instantiate the requisite classes --*/
  params = new parameters( (char *)"ga.rcp" );

  int verbose = params->getInt("VERBOSE");
  int accuracy = params->getInt("ACCURACY");
  int cpu_usage_limit = params->getInt("CPU_USAGE_LIMIT");
  int maximum_generations = params->getInt("MAXIMUM_GENERATIONS");
  float exit_limit = params->getFloat("EXIT_LIMIT");
  int dumpn_top = params->getInt("DUMP_N_TOP");
  bool show_plot = params->getBool("SHOW_PLOT");
  //int plot_freq = params->getInt("PLOT_FREQ");

  unsigned long int seed = params->getULong("SEED");

  // Initialize the random number generator
  randomize(seed);

  // Initialize the function mapping for the fitness library
  initialize_fitness_library();

  // Create a new population of params->INITIAL_POPULATION individuals
  population *society = new population();

  /*-- Since this is a CPU intensive process, renice it to low priority --*/
  setpriority( PRIO_PROCESS, 0, renice_priority );

  if ( verbose == 2 )
    gettimeofday(&tv1, NULL);

  while (!STOPNOW) {

    /*-- Start the mating dance (it must be springtime!) --*/
    society->mate();

    /*-- Allow a clean exit on <CTRL>-C (SIGINT) --*/
    if(signal(SIGINT, sig_stop) == SIG_ERR)
      perror("error catching signal: ");

    // Dump out the population status
    society->print();

    if ( verbose == 2 ) {
      if ( elapsed_time > 0 )
	printf("Gen/s = %.0f     \r", (double)(society->generation/elapsed_time));
      else
	printf("Gen/s = x.xx     \r");
    }

    /***
    // pop a plot into a gnuplot window
    if ( show_plot && !(society->generation % plot_freq) ) {
      gplot->gnuplot_resetplot();

      float avg = society->average;
      float stdev = society->stdev;

      sprintf(range, "set xrange [%f:%f];set yrange[0:50]", 
	      exit_limit, avg+stdev);
      gplot->gnuplot_cmd(range);

      float plot_max = avg + 2*stdev;
      if ( plot_max < 0.25 )
	plot_max = 0.25;
      for ( unsigned int i=0; i<nbins; i++ )
	ordinate[i] = i*plot_max/nbins;

      fitness_array = society->get_population_fitness();
      gplot->gnuplot_plot_histogram( ordinate, fitness_array, nbins, 0, (char *)"Population Fitness");
    }
    ***/

    // Take a little siesta to reduce CPU consumption
    if ( cpu_usage_limit < 100 )
      naptime(society->generation);

    STOPNOW = 
      STOPNOW || 
      society->mostfit->fitness <= exit_limit || 
      (maximum_generations > 0 && 
       (int)society->generation >= maximum_generations);

    if ( verbose == 2 ) {
      gettimeofday(&tv2, NULL);
      elapsed_time = (tv2.tv_sec - tv1.tv_sec);
    }
  } // End while (!STOPNOW)

  // Dump out the results
  printf("\n\nGeneration %i Most fit = %0.*f\n",
	 society->generation, accuracy, society->mostfit->fitness);

  outputIndividual(society->mostfit);

  // And the top N individuals
  if ( dumpn_top > 0 )
    society->dump(dumpn_top);

  if ( show_plot ) {
    char temp;
    printf("Hit <RET> to finish: ");
    temp = getc(stdin);
    temp++;
  }

  delete society;
  delete params;
  //delete gplot;
  delete [] ordinate;

  return 0;
}

static int naptime( int counter ) {

  static pid_t PID = 0;
  static glibtop_proc_time proc_time;
  static glibtop_cpu       cpu_data;

  static unsigned long this_cpu,  last_cpu,  total_cpu;
  static unsigned long this_proc, last_proc, total_proc;
  static unsigned long orig_cpu,  orig_proc;
  static float proc_load;
  static int nap_time = 10, sleep_time = 500;
  static int proc_percent = 0;

  static int cpu_limit = params->getInt("CPU_USAGE_LIMIT");
  static int verbose = params->getInt("VERBOSE");

  /*-- First, figure out who we are this time --*/
  if ( !PID ) {
    PID = getpid();

    /*-- Get the initial data from libgtop --*/
    glibtop_get_cpu        (&cpu_data);
    glibtop_get_proc_time  (&proc_time, PID);

    orig_cpu  = last_cpu  = cpu_data.total;
    orig_proc = last_proc = proc_time.utime;
  }

  if ( !(counter%nap_time) && cpu_limit < 100 ) {

    /*-- Grab the data from GLibTop --*/
    glibtop_get_cpu        (&cpu_data);
    glibtop_get_proc_time  (&proc_time, PID);

    this_cpu  = cpu_data.total;
    this_proc = proc_time.utime;

    total_cpu  = this_cpu - last_cpu;
    total_proc = this_proc - last_proc;

    if ( total_cpu != 0 )
      proc_load = (float)total_proc/(float)total_cpu;
    else
      proc_load = 1.0f;
      
    proc_percent = (int)(100.0f*proc_load);
  
    /*-- Figure out how often to go to sleep, based on CPU usage --*/
    if ( proc_percent > cpu_limit )
      nap_time = (int) ( (float)nap_time * ( (float)cpu_limit/(float)proc_percent ) );
    else if ( proc_percent < cpu_limit )
      nap_time = (int)((float)nap_time / (float)(cpu_limit/100.0f)/proc_load);

    if (nap_time < 1)
      nap_time = 1;

    /* Use an overall average to balence the CPU load
     * if we're relying on how often to enter the sleep cycle
     * because no one cycle uses enough CPU to get an
     * accurate read
     */
    last_cpu = orig_cpu;last_proc = orig_proc;

    /* If we're sleeping on every cycle and still overusing the CPU it's hopeless.
     * Start jinking about with the sleep time to bring the CPU under control.
     */
    if ( nap_time == 1 && proc_percent/10 != cpu_limit/10 ) {
      if ( proc_percent > cpu_limit ) {

	sleep_time = (int)((float)sleep_time/((float)cpu_limit/(float)proc_percent ) );
	if ( sleep_time < 50 )
	  sleep_time = 50;
      
      } else if ( proc_percent < cpu_limit ) {

	sleep_time = (int)((float)sleep_time*(float)(cpu_limit/100.0f)/proc_load);
	if ( sleep_time < 50 )
	  sleep_time = 50;

      }
      /* In this case, single cycles use enough CPU to measure,
       * so measure our usage on the fly instead of with an
       * overall average.
       */ 
      last_cpu = this_cpu;last_proc = this_proc;

    }

    fflush(stdout);

    /*-- Pause for a breath, now & then --*/
    usleep(sleep_time);
  }

  if ( verbose == 2 )
    printf(" Load = %i     \r", proc_percent);

  return proc_percent;
} // End static void naptime()

/***
void randomize( void ) {
  unsigned long int seed = params->getLong("SEED");
  int filedes = 0;
  extern char state[512];

  if ( !seed ) {
    errno = 0;
    if ( (filedes=open("/dev/urandom", O_RDONLY)) > 0 ) {
      ssize_t res = read(filedes, (void *)&seed, sizeof(unsigned long int));
      if ( !res ) {
	perror("Unable to seed random number generator!\n");
	exit(errno);
      }
      close(filedes);
    } else {
      perror("Unable to seed random number generator!\n");
      exit(errno);
    }
    // Dump out the seed in case we'd like to do this exact run again
    fprintf(stderr, "Starting run with random seed: %lu\n", seed);
  }

  initstate( seed, (char *)state, 256);
  srandom(seed);

  return;
}
***/
