#ifndef __POPULATION_H
#define __POPULATION_H

#include "individual.h"
#include "global.h"
#include "parameters.h"
#include "utilities.h"
#include "fitness.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <string.h>

class population {

 public:
  population( void );
  ~population( void );

  float get_avg_fitness( void );
  float get_stdev_fitness( void );
  int   get_count( void );
  void  dump( int max = 0 );
  void  mate( void );
  unsigned int generation;
  double *get_population_fitness(void);
  void print( void );

  unsigned int clones;
  unsigned int count;
  double stdev;
  double average;
  double variation;

  individual *mostfit;

 protected:

 private:
  void copy_elites( void );
  void roulette_fill( void );
  void check_for_clones( void );
  void recount( void );
  void copy( population * );
  void flush( void );
  void trim ( int );
  void get_fittest( void );
  void mutation_gain( void );
  void get_statistics( void );

  void sort( void );
  void heap_sort( void ** );
  void quick_sort( void **, uint );
  void quick_sort( void ** );

  double *fitness_array;
  unsigned int allocation;
  bool mating_in_progress;

  individual * get_individual( int );
  individual * push( void );
  individual * push( individual * );
  individual * pop( void );
  individual * shift( void );
  individual * unshift( individual * );
  individual * unshift( void );

  individual *first;
  individual *last;
  individual *person;

};

#endif
