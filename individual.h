#ifndef __INDIVIDUAL_H
#define __INDIVIDUAL_H

#include "global.h"

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <limits>
#include <math.h>

using namespace std;

class individual {

 public:
  individual( void );
  individual( individual **, bool=true, bool=false );
  ~individual( void );

  void testFitness( void );
  void mutate_simple( void );
  void mutate( void );
  void copy ( individual *, bool=false );
  void set_genes( void );
  void output( bool=false );
  bool isClone( individual * );

  individual *make_baby(individual *);
  individual *get_mate( int, individual ** );

  int count;
  int nGenes;
  float fitness;
  float *gene;
  int progeny;
  int generation;
  float max_fitness;
  int accuracy;

  individual *previous;
  individual *next;

  individual *operator* (individual & );
  individual *operator= (individual & );

  //ostream &operator<< ( individual * );

 protected:
 
 private:

};

#endif
