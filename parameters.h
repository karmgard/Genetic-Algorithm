#ifndef __PARAMETERS_H
#define __PARAMETERS_H

#include <iostream>
#include <string>
#include <unordered_map>

#include <fcntl.h>
#include <errno.h>
#include <string.h>

#define MAX_UL_INT 0xffffffff
#define MAX_INT    0x7fffffff

using namespace std;

union Values {

  /*
  char asChar;
  unsigned char asUChar;
  short asShort;
  unsigned short asUShort;
  */

  long asLong;
  unsigned long asULong; 
  bool asBool;
  int asInt;
  float asFloat;
  double asDouble;
  char* asCString;
  unsigned int asUInt;

  Values() { asULong = 0; }
  /*Values(char in) { asUChar = in; }
  Values(unsigned char in) { asChar = in; }
  Values(short in) { asShort = in; }
  Values(unsigned short in) { asUShort = in; }
  */
  Values(long in) { asLong = in; }
  Values(unsigned int in) { asUInt = in; }
  Values(unsigned long in) { asULong = in; }
  Values(bool in) { asBool = in; }
  Values(int in) { asInt = in; }
  Values(float in) { asFloat = in; }
  Values(double in) { asDouble = in; }
  Values(char* in) { asCString = in; }

  /*operator char() { return asChar; }
  operator unsigned char() { return asUChar; }
  operator short() { return asShort; }
  operator unsigned short() { return asUShort; }
  */
  operator long() { return asLong; }
  operator unsigned int() { return asUInt; }

  operator unsigned long() { return asULong; }
  operator bool() { return asBool; }
  operator int() { return asInt; }
  operator float() { return asFloat; }
  operator double() { return asDouble; }
  operator char*() { return asCString; }
};

typedef unordered_map <string, Values> paramMap;

class parameters {

 public:

  parameters( char * );
  ~parameters( void );

  ulong  getULong   ( string );
  bool   getBool    ( string );
  int    getInt     ( string );
  float  getFloat   ( string );
  double getDouble  ( string );
  char * getString  ( string );
  uint   getLong    ( string );
  uint   getUInt    ( string );

  /*
  char   getChar    ( string );
  short  getShort   ( string );
  ushort getUShort  ( string );
  */

  double *pLO;
  double *pHI;

  // Internal values 
  int VERBOSE;
  bool SHOW_PLOT;
  int PLOT_FREQ;
  float EXIT_LIMIT;
  int CPU_USAGE_LIMIT;
  int MAXIMUM_GENERATIONS;
  int ACCURACY;
  int DUMP_N_TOP;

  int ELITISM_GENERATIONS;
  float PERCENT_ELITES_KEPT;
  bool KEEP_STABLE_POPULATION;
  float MAX_FITNESS;
  int INITIAL_POPULATION;
  int NUMBER_OF_GENES;
  float MUTATION_RATE;
  float MUTATION_GAIN;
  string SORT_TYPE;
  uint NUM_THREADS;
  bool MUTATE_SIMPLE;
  unsigned long SEED;

 protected:

 private:
  
  paramMap pMap;

};

#endif
