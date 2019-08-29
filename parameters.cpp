#include "parameters.h"

parameters::parameters( char * filename ) {

  FILE *pFile;
  char inStream[256];

  /*-- Now, try to read in the parameter file --*/
  errno = 0;
  if ( (pFile=fopen(filename,"r")) == NULL) {
    perror(filename);
    exit(errno);
  }

  char key[256], value[256], type[32];

  memset( (void *)key,   0, 256*sizeof(char));
  memset( (void *)value, 0, 256*sizeof(char));
  memset( (void *)type,  0,   32*sizeof(char));

  while ( fgets(inStream, 255, pFile) != NULL ) {

    if ( inStream[0] == '#' )
      continue;

    if ( !strncmp(inStream, "EOF", 3) )
      break;

    sscanf(inStream, "%s %s = %s %*s", type, key, value);

    if ( !strncmp(type, "EOF", 3) )
      break;

    if ( !type[0] || !key[0] || !value[0] || type[0] == '#' || type[0] == ' ' )
      continue;

    if ( !strcmp(type, "int") ) {
      if ( strstr( (const char *)inStream, "inf" ) || strstr( (const char *)inStream, "INF" ) )
	pMap[key] = MAX_INT;
      else
	pMap[key] = atoi(value);
    }

    if ( !strcmp(type, "long") ) {
      if ( strstr( (const char *)inStream, "inf" ) || strstr( (const char *)inStream, "INF" ) )
	pMap[key] = MAX_INT;
      else
	pMap[key] = atol(value);
    }

    if ( !strcmp(type, "ulong") ) {
      if ( strstr( (const char *)inStream, "inf" ) || strstr( (const char *)inStream, "INF" ) )
	pMap[key] = MAX_INT;
      else
	pMap[key] = abs(atol(value));
    }

    else if ( !strcmp(type, "uint") ) {
      if ( strstr( (const char *)inStream, "inf" ) || strstr( (const char *)inStream, "INF" ) )
	pMap[key] = MAX_INT;
      else
	pMap[key] = abs(atoi(value));
    }

    else if ( !strcmp(type, "string") ) {
      pMap[key] = new char [strlen(value)+2];
      strcpy(pMap[key], value);
    }

    else if ( !strcmp(type, "float") || !strcmp(type, "double") )
      pMap[key] = strtod(value, (char **)NULL);

    else if ( !strcmp(type, "bool") )
      pMap[key] = (!strcmp(value, "true")) ? true : false;

    else if ( !strcmp(type, "float_array") ) {

      char temp[11 + 4*pMap["NUMBER_OF_GENES"].asInt];
      inStream[strlen(inStream)-1] = '\0';
      double number = 0.0f;

      for (int i=0; i<pMap["NUMBER_OF_GENES"].asInt; i++) {
	strcpy( temp, "%*s %*s = " );

	for (int j=0; j<i; j++)
	  strcat(temp, "%*f, ");
	strcat(temp, "%lf\0");
 
	if ( !strcmp(key, "LOWER_LIMITS") ) {
	  sscanf(inStream, temp, &number);
	  pLO[i] = number;
	}
	else if ( !strcmp(key, "UPPER_LIMITS") ) {
	  sscanf(inStream, temp, &number);
	  pHI[i] = (double)number;
	}
      }
    }

    if ( !strcmp(key, "NUMBER_OF_GENES") ) {
      pLO = new double[pMap["NUMBER_OF_GENES"].asInt];
      pHI = new double[pMap["NUMBER_OF_GENES"].asInt];
    }

    memset( (void *)key,   0, 128*sizeof(char));
    memset( (void *)value, 0, 128*sizeof(char));
    memset( (void *)type,  0,   16*sizeof(char));
  }

  if ( pMap.count("LOWER_LIMIT_ALL") ) {
    for ( int i=0; i<pMap["NUMBER_OF_GENES"].asInt; i++ )
      pLO[i] = pMap["LOWER_LIMIT_ALL"].asDouble;
  }
  if ( pMap.count("UPPER_LIMIT_ALL") ) {
    for ( int i=0; i<pMap["NUMBER_OF_GENES"].asInt; i++ )
      pHI[i] = pMap["UPPER_LIMIT_ALL"].asDouble;
  }

  // Hard wire all the internal stuff we know we'll need so that access is faster
  // User defined parameters for the fitness function should be accessed through 
  // the accessors to the map
  VERBOSE                = getInt("VERBOSE");
  SHOW_PLOT              = getBool("SHOW_PLOT");
  PLOT_FREQ              = getInt("PLOT_FREQ");
  EXIT_LIMIT             = getFloat("EXIT_LIMIT");
  CPU_USAGE_LIMIT        = getInt("CPU_USAGE_LIMIT");
  MAXIMUM_GENERATIONS    = getInt("MAXIMUM_GENERATIONS");
  ACCURACY               = getInt("ACCURACY");
  DUMP_N_TOP             = getInt("DUMP_N_TOP");
  ELITISM_GENERATIONS    = getInt("ELITISM_GENERATIONS");
  PERCENT_ELITES_KEPT    = getFloat("PERCENT_ELITES_KEPT");
  KEEP_STABLE_POPULATION = getBool("KEEP_STABLE_POPULATION");
  MAX_FITNESS            = getFloat("MAX_FITNESS");
  INITIAL_POPULATION     = getInt("INITIAL_POPULATION");
  NUMBER_OF_GENES        = getInt("NUMBER_OF_GENES");
  MUTATION_RATE          = getFloat("MUTATION_RATE");
  MUTATION_GAIN          = getFloat("MUTATION_GAIN");
  SORT_TYPE              = getString("SORT_TYPE");
  NUM_THREADS            = getUInt("NUM_THREADS");
  SEED                   = getULong("SEED");

  return;
}

parameters::~parameters( void ) {
  delete [] pLO;
  delete [] pHI;
  delete [] pMap["FITNESS_FUNCTION"];
  return;
}

unsigned long parameters::getULong( string key ) {
    if ( this->pMap.count(key) )
    return this->pMap[key].asULong;
  else
    return 0;
}

bool parameters::getBool( string key ) {
    if ( this->pMap.count(key) )
      return this->pMap[key].asBool;
    else
      return false;
}

int parameters::getInt( string key ) {
    if ( this->pMap.count(key) )
      return this->pMap[key].asInt;
    else
      return 0;
}

float parameters::getFloat( string key ) {
  if ( this->pMap.count(key) )
    return (float)this->pMap[key].asDouble;
  else
    return 0.0;
}

double parameters::getDouble( string key ) {
  if ( this->pMap.count(key) )
    return this->pMap[key].asDouble;
  else
    return 0.0;
}

char * parameters::getString( string key ) {
      if ( this->pMap.count(key) )
    return this->pMap[key].asCString;
  else
    return '\0';
}

unsigned int parameters::getLong( string key ) {
    if ( this->pMap.count(key) )
    return this->pMap[key].asLong;
  else
    return 0;
}

unsigned int parameters::getUInt( string key ) {
  if ( this->pMap.count(key) )
    return this->pMap[key].asUInt;
  else
    return 0;
}

/*

char parameters::getChar( string key ) {
    if ( this->pMap.count(key) )
    return this->pMap[key].asChar;
  else
    return '\0';
}

short parameters::getShort( string key ) {
    if ( this->pMap.count(key) )
    return this->pMap[key].asShort;
  else
    return 0;
}

unsigned short parameters::getUShort( string key ) {
    if ( this->pMap.count(key) )
    return this->pMap[key].asUShort;
  else
    return 0;
}

*/

