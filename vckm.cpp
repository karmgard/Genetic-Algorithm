#include "vckm.h"

static bool  proc_matrix_element[9] = {true, true, true, true, true, true, true, true, true};

/*-- Prototype the parameterized matrices functions --*/
static void complexroots( complex *M, 
			  float a, float b, float c, 
			  int sign_a, int sign_b, int sign_c, 
			  float d1, float d2 );
static void calculate_general_matrix( complex *m, float Cu, float Cd );
static void calculate_D_matrix( complex *m, 
				float D1u, float D1d, 
				float D2u, float D2d,
				float o1,  float o2 );
static void calculate_A_matrix(complex *m, float Au, float Ad, float d1, float d2);
static void calculate_Y_matrix(complex *m, float Au, float Ad, float d1, float d2);

#ifndef Mu
static float Mu, Md, Ms, Mc, Mb, Mt;

void Vckm(char which, complex *m, float C1u, float C1d, float C2u, float C2d, float d1, float d2,
	  float mu, float md, float ms, float mc, float mb, float mt) {

  /*-- Set the quark masses for this try --*/
  Mu = mu;
  Md = md;
  Ms = ms;
  Mc = mc;
  Mb = mb;
  Mt = mt;
#else
  void Vckm(char which, complex *m, float C1u, float C1d, float C2u, float C2d, float d1, float d2) {
#endif

    switch (which) {
    case 'D':
      calculate_D_matrix(m, C1u, C1d, C2u, C2d, d1, d2);
      break;
    case 'A':
      calculate_A_matrix(m, C1u, C1d, d1, d2);
      break;
    case 'Y':
      calculate_Y_matrix(m, C1u, C1d, d1, d2);
      break;
    case 'G':
      calculate_general_matrix( m, C1u, C1d );
      break;
    }
    return;
  }

  static void complexroots( complex *M, 
			    float a, float b, float c, 
			    int sign_a, int sign_b, int sign_c, 
			    float d1, float d2 ) {
    if ( a < 0 ) {
      M->r = 0.0f;
      M->i = +sign_a*sqrt(-a);
    } else {
      M->r = +sign_a*sqrt(+a);
      M->i = 0.0f;
    }
    if ( b < 0 ) {
      M->r += -sign_b*sin(d1)*sqrt(-b);
      M->i += +sign_b*cos(d1)*sqrt(-b);
    } else {
      M->r += +sign_b*cos(d1)*sqrt(+b);
      M->i += +sign_b*sin(d1)*sqrt(+b);
    }
    if ( c < 0 ) {
      M->r += -sign_c*sin(d2)*sqrt(-c);
      M->i += +sign_c*cos(d2)*sqrt(-c);
    } else {
      M->r += +sign_c*cos(d2)*sqrt(+c);
      M->i += +sign_c*sin(d2)*sqrt(+c);
    }
  }

  /*-- Print the matrix resulting from the specified parameters --*/
  void dumpMatrix( void *p ) {
    complex *M = (complex *)calloc((size_t) 9, sizeof(complex));
    float V[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    individual *person = (individual *)p;
    
    if ( !VCKM_MATRIX )
      VCKM_MATRIX = params->getString("FITNESS_FUNCTION")[3];

    printf("\n\nThe Point (");
    for (int i=0; i<person->nGenes; i++) {
      if ( i < person->nGenes-1 )
	printf("%f, ", person->gene[i]);
      else
	printf("%f) has fitness %f\n ", person->gene[i], person->fitness);
    }

    errno = 0;
    switch (VCKM_MATRIX) {
    case 'A': case 'a': case 'Y': case 'y':
      Vckm( VCKM_MATRIX, M, 
	    person->gene[0], person->gene[1], 
	    0.0f,            0.0f, 
	    person->gene[2], person->gene[3] );
      break;
    case 'd': case 'D':
      Vckm( VCKM_MATRIX, M, 
	    person->gene[0], person->gene[1], person->gene[2], 
	    person->gene[3], person->gene[4], person->gene[5] );
      break;
    default:
      printf("That matrix is not yet defined! Bailing...\n");
      exit(0);
    }

    if ( errno ) {
      perror("vckm()");
      printf("Failed with error %i\n", errno);
      exit(errno);
    }

    errno = 0;
    for (int i=0; i<9; i++)
      V[i] = Cabs(M[i]);
    
    if ( errno ) {
      perror("Cabs(M)");
      printf("Failed with error %i\n", errno);
      exit(errno);
    }

    /*-- Free up the memory allocated to the matrix --*/
    free(M);

    printf("Which results in\n");
    printf( "         [ %1.4f %1.3f  %1.4f ]\n",  V[0],V[1],V[2]);
    printf( "|Vckm| = [ %1.3f  %1.4f %1.3f  ]\n", V[3],V[4],V[5]);
    printf( "         [ %1.3f  %1.3f  %1.4f ]\n", V[6],V[7],V[8]);

    if ( person->fitness > 0.0f ) {

      float errors[9] = {0.0f};
      for ( int i=0; i<9; i++ ) {
	if ( V[i] < V_Min[i] )
	  errors[i] = V[i] - V_Min[i];
	else if ( V[i] > V_Max[i] )
	  errors[i] = V[i] - V_Max[i];
      }

      printf("and has errors\n");
      printf( "         [ %+1.4f %+1.3f  %+1.4f ]\n",    errors[0],errors[1],errors[2]);
      printf( "|Vckm| = [ %+1.3f  %+1.4f %+1.3f  ]\n",   errors[3],errors[4],errors[5]);
      printf( "         [ %+1.3f  %+1.3f  %+1.4f ]\n\n", errors[6],errors[7],errors[8]);
    }
    return;
  }

  /*-- Entry point into this file --*/
  void ckm_fitness( void *p ) {
    if ( !VCKM_MATRIX )
      VCKM_MATRIX = params->getString("FITNESS_FUNCTION")[3];

    if ( !VCKM_MAX_FITNESS )
      VCKM_MAX_FITNESS = params->getFloat("MAX_FITNESS");

    if ( !VCKM_EXIT_LIMIT )
      VCKM_EXIT_LIMIT = params->getFloat("EXIT_LIMIT");

    complex *M = (complex *)calloc((size_t) 9, sizeof(complex));
    float V[9];
    float test = 0.0f, result = VCKM_MAX_FITNESS;

    individual *person = (individual *)p;

    for ( int i=0; i<9; i++)
      V[i] = 10.0f;

    errno = 0;
    switch (VCKM_MATRIX) {
    case 'A': case 'a': case 'Y': case 'y':
      Vckm( VCKM_MATRIX, M, 
	    person->gene[0], person->gene[1], 
	    0.0f,            0.0f, 
	    person->gene[2], person->gene[3] );
      break;
    case 'd': case 'D':
      Vckm( VCKM_MATRIX, M, 
	    person->gene[0], person->gene[1], person->gene[2], 
	    person->gene[3], person->gene[4], person->gene[5] );
      break;
    default:
      printf("That matrix is not yet defined! Bailing...\n");
      exit(1);
    }

    if ( errno ) {
      perror("vckm()");
      printf("Failed with error %i\n", errno);
      exit(errno);
    }

    errno = 0;
    for (int i=0; i<9; i++)
      V[i] = Cabs(M[i]);
    
    if ( errno ) {
      perror("Cabs(M)");
      printf("Failed with error %i\n", errno);
      exit(errno);
    }

    /*-- Free up the memory allocated to the matrix --*/
    free(M);

    /*-- Now, see how well this matrix works.... --*/
    for ( int i=0; i<9; i++ ) {
      if ( proc_matrix_element[i] ) {
	if ( proc_matrix_element[i] != proc_matrix_element[i] || fabs(proc_matrix_element[i]) == INFINITY )
	  test =  VCKM_MAX_FITNESS;
	else if ( V[i] > 2.0000 )
	  test = VCKM_MAX_FITNESS;

	if ( i == 0 || i == 2 || i == 4 || i == 8 ) {
	  if ( fround(V[i], 4) > V_Max[i] || fround(V[i], 4) < V_Min[i] )
	    test += 
	      (V[i] - ( 0.5*(V_Max[i]+V_Min[i]) ))*
	      (V[i] - ( 0.5*(V_Max[i]+V_Min[i]) ));
	} else {
	  if ( fround(V[i], 3) > V_Max[i] || fround(V[i], 3) < V_Min[i] )
	    test += 
	      (V[i] - ( 0.5*(V_Max[i]+V_Min[i]) ))*
	      (V[i] - ( 0.5*(V_Max[i]+V_Min[i]) ));
	}
      }
    }

    if ( test >= 0 )
      result = sqrt(test);
    else
      perror("Sqrt of -1!\n");

    if ( result > VCKM_MAX_FITNESS )
      result = VCKM_MAX_FITNESS;

    person->fitness = result;

    // If this is a no-exit run.... dump out the results if they're perfect
    if ( VCKM_EXIT_LIMIT < 0.0f && result == 0.0f )
      dumpMatrix( person );

    return;
  }

  static void calculate_general_matrix( complex *m, float Cu, float Cd ) {

    printf("\n\n==============> Error! Not set up yet. Bailing <====================\n\n");
    exit(1);

    return;
  }

  static void calculate_D_matrix( complex *m, 
				  float D1u, float D1d, 
				  float D2u, float D2d,
				  float o1,  float o2 ) {
    float a, b, c;

    //o1 = 80.1792;
    //o2 = 174.9;
    //o2 = 89.468529 

    if ( proc_matrix_element[0] ) {

      /*-- Calculate M(1,1) --*/

      a = ((-D1d + Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)*(D1u + Mc)*
	   (D1u + D2u + Mc - Mt)*(-D1u + Mt))/
	((Mb - Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)*
	 (2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)*(Mc + Mu));
      
      b = ((-D1d + Md)*(-D1u + Mu)*(-D1d + Md)*(-D1u + Mu))*
	((2*D1d + D2d - Mb - Md + Ms)/(-D1d + Mb)*(D1d - Md)*(D1d + Ms))*
	(((D1d - Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms))/
	 ((-Mb + Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)))*
	(((2*D1u + D2u + Mc - Mt - Mu))/((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)))*
	(((2*D1u + D2u + Mc - Mt - Mu)*(Mc + Mu)*(-Mt + Mu))/
	 ((D1u + Mc)*(D1u - Mt)*(D1u + D2u + Mc - Mt)));
      
      c = ((-D1d + Md)*(-D1d + Md)*(-D1u + Mu)*(-D1u + Mu)*
	   ((-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	    (D1d + D2d - Md + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	   ((D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mu)*(-D1u - D2u + Mt + Mu)/
	    (2*D1u + D2u + Mc - Mt - Mu)) ) / 
	( (D1d + D2d - Mb + Ms)*(D1d + D2d - Mb + Ms)*(D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mt)*
	  ((-D1d + Mb)*(D1d - Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	  ((-Mb + Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)/
	   ((D1d - Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)))*
	  ((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)/(2*D1u + D2u + Mc - Mt - Mu))*
	  ((2*D1u + D2u + Mc - Mt - Mu)*(Mc + Mu)*(-Mt + Mu)/
	   ((D1u + Mc)*(D1u - Mt)*(D1u + D2u + Mc - Mt))) );
      
      complexroots(&m[0], a, b, c, +1, +1, +1, o1, o2);
    }
    
    if ( proc_matrix_element[1] ) {
      
      /*-- Calculate M(1,2) --*/
      
      a = ((-D1d + Mb)*(D1d - Md)*(D1d + Ms)*(-D1d - D2d + Mb + Md)*(D1d + Ms)*
	   (D1u + Mc)*(D1u + D2u + Mc - Mt)*(-D1u + Mt)) /
	((2*D1d + D2d - Mb - Md + Ms)*(Mb + Ms)*(Md + Ms)*(2*D1u + D2u + Mc - Mt - Mu)*
	 (Mt - Mu)*(Mc + Mu)*(D1d + Ms)*(D1d + Ms));
      
      b = ( (-D1u + Mu)*(-D1u + Mu)*(-D1d - D2d + Mb + Md)*(D1d + Ms)*
	    (2*D1u + D2u + Mc - Mt - Mu)*(D1u + Mc)*(D1u - Mt)*(D1u + D2u + Mc - Mt) ) /
	( (Mb + Ms)*(Md + Ms)*(D1u + Mc)*(-D1u + Mt)*(D1u - Mu)*
	  (2*D1u + D2u + Mc - Mt - Mu)*(Mc + Mu)*(-Mt + Mu) );
      
      c = ( (-D1u + Mu)*(-D1u + Mu)*(D1d + D2d - Mb - Md)*(-D1d - D2d + Md - Ms)*
	    (D1d + D2d - Mb + Ms)/(2*D1d + D2d - Mb - Md + Ms)*(D1u + D2u + Mc - Mt)*
	    (D1u + D2u + Mc - Mu)*(-D1u - D2u + Mt + Mu) ) /
	( (2*D1u + D2u + Mc - Mt - Mu)*(2*D1u + D2u + Mc - Mt - Mu)*(D1u + D2u + Mc - Mt)*
	  (D1u + D2u + Mc - Mt)*((Mb + Ms)*(Md + Ms)/((-D1d - D2d + Mb + Md)*(D1d + Ms)))*
	  ((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)/(2*D1u + D2u + Mc - Mt - Mu))*
	  ((2*D1u + D2u + Mc - Mt - Mu)*(Mc + Mu)*(-Mt + Mu)/
	   ((D1u + Mc)*(D1u - Mt)*(D1u + D2u + Mc - Mt))) );
      
      complexroots(&m[1], a, b, c, -1, +1, +1, o1, o2);
    }
    
    if ( proc_matrix_element[2] ) {
      
      /*-- Calculate M(1,3) --*/
      
      a = ( (2*D1d + D2d - Mb - Md + Ms)*(2*D1d + D2d - Mb - Md + Ms)*
	    ((D1d - Mb)*(-D1d + Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	    ((-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*(D1d + D2d - Md + Ms)/
	     (2*D1d + D2d - Mb - Md + Ms)) ) /
	( (-D1d + Mb)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  (-D1d + Mb)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)*(Mc + Mu)/
	   ((D1u + Mc)*(D1u + D2u + Mc - Mt)*(-D1u + Mt))) );
      
      b = ( (2*D1d + D2d - Mb - Md + Ms)*(-D1u + Mu)*(2*D1d + D2d - Mb - Md + Ms)*(-D1u + Mu)*
	    ((D1d + D2d - Mb - Md)*(-D1d - D2d + Md - Ms)*(D1d + D2d - Mb + Ms)/
	     (2*D1d + D2d - Mb - Md + Ms)) )/
	( (-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)/(2*D1u + D2u + Mc - Mt - Mu))*
	  ((2*D1u + D2u + Mc - Mt - Mu)*(Mc + Mu)*(-Mt + Mu)/
	   ((D1u + Mc)*(D1u - Mt)*(D1u + D2u + Mc - Mt))) );

      c = ( (-D1u + Mu)*(-D1u + Mu)*
	    ((D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mu)*
	     (-D1u - D2u + Mt + Mu)/(2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mt)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)/(2*D1u + D2u + Mc - Mt - Mu))*
	  ((2*D1u + D2u + Mc - Mt - Mu)*(Mc + Mu)*(-Mt + Mu)/
	   ((D1u + Mc)*(D1u - Mt)*(D1u + D2u + Mc - Mt))) );
      
      complexroots(&m[2], a, b, c, +1, +1, +1, o1, o2);
    }

    if ( proc_matrix_element[3] ) {
      
      /*-- Calculate M(2,1) --*/
      
      a = ((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)/(2*D1u + D2u + Mc - Mt - Mu))/
	((D1u + Mc)*(D1u + Mc)*((Mb - Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)/
				((-D1d + Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)))*
	 ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))));
      
      b = ((-D1d + Md)*(-D1d + Md)) /
	(((-D1d + Mb)*(D1d - Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	 ((-Mb + Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)/
	  ((D1d - Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)))*
	 ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))));
      
      c = ((-D1d + Md)*(-D1d + Md)*((-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
				    (D1d + D2d - Md + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	   ((D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu)*(-D1u - D2u - Mc + Mu)/
	    (2*D1u + D2u + Mc - Mt - Mu)))/
	((D1d + D2d - Mb + Ms)*(D1d + D2d - Mb + Ms)*
	 ((-D1d + Mb)*(D1d - Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	 ((-Mb + Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)/
	  ((D1d - Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)))*
	 (D1u + D2u - Mt - Mu)*(D1u + D2u - Mt - Mu)*
	 ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))));
      
      complexroots(&m[3], a, b, c, -1, +1, +1, o1, o2);
    }
    
    if ( proc_matrix_element[4] ) {
      
      /*-- Calculate M(2,2) --*/
      a = (((-D1d + Mb)*(D1d - Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	   ((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)/(2*D1u + D2u + Mc - Mt - Mu))) /
	((D1u + Mc)*(D1d + Ms)*(D1u + Mc)*(D1d + Ms)*
	 ((Mb + Ms)*(Md + Ms)/((-D1d - D2d + Mb + Md)*(D1d + Ms)))*
	 ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))));
      
      b = ((-D1d - D2d + Mb + Md)*(D1d + Ms)*(D1u + Mc)*(-D1u - D2u + Mt + Mu)) / 
	((Mb + Ms)*(Md + Ms)*(Mc + Mt)*(Mc + Mu));
      
      c = (((D1d + D2d - Mb - Md)*(-D1d - D2d + Md - Ms)*(D1d + D2d - Mb + Ms)/
	    (2*D1d + D2d - Mb - Md + Ms))*((D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu)*
					   (-D1u - D2u - Mc + Mu)/(2*D1u + D2u + Mc - Mt - Mu))) /
	((D1d + D2d - Mb - Md)*(D1d + D2d - Mb - Md)*(D1u + D2u - Mt - Mu)*(D1u + D2u - Mt - Mu)*
	 ((Mb + Ms)*(Md + Ms)/((-D1d - D2d + Mb + Md)*(D1d + Ms)))*
	 ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))));
      
      complexroots(&m[4], a, b, c, +1, +1, +1, o1, o2);
    }
    
    if ( proc_matrix_element[5] ) {
      
      /*-- Calculate M(2,3) --*/
      
      a = ( (2*D1d + D2d - Mb - Md + Ms)*(2*D1d + D2d - Mb - Md + Ms)*
	    ((D1d - Mb)*(-D1d + Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	    ((-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*(D1d + D2d - Md + Ms)/
	     (2*D1d + D2d - Mb - Md + Ms))*((D1u + Mc)*(-D1u + Mt)*(D1u - Mu)/
					    (2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (-D1d + Mb)*(D1u + Mc)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  (-D1d + Mb)*(D1u + Mc)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))) );
      
      b = ( (2*D1d + D2d - Mb - Md + Ms)*(2*D1d + D2d - Mb - Md + Ms)*
	    ((D1d + D2d - Mb - Md)*(-D1d - D2d + Md - Ms)*(D1d + D2d - Mb + Ms)/
	     (2*D1d + D2d - Mb - Md + Ms)) ) /
	( (-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))) );
      
      c = ( (D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu)*(-D1u - D2u - Mc + Mu)/
	    (2*D1u + D2u + Mc - Mt - Mu) ) /
	( (D1u + D2u - Mt - Mu)*(D1u + D2u - Mt - Mu)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(Mc + Mu)/((D1u + Mc)*(-D1u - D2u + Mt + Mu))) );
      
      complexroots(&m[5], a, b, c, -1, +1, +1, o1, o2);
    }
    
    if ( proc_matrix_element[6] ) {
      
      /*-- Calculate M(3,1) --*/
      
      a = ( (2*D1u + D2u + Mc - Mt - Mu)*(2*D1u + D2u + Mc - Mt - Mu)*
	    ((D1u + Mc)*(D1u - Mt)*(-D1u + Mu)/(2*D1u + D2u + Mc - Mt - Mu))*
	    ((D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mu)*(-D1u - D2u + Mt + Mu)/
	     (2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (D1u + D2u + Mc - Mt)*(-D1u + Mt)*(D1u + D2u + Mc - Mt)*(-D1u + Mt)*
	  (-D1u - D2u + Mt + Mu)*(-D1u - D2u + Mt + Mu)*
	  ((Mb - Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)/
	   ((-D1d + Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))) );

      b = ( (-D1d + Md)*(2*D1u + D2u + Mc - Mt - Mu)*(-D1d + Md)*(2*D1u + D2u + Mc - Mt - Mu)*
	    ((D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu)*(-D1u - D2u - Mc + Mu)/
	     (2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mt)*(-D1u - D2u + Mt + Mu)*(-D1u - D2u + Mt + Mu)*
	  ((-D1d + Mb)*(D1d - Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	  ((-Mb + Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)/
	   ((D1d - Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))) );
      
      c = ( (-D1d + Md)*(-D1d + Md)*((-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
				     (D1d + D2d - Md + Ms)/(2*D1d + D2d - Mb - Md + Ms)) ) /
	( (D1d + D2d - Mb + Ms)*(D1d + D2d - Mb + Ms)*
	  ((-D1d + Mb)*(D1d - Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	  ((-Mb + Md)*(2*D1d + D2d - Mb - Md + Ms)*(Md + Ms)/
	   ((D1d - Mb)*(D1d + Ms)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))));
      
      complexroots(&m[6], a, b, c, +1, +1, +1, o1, o2);
    }
    
    if ( proc_matrix_element[7] ) {
      
      /*-- Compute M(3,2) --*/
      a = ( (2*D1u + D2u + Mc - Mt - Mu)*(2*D1u + D2u + Mc - Mt - Mu)*
	    ((-D1d + Mb)*(D1d - Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	    ((D1u + Mc)*(D1u - Mt)*(-D1u + Mu)/(2*D1u + D2u + Mc - Mt - Mu))*
	    ((D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mu)*(-D1u - D2u + Mt + Mu)/
	     (2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (D1d + Ms)*(D1u + D2u + Mc - Mt)*(-D1u + Mt)*(-D1u - D2u + Mt + Mu)*
	  (D1d + Ms)*(D1u + D2u + Mc - Mt)*(-D1u + Mt)*(-D1u - D2u + Mt + Mu)*
	  ((Mb + Ms)*(Md + Ms)/((-D1d - D2d + Mb + Md)*(D1d + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))));
      
      b = ( (2*D1u + D2u + Mc - Mt - Mu)*(2*D1u + D2u + Mc - Mt - Mu)*
	    ((D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu)*(-D1u - D2u - Mc + Mu)/
	     (2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mt)*(-D1u - D2u + Mt + Mu)*(-D1u - D2u + Mt + Mu)*
	  ((Mb + Ms)*(Md + Ms)/((-D1d - D2d + Mb + Md)*(D1d + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))) ); 
      
      c = ( (D1d + D2d - Mb - Md)*(-D1d - D2d + Md - Ms)*(D1d + D2d - Mb + Ms)/
	    (2*D1d + D2d - Mb - Md + Ms) ) /
	( (D1d + D2d - Mb - Md)*(D1d + D2d - Mb - Md)*
	  ((Mb + Ms)*(Md + Ms)/((-D1d - D2d + Mb + Md)*(D1d + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))) );
      
      complexroots(&m[7], a, b, c, -1, +1, +1, o1, o2);
    }
    
    if ( proc_matrix_element[8] ) {
      
      /*-- Calculate M(3,3) --*/
      
      a = ( (2*D1d + D2d - Mb - Md + Ms)*(2*D1u + D2u + Mc - Mt - Mu)*
	    (2*D1d + D2d - Mb - Md + Ms)*(2*D1u + D2u + Mc - Mt - Mu)*
	    ((D1d - Mb)*(-D1d + Md)*(D1d + Ms)/(2*D1d + D2d - Mb - Md + Ms))*
	    ((-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*(D1d + D2d - Md + Ms)/
	     (2*D1d + D2d - Mb - Md + Ms))*
	    ((D1u + Mc)*(D1u - Mt)*(-D1u + Mu)/(2*D1u + D2u + Mc - Mt - Mu))*
	    ((D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mu)*(-D1u - D2u + Mt + Mu)/
	     (2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (-D1d + Mb)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  (-D1d + Mb)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  (D1u + D2u + Mc - Mt)*(-D1u + Mt)*(D1u + D2u + Mc - Mt)*(-D1u + Mt)*
	  (-D1u - D2u + Mt + Mu)*(-D1u - D2u + Mt + Mu)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))) );
      
      b = ( (2*D1d + D2d - Mb - Md + Ms)*(2*D1u + D2u + Mc - Mt - Mu)*
	    (2*D1d + D2d - Mb - Md + Ms)*(2*D1u + D2u + Mc - Mt - Mu)*
	    ((D1d + D2d - Mb - Md)*(-D1d - D2d + Md - Ms)*(D1d + D2d - Mb + Ms)/
	     (2*D1d + D2d - Mb - Md + Ms))*
	    ((D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu)*(-D1u - D2u - Mc + Mu)/
	     (2*D1u + D2u + Mc - Mt - Mu)) ) /
	( (-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*(-D1d - D2d + Mb + Md)*(D1d + D2d - Mb + Ms)*
	  (D1u + D2u + Mc - Mt)*(D1u + D2u + Mc - Mt)*(-D1u - D2u + Mt + Mu)*(-D1u - D2u + Mt + Mu)*
	  ((Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)/
	   ((D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)))*
	  ((Mc + Mt)*(2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu)/
	   ((D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu))) );
      
      c = ( (D1d - Mb)*(D1d + D2d - Mb - Md)*(D1d + D2d - Mb + Ms)*
	    (D1u - Mt)*(D1u + D2u + Mc - Mt)*(D1u + D2u - Mt - Mu) ) /
	( (Mb - Md)*(Mb + Ms)*(2*D1d + D2d - Mb - Md + Ms)*(Mc + Mt)*
	  (2*D1u + D2u + Mc - Mt - Mu)*(Mt - Mu) );
      
      complexroots(&m[8], a, b, c, +1, +1, +1, o1, o2);
    }
    
    return;
  }

  static void calculate_A_matrix(complex *m, float Au, float Ad, float d1, float d2) {
    float a, b, c;

    if ( proc_matrix_element[0] ) {

      /*-- Calculate the various terms in the matrix element M11 --*/
      a = Ms*Mc*Mt*Mb*(-Ad + Mb - Ms)*(-Au - Mc + Mt)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Md*Mu*(-Ad + Mb - Ms)*(-Au - Mc + Mt)/(((Mb - Md)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)));
      
      c = Mu*Md*(-Ad + Mb + Md)*(Ad - Md + Ms)*(Au + Mc - Mu)*(-Au + Mt + Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[0], a, b, c, +1, +1, +1, d1, d2);
    }
    
    if ( proc_matrix_element[1] ) {
      
      /*-- Calculate the various terms in the matrix element M12 --*/
      a = Mc*Mt*Mb*Md*(-Ad + Mb + Md)*(-Au - Mc + Mt)/
	(((-Ad + Mb + Md - Ms)*(Mb + Ms)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Mu*Ms*(-Ad + Mb + Md)*(-Au - Mc + Mt)/(((Mb + Ms)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)));
      
      c = Ms*(-Ad + Mb - Ms)*(Ad - Md + Ms)*(Au + Mc - Mu)*Mu*(-Au + Mt + Mu)/
	(((-Ad + Mb + Md - Ms)*(Mb + Ms)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[1], a, b, c, -1, +1, +1, d1, d2);
    }
    
    if ( proc_matrix_element[2] ) {
      
      /*-- Calculate the various terms in the matrix element M13 --*/
      a = Md*Ms*Mc*Mt*(Ad - Md + Ms)*(-Au - Mc + Mt)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Mb + Ms))*((Mt - Mu)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Mb*Mu*(Ad - Md + Ms)*(-Au - Mc + Mt)/(((Mb - Md)*(Mb + Ms))*((Mt - Mu)*(Mc + Mu)));
      
      c = Mu*Mb*(-Ad + Mb + Md)*(-Ad + Mb - Ms)*(Au + Mc - Mu)*(-Au + Mt + Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Mb + Ms))*((Mt - Mu)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[2], a, b, c, +1, +1, -1, d1, d2);
    }
    
    if ( proc_matrix_element[3] ) {
      
      /*-- Calculate the various terms in the matrix element M21 --*/
      a = Mt*Mu*Mb*Ms*(-Ad + Mb - Ms)*(-Au + Mt + Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Md*Mc*(-Ad + Mb - Ms)*(-Au + Mt + Mu)/(((Mb - Md)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)));
      
      c = Md*Mc*(-Ad + Mb + Md)*(Ad - Md + Ms)*(-Au - Mc + Mt)*(Au + Mc - Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[3], a, b, c, -1, +1, +1, d1, d2);
    }
    
    if ( proc_matrix_element[4] ) {
      
      /*-- Calculate the various terms in the matrix element M22 --*/
      a = Mt*Mu*Mb*Md*(-Ad + Mb + Md)*(-Au + Mt + Mu)/
	(((-Ad + Mb + Md - Ms)*(Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Ms*(-Ad + Mb + Md)*Mc*(-Au + Mt + Mu)/(((Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)));
      
      c = Ms*Mc*(-Ad + Mb - Ms)*(Ad - Md + Ms)*(-Au - Mc + Mt)*(Au + Mc - Mu)/
	(((-Ad + Mb + Md - Ms)*(Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[4], a, b, c, +1, +1, +1, d1, d2);
    }
    
    if ( proc_matrix_element[5] ) {
      
      /*-- Calculate the various terms in the matrix element M23 --*/
      a = Mt*Mu*Md*Ms*(Ad - Md + Ms)*(-Au + Mt + Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Mb + Ms))*((Mc + Mt)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Mc*Mb*(Ad - Md + Ms)*(-Au + Mt + Mu)/(((Mb - Md)*(Mb + Ms))*((Mc + Mt)*(Mc + Mu)));
      
      c = Mc*Mb*(-Ad + Mb + Md)*(-Ad + Mb - Ms)*(-Au - Mc + Mt)*(Au + Mc - Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Mb + Ms))*((Mc + Mt)*(Mc + Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[5], a, b, c, -1, +1, -1, d1, d2);
    }
    
    if ( proc_matrix_element[6] ) {
      
      /*-- Calculate the various terms in the matrix element M31 --*/
      a = Mu*Ms*Mc*Mb*(-Ad + Mb - Ms)*(Au + Mc - Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Md*Mt*(-Ad + Mb - Ms)*(Au + Mc - Mu)/(((Mb - Md)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)));
      
      c = Md*Mt*(-Ad + Mb + Md)*(Ad - Md + Ms)*(-Au - Mc + Mt)*(-Au + Mt + Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[6], a, b, c, +1, +1, -1, d2, d2);
    }
    
    if ( proc_matrix_element[7] ) {
      
      /*-- Calculate the various terms in the matrix element M32 --*/
      a = Mu*Mc*Mb*Md*(-Ad + Mb + Md)*(Au + Mc - Mu)/
	(((-Ad + Mb + Md - Ms)*(Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Ms*Mt*(-Ad + Mb + Md)*(Au + Mc - Mu)/(((Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)));
      
      c = Ms*Mt*(-Ad + Mb - Ms)*(Ad - Md + Ms)*(-Au - Mc + Mt)*(-Au + Mt + Mu)/
	(((-Ad + Mb + Md - Ms)*(Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[7], a, b, c, -1, +1, -1, d1, d2);
    }

    if ( proc_matrix_element[8] ) {
      
      /*-- Calculate the various terms in the matrix element M33 --*/ 
      a = Mu*Md*Ms*Mc*(Ad - Md + Ms)*(Au + Mc - Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Mb + Ms))*((Mc + Mt)*(Mt - Mu)*(-Au - Mc + Mt + Mu)));
      
      b = Mb*Mt*(Ad - Md + Ms)*(Au + Mc - Mu)/(((Mb - Md)*(Mb + Ms))*((Mc + Mt)*(Mt - Mu)));
      
      c = Mb*Mt*(-Ad + Mb + Md)*(-Ad + Mb - Ms)*(-Au - Mc + Mt)*(-Au + Mt + Mu)/
	(((Mb - Md)*(-Ad + Mb + Md - Ms)*(Mb + Ms))*((Mc + Mt)*(Mt - Mu)*(-Au - Mc + Mt + Mu)));
      
      complexroots(&m[8], a, b, c, +1, +1, +1, d1, d2);
    }
    
    return;
  }

  static void calculate_Y_matrix(complex *m, float Yu, float Yd, float d1, float d2) {

    float a,b,c;

    if ( proc_matrix_element[0] ) {

      /*-- Compute the value of M11 --*/
      a = Mb*Ms*Mc*Mt*(Mb - Ms - Yd)*(-Mc + Mt - Yu)/
	(((Mb - Md)*(Md + Ms)*(Mb + Md - Ms - Yd))*
	 ((Mt - Mu)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));

      b = Md*(Mb - Ms - Yd)*Mu*(-Mc + Mt - Yu)/(((Mb - Md)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)));

      c = Md*(Mb + Md - Yd)*(-Md + Ms + Yd)*Mu*(Mt + Mu - Yu)*(Mc - Mu + Yu)/
	(((Mb - Md)*(Md + Ms)*(Mb + Md - Ms - Yd))*
	 ((Mt - Mu)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      
      complexroots( &m[0], a, b, c, +1, +1, +1, d1, d2);
    }
    
    if ( proc_matrix_element[1] ) {
      
      /*-- Compute the value of M12 --*/
      a =  Mb*Md*Mc*Mt*(Mb + Md - Yd)*(-Mc + Mt - Yu)/
	(((Mb + Ms)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mt - Mu)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      
      b = Ms*(Mb + Md - Yd)*Mu*(-Mc + Mt - Yu)/(((Mb + Ms)*(Md + Ms))*((Mt - Mu)*(Mc + Mu)));
      
      c = Ms*(Mb - Ms - Yd)*(-Md + Ms + Yd)*Mu*(Mt + Mu - Yu)*(Mc - Mu + Yu)/
	(((Mb + Ms)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mt - Mu)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      
      complexroots( &m[1], a, b, c, -1, +1, +1, d1, d2);
    }
    
    
    if ( proc_matrix_element[2] ) {
      
      /*-- Compute the value of M13 --*/
      a = Md*Ms*Mc*Mt*(-Md + Ms + Yd)*(-Mc + Mt - Yu)/
	(((Mb - Md)*(Mb + Ms)*(Mb + Md - Ms - Yd))*
	 ((Mt - Mu)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      
      b = Mu*Mb*(-Md + Ms + Yd)*(-Mc + Mt - Yu)/
	(((Mb - Md)*(Mb + Ms))*((Mt - Mu)*(Mc + Mu)));
      
      c = Mu*Mb*(Mb + Md - Yd)*(Mb - Ms - Yd)*(Mt + Mu - Yu)*(Mc - Mu + Yu)/
	(((Mb - Md)*(Mb + Ms)*(Mb + Md - Ms - Yd))*
	 ((Mt - Mu)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));;
      
      complexroots( &m[2], a, b, c, +1, +1, -1, d1, d2);
    }
    
    if ( proc_matrix_element[3] ) {
      
      /*-- Compute the value of M21 --*/
      a = Mt*Mu*Mb*Ms*(Mb - Ms - Yd)*(Mt + Mu - Yu)/
	(((Mb - Md)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      b = Md*Mc*(Mb - Ms - Yd)*(Mt + Mu - Yu)/(((Mb - Md)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)));
      c = Md*(Mb + Md - Yd)*(-Md + Ms + Yd)*Mc*(-Mc + Mt - Yu)*(Mc - Mu + Yu)/
	(((Mb - Md)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      
      complexroots( &m[3], a, b, c, -1, +1, +1, d1, d2 );
    }
    
    if ( proc_matrix_element[4] ) {
      
      /*-- Calculate the value of M22 --*/
      a = Mb*Md*Mt*Mu*(Mb + Md - Yd)*(Mt + Mu - Yu)/
	(((Mb + Ms)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      b = Ms*Mc*(Mb + Md - Yd)*(Mt + Mu - Yu)/(((Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mc + Mu)));
      c = Ms*Mc*(Mb - Ms - Yd)*(-Md + Ms + Yd)*(-Mc + Mt - Yu)*(Mc - Mu + Yu)/
	(((Mb + Ms)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mc + Mu)*(-Mc + Mt + Mu - Yu)));
      
      complexroots(&m[4], a, b, c, +1, +1, +1, d1, d2);
    }
    
    if ( proc_matrix_element[5] ) {
    
      /*-- compute the value of M23 --*/
      a = Md*Ms*(-Md + Ms + Yd)*Mt*Mu*(Mt + Mu - Yu)/
	((Mb - Md)*(Mb + Ms)*(Mb + Md - Ms - Yd)*
	 (Mc + Mt)*(Mc + Mu)*(-Mc + Mt + Mu - Yu));
      
      b = Mc*Mb*(-Md + Ms + Yd)*(Mt + Mu - Yu)/
	((Mb - Md)*(Mb + Ms)*(Mc + Mt)*(Mc + Mu));
      
      c = Mc*Mb*(Mb + Md - Yd)*(Mb - Ms - Yd)*(-Mc + Mt - Yu)*(Mc - Mu + Yu)/
	((Mb - Md)*(Mb + Ms)*(Mb + Md - Ms - Yd)*
	 (Mc + Mt)*(Mc + Mu)*(-Mc + Mt + Mu - Yu));
      
      complexroots( &m[5], a, b, c, -1, +1, -1, d1, d2);
    }

    
    if ( proc_matrix_element[6] ) {

      /*-- Compute the value of M31 --*/
      a = Mb*Ms*Mc*Mu*(Mb - Ms - Yd)*(Mc - Mu + Yu)/
	(((Mb - Md)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mt - Mu)*(-Mc + Mt + Mu - Yu)));
      b = Md*Mt*(Mb - Ms - Yd)*(Mc - Mu + Yu)/(((Mb - Md)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)));
      c = Md*Mt*(Mb + Md - Yd)*(-Md + Ms + Yd)*(-Mc + Mt - Yu)*(Mt + Mu - Yu)/
	(((Mb - Md)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mt - Mu)*(-Mc + Mt + Mu - Yu)));

      complexroots(&m[6], a, b, c, +1, +1, -1, d1, d2);
    }
      

    if ( proc_matrix_element[7] ) {
      
      /*-- Compute the value of M32 --*/
      a = Mc*Mu*Mb*Md*(Mb + Md - Yd)*(Mc - Mu + Yu)/
	(((Mb + Ms)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mt - Mu)*(-Mc + Mt + Mu - Yu)));
      b = Mt*Ms*(Mb + Md - Yd)*(Mc - Mu + Yu)/(((Mb + Ms)*(Md + Ms))*((Mc + Mt)*(Mt - Mu)));
      c = Ms*Mt*(Mb - Ms - Yd)*(-Md + Ms + Yd)*(-Mc + Mt - Yu)*(Mt + Mu - Yu)/
	(((Mb + Ms)*(Md + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mt - Mu)*(-Mc + Mt + Mu - Yu)));
      
      complexroots(&m[7], a, b, c, -1, +1, -1, d1, d2);
    }
    
    if ( proc_matrix_element[8] ) {
      
      /*-- Calculate the value of M33 --*/
      a = Md*Ms*Mc*Mu*(-Md + Ms + Yd)*(Mc - Mu + Yu)/
	(((Mb - Md)*(Mb + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mt - Mu)*(-Mc + Mt + Mu - Yu)));
      b = Mb*Mt*(-Md + Ms + Yd)*(Mc - Mu + Yu)/(((Mb - Md)*(Mb + Ms))*((Mc + Mt)*(Mt - Mu)));
      c = Mb*Mt*(Mb + Md - Yd)*(Mb - Ms - Yd)*(-Mc + Mt - Yu)*(Mt + Mu - Yu)/
	(((Mb - Md)*(Mb + Ms)*(Mb + Md - Ms - Yd))*((Mc + Mt)*(Mt - Mu)*(-Mc + Mt + Mu - Yu)));
      
      complexroots(&m[8], a, b, c, +1, +1, +1, d1, d2);
    }
    
    return;
  }
