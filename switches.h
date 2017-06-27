// The top most header - contains switches and constants
// Inherited in the rod code as well
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "mpi.h"

#define resume

#define MPIreal MPI_DOUBLE
#define real double
#define Pi   acos(-1.0)

#define pinning
//#define pinnedMTfrac 0.5 // for circular cell:- 1-> all pinned, 0-> none pinned
#define frictionattip 1

#define centrofactor 25
#define Nmax 256 // Max segments allowed

// Cell shape and size

#define square
#define thickness 1  // Cell thickness in microns

//#define circle
#define Dia       40 // Cell diameter or square side in microns

//#define motorforceON
#define excludedvolON
//#define stallforceON

// Excluded volume parameters - in rod units
#define U0         500 //units of 1e-21J
#define lambda     50  //nanometer
#define cutoffdist 100  //nanometer


// clamped or hinged MTs BC


// variables shared with the rod code
int iprint, iter;
char datadir[64];
