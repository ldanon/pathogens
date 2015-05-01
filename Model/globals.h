#ifndef _GLOBALS_H
#define _GLOBALS_H

//#define VTKNET
#define _CRT_SECURE_NO_DEPRECATE 

#define N_NODE 10000
#define NUM_COM 2
#define NUM_COM_SIZE 256
#define ITERATIONS 10000
#define NUM_LOCI 2
#define NUM_STRAINS 4
#define REALS 1

#define max_size 30000

#define TRACK 35 // node to track

extern long seed;


#include <gsl/gsl_statistics_double.h>   // gnu scientific library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

#endif
