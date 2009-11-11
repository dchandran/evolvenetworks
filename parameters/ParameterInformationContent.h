#ifndef PARAMETER_INFORMATION_CONTENT_MONTE_CARLO
#define PARAMETER_INFORMATION_CONTENT_MONTE_CARLO
#include "mtrand.h"
#include "cvodesim.h"

typedef double (*OutputFunction)(double*);

double * parameterMI(OutputFunction f, int numParams,double * low, double * hi);

#endif
