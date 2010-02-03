#ifndef PARAMETER_DISTRIBUTION_MONTECARLO_H
#define PARAMETER_DISTRIBUTION_MONTECARLO_H

/*! \brief generate random parameters
* \param double* return value where random parameters are assigned
*/
typedef void (*RandomParameterFunc)(double*);
/*! \brief return GOOD (1) or BAD (0) for the given set of parameters
* \param double* parameters
* \return int 1 or 0 for accept or reject this array of parameter
*/
typedef int (*TargetFunc)(double*);

/*! \brief generate random parameters and run the given function
* \param int number of samples to generate
* \param int number of parameters to randomize
* \param RandomParameterFunc the function for generating random parameter arrays
* \param TargetFunc the function for separating good and bad parameters
* \param double* (pre-allocated return value) mean vector for the distribution of acceptable parameters
* \param double** (pre-allocated return value) covariance matrix for the distribution of acceptable parameters
*/
void sampleParameterSpace(int numSamples, int numParameters, RandomParameterFunc, TargetFunc, double * , double **);

/*! \brief set the number of threads to generate when running the sampling algorithm
* \param int number of threads
*/
void setNumThreads(int);

#endif
