%module sim
%{
#include "cvodesim.h"
#include "ssa.h"
#include "langevin.h"
#include "cells_ssa.h"
%}

typedef void (*ODEFunction)(double time,double y[],double* dydt,void* params);
typedef void (*PropensityFunction)(double time,double y[],double* rates,void* params);
typedef double (*EventFunction)(double time,double y[],void* params);
extern void ODEevents(int numEvents, EventFunction * eventFunctions);
extern void ODEflags(int);
extern void ODEtolerance(double,double);
extern double* ODEsim(int N, double iv[], ODEFunction function, double startTime, double endTime, double stepSize, void * params);
extern double * ODEsim2(int, int, double iv[], PropensityFunction f, double point[], double, double, double,void*);
extern double* jacobian(int N, double iv[],  ODEFunction function, void * params);
extern double* jacobian2(int m, int n, double N[], PropensityFunction f, double point[], void * params);
extern double* steadyState(int N, double iv[], ODEFunction function, void * params, double minerr, double maxtime, double delta);
extern double* steadyState2(int m, int n, double N[], PropensityFunction f, double iv[], void * params, double minerr, double maxtime, double delta);
extern double* getDerivatives(int N, double iv[], ODEFunction function, double startTime, double endTime, double stepSize, void * params);
extern double* getDerivatives2(int m, int n, double N[], PropensityFunction f, double iv[], double startTime, double endTime, double stepSize, void * params);
extern double * SSA(int, int, double iv[], PropensityFunction, double*, double, double, int, int*,void*);
extern double * getRatesFromSimulatedData(double data[], int rows, int cols1, int cols2, int skip, PropensityFunction, void* param);
extern double * Langevin(int n, int m, double N[], PropensityFunction propensity, double iv[], double endTime, double dt, void * params);
extern double ** cells_ssa(int, int, double N[], PropensityFunction, double iv[], double, int, void*, int, double, double, double);
