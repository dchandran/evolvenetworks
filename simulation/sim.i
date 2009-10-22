%module sim
%{
#include "cvodesim.h"
#include "ssa.h"
#include "langevin.h"
#include "cells_ssa.h"
%}

typedef void (*ODEFunction)(double time,double * y,double* dydt,void* params);
typedef void (*PropensityFunction)(double time,double * y,double* rates,void* params);
typedef double (*EventFunction)(double time,double * y,void* params);
extern void ODEevents(int numEvents, EventFunction * eventFunctions);
extern void ODEflags(int);
extern void ODEtolerance(double,double);
extern double* ODEsim(int N, double * iv, ODEFunction function, double startTime, double endTime, double stepSize, void * params);
extern double * ODEsim2(int, int, double * iv, PropensityFunction f, double * point, double, double, double,void*);
extern double* jacobian(int N, double * iv,  ODEFunction function, void * params);
extern double* jacobian2(int m, int n, double * N, PropensityFunction f, double * point, void * params);
extern double* steadyState(int N, double * iv, ODEFunction function, void * params, double minerr, double maxtime, double delta);
extern double* steadyState2(int m, int n, double * N, PropensityFunction f, double * iv, void * params, double minerr, double maxtime, double delta);
extern double* getDerivatives(int N, double * iv, ODEFunction function, double startTime, double endTime, double stepSize, void * params);
extern double* getDerivatives2(int m, int n, double * N, PropensityFunction f, double * iv, double startTime, double endTime, double stepSize, void * params);
extern double * SSA(int, int, double * iv, PropensityFunction, double*, double, double, int, int*,void*);
extern double * getRatesFromSimulatedData(double * data, int rows, int cols1, int cols2, int skip, PropensityFunction, void* param);
extern double * Langevin(int n, int m, double * N, PropensityFunction propensity, double * iv, double endTime, double dt, void * params);
extern double ** cells_ssa(int, int, double * N, PropensityFunction, double * iv, double, int, void*, int, double, double, double);

%inline %{
/* Create a new array */
double * new_array(int n) {
	return (double *) malloc(n * sizeof(double));
}
%}
%inline %{
/* get a value from a 1D array */
double get_value(double * array, int i) {
	return array[i];
}
%}
%inline %{
/* get a value from a 2D array */
double get_value2D(double * array, int N, int i, int j) {
	return array[ (((i)*(N)) + (j)) ];
}
%}
%inline %{
/* set a value in a 1D array */
void set_value(double * array, int i, double value) {
	array[i] = value;
}
%}
%inline %{
/* set a value in a 2D array */
void set_value2D(double * array, int N, int i, int j, double value) {
	array[ (((i)*(N)) + (j)) ] = value;
}
%}

typedef void (*ODEFunction)(double time,double * y,double* dydt,void* params);
typedef void (*PropensityFunction)(double time,double * y,double* rates,void* params);
typedef double (*EventFunction)(double time,double * y,void* params);

%inline %{


/* convert a python ODE callback function to C*/
static PyObject * SIM_PYTHON_ODE_FUNCTION_CALLBACK = NULL;

void SIM_PYTHON_ODE_FUNCTION(double time,double * y,double* dydt,void* params)
{
	PyObject *arglist;
	if (SIM_PYTHON_ODE_FUNCTION_CALLBACK)
	{		
		/* call the callback */
		arglist = Py_BuildValue("(dO&O&O&)", time, y, dy, params);
		PyEval_CallObject(SIM_PYTHON_ODE_FUNCTION_CALLBACK, arglist);
		Py_DECREF(arglist);
	}
}

static PyObject * sim_create_callback(PyObject *dummy, PyObject *args)
{
    PyObject *temp;

    if (PyArg_ParseTuple(args, "O:set_callback", &temp)) {
        if (!PyCallable_Check(temp)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return NULL;
        }
        Py_XINCREF(temp);         /* Add a reference to new callback */
		if (SIM_PYTHON_ODE_FUNCTION_CALLBACK)
			Py_XDECREF(SIM_PYTHON_ODE_FUNCTION_CALLBACK);  /* Dispose of previous callback */
        SIM_PYTHON_ODE_FUNCTION_CALLBACK = temp;       /* Remember new callback */
    }
    return Py_BuildValue("i",(int)(&SIM_PYTHON_ODE_FUNCTION));
}

static void sim_delete_callback()
{
	if (SIM_PYTHON_ODE_FUNCTION_CALLBACK)
	{
		Py_XDECREF(SIM_PYTHON_ODE_FUNCTION_CALLBACK);
	}
}

%}

