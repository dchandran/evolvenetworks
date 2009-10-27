%module ga
%{
#include "ga.h"
#include "reactionNetwork.h"
%}
typedef void* GAindividual;
typedef GAindividual* GApopulation;
typedef void (*GADeleteFnc)(GAindividual);
typedef GAindividual (*GACloneFnc)(GAindividual);
typedef double(*GAFitnessFnc)(GAindividual);
typedef GAindividual (*GACrossoverFnc)(GAindividual, GAindividual);
typedef GAindividual (*GAMutateFnc)(GAindividual);
typedef int(*GASelectionFnc)(GApopulation , double * , double , int );
typedef int(*GACallbackFnc)(int iter,GApopulation,int popSz);
extern void GAinit(GADeleteFnc, GACloneFnc ,GAFitnessFnc, GACrossoverFnc, GAMutateFnc, GASelectionFnc);
extern GApopulation GArun(GApopulation,int sz0,int sz1,int maxIter, GACallbackFnc);
extern int GArouletteWheelSelection(GApopulation , double * , double , int );
extern int GAtournamentSelection(GApopulation , double * , double , int );
extern int GAeliteSelection(GApopulation , double * , double , int );
extern void GAsetupNewStruct(GADeleteFnc, GACloneFnc);
extern void GAsetFitnessFunction(GAFitnessFnc);
extern void GAsetCrossoverFunction(GACrossoverFnc);
extern void GAsetMutationFunction(GAMutateFnc);
extern void GAsetSelectionFunction(GASelectionFnc);
extern GAFitnessFnc GAgetFitnessFunction();
extern GACrossoverFnc GAgetCrossoverFunction();
extern GAMutateFnc GAgetMutationFunction();
extern GASelectionFnc GAgetSelectionFunction();
extern GApopulation GAnextGen(int,GApopulation,int,int,short);
extern void GAsort(GApopulation, GAFitnessFnc, int);
extern void GAfree(GApopulation population);
typedef struct 
{
	int type;
	GAindividual network;
	double * initialValues;
}
ReactionNetwork;
extern void printNetwork(FILE * stream, GAindividual);
extern void printNetworkToFile( char * filename, GAindividual);
extern int getNumSpecies(GAindividual);
extern int getNumReactions(GAindividual);
extern double* getStoichiometryMatrix(GAindividual);
extern double* getReactionRates(GAindividual, double*);
extern void lineageTrackingON();
extern void lineageTrackingOFF();
extern int* getOriginalParents(int,int);
extern int* getImmediateParents(int,int);
extern void setRatesFunction( int, PropensityFunction );
extern void setStoichiometryFunction( int, double* (*f)(GAindividual) );
extern void setInitialValues( GAindividual, double *);
extern double * getInitialValues( GAindividual );
extern double * simulateNetworkODE( GAindividual, double, double );
extern double * networkSteadyState( GAindividual );
extern double * simulateNetworkStochastically( GAindividual, double, int* );
extern void setFitnessFunction(GAFitnessFnc);
extern void setCrossoverFunction( int , GACrossoverFnc);
extern void setMutationFunction( int , GAMutateFnc );
extern void setNetworkTypeProbability(int, double);
extern void setNetworkType(int);
extern void setNetworkSize(int min_vars,int max_vars,int min_reactions,int max_reactions);
extern GApopulation evolveNetworks(int init_popSz,int final_popSz,int iterations,GACallbackFnc callback);
extern GApopulation randomNetworks(int);
extern void deleteNetwork(GAindividual);
extern GAindividual cloneNetwork(GAindividual);
extern void setCrossoverRate(double);
extern void setAverageInitialValue(double);
extern void setMutationRateOfInitialValues(double);
extern GAindividual mutateNetwork(GAindividual);
extern GAindividual crossoverNetwork(GAindividual, GAindividual);
extern double compareSteadyStates(GAindividual, double **, int , int, int, int, double ** );
extern void enableLogFile(char * filename);
extern void disableLogFile();
extern void configureContinuousLog(int bestNetworkFitness, 
							int bestNetworkScript,
							int bestNetworkSize, 
							int bestNetworkLineage,
							int allFitness,
							int allNetworkLineage );
extern void configureFinalLog(int bestNetworkFitness, 
							int bestNetworkScript,
							int bestNetworkSize, 
							int bestNetworkLineage,
							int allFitness, 
							int allNetworkLineage,
							int seeds);
extern void configureSteadyStateFunction(double tolerance, 
									double delta,
									double maxTime);
