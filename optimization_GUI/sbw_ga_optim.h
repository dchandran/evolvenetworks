#ifndef SBW_GALIB_OPTIMIZATION_H
#define SBW_GALIB_OPTIMIZATION_H

#include <vector>
#include <string>
#include "SBW/sbwenums.h"
#include "SBW/SBWException.h"
#include "SBW/DataBlockType.h"
#include "SBW/SBW.h"
#include "muParserDef.h"
#include "muParser.h"
#include "muParserInt.h"
#include "sbml_sim.h"

class SBW_GA_optimizer
{

public:

	SBW_GA_optimizer();

	~SBW_GA_optimizer();

	//API related to the optimization algorithms
	SystemsBiologyWorkbench::DataBlockWriter getTypesOfObjectiveFunctions( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );

	SystemsBiologyWorkbench::DataBlockWriter setTypeOfObjectiveFunction( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter getAlgorithmNames(	SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter setAlgorithmByName( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter setAlgorithmByIndex( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter getAlgorithmParameterNames( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter getAlgorithmParameterValues( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter setAlgorithmParameterByName( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter setAlgorithmParameterByIndex( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	//inputs

	SystemsBiologyWorkbench::DataBlockWriter loadData( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );

	SystemsBiologyWorkbench::DataBlockWriter setTargetFlux( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );

	//API related to the model

	SystemsBiologyWorkbench::DataBlockWriter loadSBML( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter getSBML( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter getModelParameterNames( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter getModelParameterValues( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter setModelParameterValueByName( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );
	
	SystemsBiologyWorkbench::DataBlockWriter setModelParameterValueByIndex( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );

	//run

	SystemsBiologyWorkbench::DataBlockWriter optimize( SystemsBiologyWorkbench::Module, SystemsBiologyWorkbench::DataBlockReader );

	static void registerMethods(MethodTable<SBW_GA_optimizer> & table);

private:
	
	//simulator and optimizer
	SBML_sim * simulator;
	
	enum AlgorithmTypes { GA = 0, Crowding };
	AlgorithmTypes algorithmType;
	
	enum ObjectiveFunctionTypes { TimeCourse = 0, SteadyState, MaximizeFlux };
	ObjectiveFunctionTypes objectiveFunction;

	void readTargetData(const std::string& filename, const std::string& delim=",");
	std::vector< std::vector<double> > targetData;
};

#endif
