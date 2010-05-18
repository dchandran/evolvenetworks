#include "sbw_ga_optim.h"
#include "muParserDef.h"
#include "muParser.h"
#include "muParserInt.h"
#include "sbml/SBMLReader.h"
#include "sbml/SBMLDocument.h"
#include "sbml/Model.h"
#include "GASimpleGA.h"
#include "GARealGenome.h"
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

extern "C"
{
	#include "cvodesim.h"
}

using namespace std;
using namespace SystemsBiologyWorkbench;

SBW_GA_optimizer::SBW_GA_optimizer() : algorithmType(GA), objectiveFunction(TimeCourse), stoichiometryMatrix(0)
{
	
}

DataBlockWriter SBW_GA_optimizer::getTypesOfObjectiveFunctions( Module , DataBlockReader )
{
	vector<string> names;
	names.push_back("TimeCourse");
	names.push_back("SteadyState");
	names.push_back("MaximizeFlux");
	return DataBlockWriter() << names;

}

DataBlockWriter SBW_GA_optimizer::setTypeOfObjectiveFunction( Module , DataBlockReader reader)
{
	string name;	
	reader >> name;
	
	if (name.compare("TimeCourse") == 0)
		 objectiveFunction = TimeCourse;
	if (name.compare("SteadyState") == 0)
		objectiveFunction = SteadyState;
	if (name.compare("MaximizeFlux") == 0)
		objectiveFunction = MaximizeFlux;
		
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::getAlgorithmNames( Module, DataBlockReader )
{
	vector<string> names;
	names.push_back("GA");
	names.push_back("GA with crowding");
	return DataBlockWriter() << names;
}

DataBlockWriter SBW_GA_optimizer::setAlgorithmByName( Module, DataBlockReader reader)
{
	string name;
	reader >> name;
	
	if (name.compare("GA") == 0)
		algorithmType = GA;
	else
		algorithmType = Crowding;
		
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::setAlgorithmByIndex( Module, DataBlockReader reader)
{
	int i;
	reader >> i;
	
	if (i == 0)
		algorithmType = GA;
	else
		algorithmType = Crowding;
		
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::getAlgorithmParameterNames( Module, DataBlockReader )
{
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::getAlgorithmParameterValues( Module, DataBlockReader )
{
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::setAlgorithmParameterByName( Module, DataBlockReader )
{
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::setAlgorithmParameterByIndex( Module, DataBlockReader )
{
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::loadSBML( Module, DataBlockReader reader)
{
	string filename = "oscil.sbml";
	//reader >> filename;
	
	SBMLReader * sbmlreader = new SBMLReader;
	SBMLDocument * doc = sbmlreader->readSBML(filename);

	if (!doc || doc->getNumErrors() > 0)
	{
		printf("error reading SBML file\n");
	}
	else
	{
		Model * model = doc->getModel();
		ListOfParameters * params = model->getListOfParameters();
		ListOfReactions * reacs = model->getListOfReactions();
		ListOfSpecies * species = model->getListOfSpecies();
		ListOfSpeciesTypes * types = model->getListOfSpeciesTypes();
		
		

		parameterNames.clear();
		parameterValues.clear();
		rateEqns.clear();
		variableNames.clear();
		reactionNames.clear();
		if (stoichiometryMatrix)
			delete stoichiometryMatrix;
		
		for (int i=0; i < species->size(); ++i)
			if (!species->get(i)->getConstant() && !species->get(i)->getBoundaryCondition())
			{
				variableNames.push_back(species->get(i)->getId());
			}

		for (int i=0; i < variableNames.size(); ++i)
			cout << "var: " << variableNames[i] << endl;

		for (int i=0; i < params->size(); ++i)
		{
			parameterNames.push_back(params->get(i)->getId());
			parameterValues.push_back(params->get(i)->getValue());
		}

		cout << endl;

		for (int i=0; i < parameterNames.size(); ++i)
			cout << "param: " << parameterNames[i] << " = " << parameterValues[i] << endl;

		int numReacs = reacs->size();

		stoichiometryMatrix = new double[ numReacs * variableNames.size() ];
		
		for (int i=0; i < numReacs; ++i)
		{
			Reaction * r = reacs->get(i);
			reactionNames.push_back(r->getId());
			rateEqns.push_back(r->getKineticLaw()->getFormula());
			ListOfSpeciesReferences * reactants = r->getListOfReactants(),
									* products  = r->getListOfProducts();

			for (int j=0; j < variableNames.size(); ++j)
			{
				stoichiometryMatrix[ j*numReacs + i ] = 0.0;

				for (int k=0; k < reactants->size(); ++k)
					if (reactants->get(k)->getSpecies() == variableNames[j])
						stoichiometryMatrix[ j*numReacs + i ] -= ((SpeciesReference*)(reactants->get(k)))->getStoichiometry();
					
				for (int k=0; k < products->size(); ++k)
					if (products->get(k)->getSpecies() == variableNames[j])
						stoichiometryMatrix[ j*numReacs + i ] += ((SpeciesReference*)(reactants->get(k)))->getStoichiometry();;
			}
		}

		cout << endl;

		for (int i=0; i < rateEqns.size(); ++i)
		{
			cout << "reac: ";
			for (int j=0; j < variableNames.size(); ++j)
				cout << stoichiometryMatrix[ j*numReacs + i ] << variableNames[j] << " ";
			cout << ";" << rateEqns[i] << endl;
		}
		
		//delete params;
		//delete reacs;
		delete doc;
	}
	return DataBlockWriter();
}

SBW_GA_optimizer::~SBW_GA_optimizer()
{
}

DataBlockWriter SBW_GA_optimizer::getSBML( Module, DataBlockReader )
{
	string s;
	return DataBlockWriter() << s;
}

void SBW_GA_optimizer::readTargetData(const string& filename, const string& delim)
{	
	vector<string> strs;
	string line;

	targetDataHeaders.clear();
	
	ifstream fin(filename.c_str());
	getline(fin,line);
	boost::split(targetDataHeaders, line, boost::is_any_of(delim));

	int numCols = targetDataHeaders.size();
	
	vector<double> row,col;
	vector< vector<double> > dat;
	
	while (!fin.eof())
	{
		getline(fin,line);
		boost::split(strs, line, boost::is_any_of(delim));
		
		if (strs.size() == numCols)
		{
			row.clear();
			for (int i=0; i < strs.size(); ++i)
				row.push_back(atof(strs[i].c_str()));
			dat.push_back(row);
		}
	}

	targetData.clear();
	for (int i=0; i < numCols; ++i)
	{
		col.clear();
		for (int j=0; j < dat.size(); ++j)
			col.push_back(dat.at(j).at(i));
		targetData.push_back(col);
	}

	fin.close();

	ofstream fout("temp.dat");
	for (int i=0; i < targetData.size(); ++i)
	{
		cout << targetDataHeaders[i] << "\t";
		col = targetData[i];
		for (int j=0; j < col.size(); ++j)
			fout << col[j] << "\t";
		fout << "\n";
	}
	fout.close();
}

DataBlockWriter SBW_GA_optimizer::getModelParameterNames( Module , DataBlockReader )
{
	return DataBlockWriter() << parameterNames;
}
	
DataBlockWriter SBW_GA_optimizer::getModelParameterValues( Module , DataBlockReader )
{
	return DataBlockWriter() << parameterValues;
}

DataBlockWriter SBW_GA_optimizer::setModelParameterValueByName( Module, DataBlockReader reader)
{
	string s;
	reader >> s;

	double d = 0.0;
	for (int i=0; i < parameterNames.size() && i < parameterValues.size(); ++i)
	{
		if (parameterNames[i] == s)
			d = parameterValues[i];
	}
	return DataBlockWriter() << d;
}

DataBlockWriter SBW_GA_optimizer::setModelParameterValueByIndex( Module, DataBlockReader reader)
{
	int i;
	reader >> i;

	double d = 0.0;
	if (i >= 0 && i < parameterNames.size() && i < parameterValues.size())
	{
		d = parameterValues[i];
	}
	return DataBlockWriter() << d;
}

DataBlockWriter SBW_GA_optimizer::loadData( Module , DataBlockReader reader)
{
	//string filename;
	//reader >> filename;
	readTargetData("target.csv","\t");
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::setTargetFlux( Module , DataBlockReader reader)
{
	string s;
	reader >> s;
	return DataBlockWriter();
}

void * SBW_GA_optimizer::makeModelData()
{
	model_data * u = new model_data;

	vector<double> col;
	vector<vector<double> > targetData2 = targetData;
	targetData.clear();

	targetData.push_back(targetData2.at(0));
	for (int i=0; i < parameterNames.size(); ++i)
	{
		col.clear();
		for (int j=0; j < targetDataHeaders.size(); ++j)
			if (parameterNames[i] == targetDataHeaders[j])
			{
				col = targetData2.at(j);
				break;
			}
		targetData.push_back(col);
	}

	u->numParams = parameterNames.size();
	u->numReactions = rateEqns.size();
	u->numSpecies = variableNames.size();
	u->startTime = targetData[0][0];
	u->endTime = targetData[0][ targetData[0].size() - 1 ];
	u->stepSize = targetData[0][1] - targetData[0][0];
	u->targetData = &targetData;
	u->N = stoichiometryMatrix;
	vector<mu::Parser*> rateEqns;

	return (void*)u;
}

DataBlockWriter SBW_GA_optimizer::optimize( Module , DataBlockReader )
{
	float a[] = {1};
	GARealAlleleSet allele(1,a);
	GARealGenome genome( parameterNames.size(), allele, &Objective1);
	GASimpleGA ga(genome);

	model_data * u = (model_data*)makeModelData();

	ga.userData(u);

	GASharing dist(EuclideanDistance);
	ga.scaling(dist);
	ga.minimize();
	ga.populationSize(1000);
	ga.nGenerations(100);
	ga.pMutation(0.01);
	ga.pCrossover(0.9);
	ga.evolve();
	return DataBlockWriter();
}

void SBW_GA_optimizer::registerMethods(MethodTable<SBW_GA_optimizer> & table)
{
	try
	{
		table.addMethod(
			&SBW_GA_optimizer::getTypesOfObjectiveFunctions, 
			"string[] getTypesOfObjectiveFunctions()",
			false,
			"get a list of the different types of objective functions"
			);
		table.addMethod(
			&setTypeOfObjectiveFunction, 
			"void setTypeOfObjectiveFunction(string)",
			false,
			"set the type of objective function to optimize"
			);
		table.addMethod(
			&getAlgorithmNames, 
			"string[] getAlgorithmNames()",
			false,
			"get all the available algorithm names"
			);
		table.addMethod(
			&setAlgorithmByName, 
			"void setAlgorithmByName(string)",
			false,
			"set the algorithm from one of the strings returned by getAlgorithmNames"
			);
		table.addMethod(
			&setAlgorithmByIndex, 
			"void setAlgorithmByIndex(int)",
			false,
			"set the algorithm an index corresponding to values returned by getAlgorithmNames"
			);
		table.addMethod(
			&getAlgorithmParameterNames, 
			"string[] getAlgorithmParameterNames()",
			false,
			"get all the parameter names used for controlling the algorithm (not the model)"
			);
		table.addMethod(
			&getAlgorithmParameterValues, 
			"double[] getAlgorithmParameterValues()",
			false,
			"get all the parameter values corresponding to getAlgorithmParameterNames"
			);
		table.addMethod(
			&setAlgorithmParameterByName, 
			"void setAlgorithmParameterByName(string,double)",
			false,
			"set one of the algorithm parameter values"
			);
		table.addMethod(
			&setAlgorithmParameterByIndex, 
			"void setAlgorithmParameterByIndex(int, double)",
			false,
			"set one of the algorithm parameter values"
			);
		table.addMethod(
			&loadSBML, 
			"void loadSBML(string)",
			false,
			"load the model in SBML file"
			);
		table.addMethod(
			&getSBML, 
			"string getSBML()",
			false,
			"get the current model in SBML format"
			);
		table.addMethod(
			&getModelParameterNames, 
			"string[] getModelParameterNames()",
			false,
			"get all the parameter names in the model"
			);
		table.addMethod(
			&getModelParameterValues, 
			"double[] getModelParameterValues()",
			false,
			"get all the parameter values in the model"
			);
		table.addMethod(
			&setModelParameterValueByName, 
			"void setModelParameterValueByName(string,double)",
			false,
			"set one of the parameter values in the model"
			);
		table.addMethod(
			&setModelParameterValueByIndex, 
			"void setModelParameterValueByIndex(int,double)",
			false,
			"set one of the parameter values in the model"
			);
		table.addMethod(
			&loadData, 
			"void loadData(string)",
			false,
			"load target data"
			);
		table.addMethod(
			&setTargetFlux, 
			"void setTargetFlux(string)",
			false,
			"set that flux to optimize. Use setTypeOfObjectiveFunction"
			);
		table.addMethod(
			&optimize, 
			"void optimize()",
			false,
			"optimize the function set bu using setTypeOfObjectiveFunction"
			);
	}
	catch(...)
	{
		printf("error: could not add methods\n");
	}
}

int main(int argc, char* argv[])
{
	/*try
	{
		ModuleImpl modImpl(
			"GAlib",
			"Genetic algorithm written in C++",
			UniqueModule);
		
		modImpl.addServiceObject(
			"optim",
			"Optimization using genetic algorithm",
			"optimization",
			new SBW_GA_optimizer());
		
		modImpl.run(argc,argv);
	}
	catch(SBWException * e)
	{
		printf("error: %s\n",e->getMessage());
	}*/

	SBW_GA_optimizer gaOpt;
	gaOpt.loadSBML(Module(), DataBlockReader());
	gaOpt.loadData(Module(), DataBlockReader());
	
	return 0;
}
