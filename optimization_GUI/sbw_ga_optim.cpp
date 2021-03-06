#include "sbw_ga_optim.h"
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

extern "C"
{
	#include "cvodesim.h"
}

using namespace std;
using namespace SystemsBiologyWorkbench;

SBW_GA_optimizer::SBW_GA_optimizer() : algorithmType(GA), objectiveFunction(TimeCourse), stoichiometryMatrix(0), simulator(0)
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
	string filename = "feedback.xml";
	reader >> filename;
	
	simulator = new SBML_sim(filename);
	
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
	readTargetData("target.csv","\t");
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::setTargetFlux( Module , DataBlockReader reader)
{
	string s;
	reader >> s;
	return DataBlockWriter();
}

DataBlockWriter SBW_GA_optimizer::optimize( Module , DataBlockReader )
{
	if (simulator)
		simulator->optimize(
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
