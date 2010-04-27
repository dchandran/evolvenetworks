#include <iostream>
#include <vector>
#include "GASStateGA.h"
#include "GA1DArrayGenome.h"
#include "sbml_sim.h"

using namespace std;

void sbml_rates_function(double t, double * y, double * rates, void * data)
{
	int i;
	SBML_sim * u = (SBML_sim*)data;
	
	for (i=0; i < u->variableValues.size(); ++i)
	{
		u->variableValues[i] = y[i];
	}

	for (i=0; i < u->assignmentEqns.size(); ++i)
	{
		u->assignmentValues[i] = u->assignmentEqns[i].Eval();
	}

	for (i=0; i < u->rateEqns.size(); ++i)
	{
		rates[i] = u->rateEqns[i].Eval();
	}
}

int sbml_event_function(int i, double t, double * y, void * data)
{
	SBML_sim * u = (SBML_sim*)data;
	return u->triggerEqns[i].Eval();
}

void sbml_response_function(int i, double * y, void * data)
{
	SBML_sim * u = (SBML_sim*)data;
	for (int j=0; j < u->responseEqns[i].size(); ++j)
		 u->responseEqns[i][j].Eval();
}

static double power(double x, double e) {	return pow(x,e); }
static double ge(double x, double y) { 	return (double)(x >= y); }
static double gt(double x, double y) { 	return (double)(x > y); }
static double le(double x, double y) { 	return (double)(x <= y); }
static double lt(double x, double y) { 	return (double)(x < y); }


static void addSBMLFunctions(mu::Parser & p)
{
	p.DefineFun("pow", &power , false);
	p.DefineFun("ge", &ge , false);
	p.DefineFun("gt", &gt , false);
	p.DefineFun("le", &le , false);
	p.DefineFun("lt", &lt , false);
}

SBML_sim::SBML_sim(string sbml_text, bool isFile)
{
	SBMLReader * sbmlreader = new SBMLReader;
	SBMLDocument * doc;
	
	if (isFile)
		doc = sbmlreader->readSBML(sbml_text);
	else
		doc = sbmlreader->readSBMLFromString(sbml_text); 

	if (!doc || doc->getNumErrors() > 0)
	{
	}
	else
	{
		Model * model = doc->getModel();
		ListOfParameters * params = model->getListOfParameters();
		ListOfReactions * reacs = model->getListOfReactions();
		ListOfSpecies * species = model->getListOfSpecies();
		ListOfSpeciesTypes * types = model->getListOfSpeciesTypes();
		ListOfEvents * events = model->getListOfEvents();
		ListOfRules * rules = model->getListOfRules();
		
		vector<string> assignmentEquations, rateEquations, eventTriggers;
		vector< vector<string> > eventResponses;

		for (int i=0; i < events->size(); ++i)
		{
			Event * e = events->get(i);
			eventTriggers.push_back( SBML_formulaToString( e->getTrigger()->getMath() ) );
			ListOfEventAssignments * eventAssn = e->getListOfEventAssignments();
			vector<string> responses;
			string s;
			for (int j=0; j < eventAssn->size(); ++j)
			{
				s = eventAssn->get(j)->getVariable();
				s.append("=");
				s.append( SBML_formulaToString( eventAssn->get(j)->getMath() ) );
				responses.push_back(s);
			}

			eventResponses.push_back( responses );
		}
		
		for (int i=0; i < rules->size(); ++i)
		{
			Rule * r = rules->get(i);
			
			if (r->isAssignment())
			{
				AssignmentRule * ar  = (AssignmentRule*)r;
				assignmentVariables.push_back(ar->getVariable());
				assignmentValues.push_back(0.0);
				assignmentEquations.push_back(ar->getFormula());
			}
		}
		
		for (int i=0; i < species->size(); ++i)
			if (!species->get(i)->getConstant() && !species->get(i)->getBoundaryCondition())
			{
				variableNames.push_back(species->get(i)->getId());
				variableValues.push_back(species->get(i)->getInitialConcentration());
			}
			else
			{
				parameterNames.push_back(species->get(i)->getId());
				parameterValues.push_back(species->get(i)->getInitialConcentration());
			}
			
		initialValues = variableValues;
		originalParameters = parameterValues;

		for (int i=0; i < params->size(); ++i)
		{
			parameterNames.push_back(params->get(i)->getId());
			parameterValues.push_back(params->get(i)->getValue());
		}

		int numReacs = reacs->size();

		stoichiometryMatrix = new double[ numReacs * variableNames.size() ];
		
		for (int i=0; i < numReacs; ++i)
		{
			Reaction * r = reacs->get(i);
			reactionNames.push_back(r->getId());
			rateEquations.push_back(r->getKineticLaw()->getFormula());
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
		
		for (int i=0; i < rateEquations.size(); ++i)
		{
			mu::Parser p;
			addSBMLFunctions(p);
			p.SetExpr(rateEquations[i]);
			
			for (int j=0; j < variableNames.size(); ++j)
				p.DefineVar(variableNames[j],&variableValues[j]);

			for (int j=0; j < parameterNames.size(); ++j)
				p.DefineVar(parameterNames[j],&parameterValues[j]);
			
			for (int j=0; j < assignmentVariables.size(); ++j)
				p.DefineVar(assignmentVariables[j],&assignmentValues[j]);

			try
			{
				rateEqns.push_back(p);
			}
			catch(...)
			{
				reactionNames.clear();
				rateEqns.clear();
				break;
			}
		}
		
		for (int i=0; i < assignmentEquations.size(); ++i)
		{
			mu::Parser p;
			addSBMLFunctions(p);
			p.SetExpr(assignmentEquations[i]);
			
			for (int j=0; j < variableNames.size(); ++j)
				p.DefineVar(variableNames[j],&variableValues[j]);

			for (int j=0; j < parameterNames.size(); ++j)
				p.DefineVar(parameterNames[j],&parameterValues[j]);
			
			for (int j=0; j < assignmentVariables.size(); ++j)
				p.DefineVar(assignmentVariables[j],&assignmentValues[j]);

			try
			{
				p.Eval();
				assignmentEqns.push_back(p);
			}
			catch(...)
			{
				//assignmentVariables.clear();
				//assignmentEqns.clear();
				break;
			}
		}

		for (int i=0; i < eventTriggers.size(); ++i)
		{
			mu::Parser p;
			addSBMLFunctions(p);
			p.SetExpr(eventTriggers[i]);
			
			for (int j=0; j < variableNames.size(); ++j)
				p.DefineVar(variableNames[j],&variableValues[j]);

			for (int j=0; j < parameterNames.size(); ++j)
				p.DefineVar(parameterNames[j],&parameterValues[j]);
			
			for (int j=0; j < assignmentVariables.size(); ++j)
				p.DefineVar(assignmentVariables[j],&assignmentValues[j]);

			try
			{
				p.Eval();
				
				//resposes for the trigger
				vector<mu::Parser> responses;
				for (int j=0; j < eventResponses[i].size(); ++j)
				{
					mu::Parser p;
					addSBMLFunctions(p);
					p.SetExpr(eventResponses[i][j]);

					try
					{
						p.Eval();
						responses.push_back(p);
					}
					catch(...) {}
				}

				if (responses.size() > 0)
				{
					triggerEqns.push_back(p);
					responseEqns.push_back(responses);
				}
			}
			catch(...)
			{
				//assignmentVariables.clear();
				//assignmentEqns.clear();
				break;
			}
		}

		//delete params;
		//delete reacs;
		delete doc;
	}
}

vector< vector<double> > SBML_sim::simulate(double time, double stepSize) const
{
	int n = variableValues.size();
	double * y0 = new double[n];
	for (int i=0; i < variableValues.size(); ++i)
		y0[i] = variableValues[i];
	
	double * y = ODEsim2(n, reactionNames.size(), stoichiometryMatrix , &sbml_rates_function, y0, 0.0, time, stepSize, (void*)this, triggerEqns.size() , sbml_event_function, sbml_response_function);
	
	vector< vector<double> > res;
	int sz = (int)(time/stepSize);
	
	if (y)
	{
		for (int j=0; j <= n; ++j)
		{
			vector<double> col(sz,0.0);
			for (int i=0; i < sz; ++i)
				col[i] = y[ i * (n+1) + j ];
			res.push_back(col);
		}
		free(y);
	}
	
	free(y0);
	return res;
}

vector<double> SBML_sim::steadyState() const
{
	int n = variableValues.size();
	double * y0 = new double[n];
	for (int i=0; i < variableValues.size(); ++i)
		y0[i] = variableValues[i];
	
	double * y = steadyState2(n, reactionNames.size(), stoichiometryMatrix , &sbml_rates_function, y0, (void*)this, 1.0E-5, 10000.0, 1.0, triggerEqns.size() , sbml_event_function, sbml_response_function);
	vector< double > res(n,0.0);
	
	if (y)
	{
		for (int j=0; j < n; ++j)
		{
			res[j] = y[j];
		}
		free(y);
	}
	
	free(y0);
	return res;	
}

vector< vector<double> > SBML_sim::ssa(double time) const
{
	int n = variableValues.size();
	double * y0 = new double[n];
	for (int i=0; i < variableValues.size(); ++i)
		y0[i] = variableValues[i];
		
	int sz;
	double * y = SSA(n, reactionNames.size(), stoichiometryMatrix , &sbml_rates_function, y0, 0.0, time, 100000, &sz, (void*)this, triggerEqns.size(), sbml_event_function, sbml_response_function);
	
	vector< vector<double> > res;
	
	if (y)
	{
		for (int j=0; j <= n; ++j)
		{
			vector<double> col(sz,0.0);
			for (int i=0; i < sz; ++i)
				col[i] = y[ i * (n+1) + j ];
			res.push_back(col);
		}
		free(y);
	}
	
	free(y0);
	return res;
}

void SBML_sim::setVariableValues( const vector<double> & v )
{
	for (int i=0; i < v.size() && i < variableValues.size(); ++i)
		variableValues[i] = v[i];
}

void SBML_sim::setParameters( const vector<double> & v )
{
	for (int i=0; i < v.size() && i < parameterValues.size(); ++i)
		parameterValues[i] = v[i];
}

vector< string > SBML_sim::getVariableNames() const
{
	return variableNames;
}

vector< string > SBML_sim::getParameterNames() const
{
	return parameterNames;
}

vector< double > SBML_sim::getVariableValues() const
{
	return variableValues;
}

void SBML_sim::reset()
{
	variableValues = initialValues;
	parameterValues = originalParameters;
}

vector< double > SBML_sim::getRateValues() const
{
	return rateValues;
}

vector< double > SBML_sim::getParameterValues() const
{
	return parameterValues;
}

/*****************************************************************
      GENETIC ALGORITHM BASED OPTIMIZATION -- HELPER FUNCTIONS
******************************************************************/

typedef GA1DArrayGenome<float> RealGenome;

static void initializeGenome(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;
	for (int i=0; i < g.size(); ++i)
		g.gene(i,0) = mtrand() * pow(10.0, 3.0*mtrand());
}

static float EuclideanDistance(const GAGenome & c1, const GAGenome & c2)
{
  const RealGenome & a = (RealGenome &)c1;
  const RealGenome & b = (RealGenome &)c2;

  float x=0.0;
  for(int i=0; i < b.length() && i < a.length(); ++i)
	  x += (a.gene(i) - b.gene(i))*(a.gene(i) - b.gene(i));

  return (float)sqrt(x);
}

static vector< vector<double> > actual;
static double end_time = 20.0;
static double dt = 0.1;

static float Objective(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;
	
	SBML_sim * model = (SBML_sim*)(g.geneticAlgorithm()->userData());
	
	vector<double> params(g.size(),0);
	for (int i=0; i < g.size(); ++i) params[i] = g.gene(i);

	model->reset();
	model->setParameters(params);

	vector< vector<double> > res = model->simulate(end_time, dt);
	
	if (res.size() < 1)
		return 100.0;

	double sumsq = 0.0, d = 0.0;
	
	for (int i=1; i < res.size(); ++i)
	{
		for (int j=0; j < res[i].size(); ++j)
		{
			d = res[i][j] - actual[i][j];
			sumsq += d*d;
		}
	}

	return (float)(sumsq/ (res[0].size()));
}

/********************************************
      GENETIC ALGORITHM BASED OPTIMIZATION
*********************************************/

vector< vector< double> > SBML_sim::optimize(const vector< vector<double> >& data, int iter)
{
	actual = data;
	vector< vector< double> > genes;

	if (actual.size() < 2 || actual[0].size() < 2) return genes;
	
	end_time = actual[0][ actual[0].size()-1 ];
	dt = actual[0][1] - actual[0][0];
	
	RealGenome genome( parameterValues.size() , &Objective);
	genome.initializer(&initializeGenome);
	GASteadyStateGA ga(genome);
	ga.userData(this);
	
	GASharing dist(EuclideanDistance);
	ga.scaling(dist);
	ga.pReplacement(1.0);
	ga.minimize();
	ga.populationSize(1000);
	ga.nGenerations(iter);
	ga.pMutation(0.2);
	ga.pCrossover(0.9);
	GAPopulation pop;
	ga.evolve();
	
	pop = ga.population();
	pop.order(GAPopulation::HIGH_IS_BEST);
	pop.sort(gaTrue);
	
	genes.resize(pop.size());
	
	for (int i=0; i < pop.size(); ++i)
	{
		RealGenome & g = (RealGenome &)(pop.individual(i));
		genes[i].resize(g.size());

		for (int j=0; j < g.size(); ++j)
			genes[i][j] = g.gene(j);
	}
	
	parameterValues = genes[0];

	return genes;
}

