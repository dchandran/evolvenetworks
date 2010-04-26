#ifndef GA_MODULAR_H
#define GA_MODULAR_H

#include <iostream>
#include <vector>
#include <math.h>
#include "sbml_sim.h"
extern "C"
{
	#include "mtrand.h"
	#include "opt.h"
}

class Genome
{
public:
	SBML_sim * model;
	std::vector<double> params;

	const Genome * parent1;
	const Genome * parent2;

	double fitness; //scaled
	double score; //raw
	
	Genome();
	~Genome();
	Genome(const Genome&);
	
	virtual	Genome * clone() const;	
	virtual Genome * crossover(Genome*,double frac = 0.5) const;
	virtual void mutate(double frac = 0.5);	
	virtual double distance(Genome*) const;
};

class ModularGA
{

public:
	ModularGA();
	
	/*! \brief selection methods*/
	enum SelectionMethod { RouletteWheel=0, Tournament=1, Elitism=2 };
	
	/*! \brief type of objective*/
	enum ObjectiveType { Minimize=0, Maximize=1 };
	
	/*! \brief set the objective function*/
	void setObjective( double (*f)(Genome*) , ObjectiveType type = Maximize);

	/*! \brief set the population size*/
	void setPopulationSize(int);

	/*! \brief set the number of parameters*/
	void setParameterSize(int, double lowerBoundLog=0.0, double upperBoundLog=3.0);

	/*! \brief set the number of generations to run*/
	void setGenerations(int);
	
	/*! \brief set the distance (<1) that is defined as the neighborhood for computing crowding*/
	void setNeighborhood(double);

	/*! \brief evolve networks satisfying the objective*/
	double evolve();
	
	/*! \brief generate the next generation*/
	void oneStep();
	
	/*! \brief how the parents in one generation are selected for producing offsprings for the next generation*/
	void setSelectionMethod(SelectionMethod);
	
	/*! \brief what fraction on the parent population are preseved in the next generation*/
	void setSelectionRate(double);
	
	/*! \brief the fraction of a genome that gets mutated*/
	void setMutationRate(double);
	
	/*! \brief the fraction of a genome that gets replaced during crossover*/
	void setCrossoverRate(double);
	
	/*! \brief Set the callback function that will get called at the end of each generation. 
	           The callback function can return false to stop the evolution*/
	void setCallback( bool (*f)() );
	
	/*! \brief get all the genomes (unsorted)*/
	std::vector<Genome*> population() const;

	/*! \brief get the best genome*/
	Genome* best() const;
	
	/*! \brief get the worst genome*/
	Genome* worst() const;
	
protected:

	/*! \brief objective*/
	double (*_calcScore)(Genome*);
	
	/*! \brief callback*/
	bool (*_callback)();

	/*! \brief selection method*/
	SelectionMethod _selectionMethod;
	
	int _popSz;
	int _numGen;
	double _neighborhood;
	double _selectionRate;
	double _pMut;
	double _pCross;
	
	ObjectiveType _objectiveType;

	/*! \brief distance between each pair of genomes*/
	std::vector< std::vector<double> > _distMatrix;
	
	/*! \brief all the genomes*/
	std::vector< Genome* > _population;
};

#endif

