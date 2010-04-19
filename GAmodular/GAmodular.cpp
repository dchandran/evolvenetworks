#include <iostream>
#include <vector>
#include <math>
#include "sbml_sim.h"
extern "C"
{
	#include "mtrand.h"
	#include "opt.h"
}

using namespace std;

class ModularGA
{
private:

	struct Genome
	{
		Genome * parent1;
		Genome * parent2;

		SBML_sim * model;
		vector<double> params;
		double fitness;
		double score;

		Genome(): 
			model(0), 
			parent1(0), 
			parent2(0) 
		{}

		Genome(SBML_sim * m): 
			model(*m), 
			parent1(0), 
			parent2(0) 
		{
			params.resize(model->getParameterValues().size());
			for (int i=0; i < params.size(); ++i)
				params[i] = mtrand() * pow(10, mtrand() * 5);
		}

		Genome(const Genome & g) : 
			model(new SBML_sim(*g.model)), 
			params(g.params), 
			parent1(0), 
			parent2(0) 	
			{}

		~Genome()
		{
			delete model;
		}
	};

	Genome * crossover(Genome * parent1, Genome * parent2)
	{
		Genome * child = new Genome(*parent1);
		child->parent1 = parent1;
		child->parent2 = parent2;

		int n = parent1->params.size();
		int i = (int)(mtrand() * n);

		for (; i < n; ++i)
			child->params[i] = parent2->params[i];

		return child;
	}

	Genome * mutate(Genome * original, double pMut)
	{
		Genome duplicate = new Genome(*original);

		int n = original->params.size();
		int k = (int)(pMut * n);

		for (int i=0; i < k; ++i)
			duplicate->params[ (int)(mtrand() * n) ] *= (mtrand() * 2.0);

		return duplicate;
	}

	double distance(Genome * g1, Genome * g2)
	{
		double d = 0.0, sum = 0.0;
		for (int i=0; i < g1->params.size(); ++i)
		{
			d = (g1->params[i] - g2->params[i]);
			sum += d*d;
		}
		return sqrt(sum);
	}

	vector<Genome*> population;

	vector< vector<double> > dists;

	double (*calcScore)(SBML_sim*);

public:

	void objective( double (*f)(SBML_sim*) )
	{
		calcScore = f;
	}

	void run(SBML_sim * model, int populationSz, int generations = 100)
	{
		dists.resize(populationSz);
		for (int i=0; i < populationSz; ++i)
		{
			Genome * g = new Genome(model);
			dists[i].resize(populationSz);
		}
	}

	void oneStep()
	{
		//score

		double max_score, min_score, score;
		int i;

		for (i=0; i < population.size(); ++i)
		{
			population[i]->model->setParameters( population[i]->params );
			population[i]->score = score = calcScore(population[i]->model);

			if (i==0)
			{
				max_score = min_score = score;
			}
			else
			{
				if (score > max_score)
					max_score = score;
				if (score < min_score)
					min_score = score;
			}
		}
		
		double range = max_score - min_score;
		if (range <= 0) range = 1.0;

		//fitness (crowding)
		double max_dist, min_dist, dist;
		for (i=0; i < (population.size()-1); ++i)
		{
			for (int j=(i+1); j < population.size(); ++j)
			{
				dist = dists[j][i] = dists[i][j] = distance(population[i],population[j]);
				if (i==0 && j == 0)
				{
					max_dist = min_dist = dist;
				}
				else
				{
					if (dist > max_dist)
						max_dist = dist;
					if (dist < min_dist)
						min_dist = dist;
				}
			}
			dists[i][i] = 0.0;
		}
		dists[i][i] = 0.0;

		for (i=0; i < population.size(); ++i)
		{

		}
	}
};

