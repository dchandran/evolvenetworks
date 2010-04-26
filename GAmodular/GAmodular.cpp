#include "GAmodular.h"

using namespace std;

/*******************
  Helper functions
********************/
int partitionLT(vector<Genome*>& list, int left,int right,int pivotIndex)
{
	Genome * pivotValue = list[pivotIndex];
	Genome * temp = list[pivotIndex];
	list[pivotIndex] = list[right];
	list[right] = temp;
	int storeIndex = left;
	
	for (int i=left; i < right; ++i)
	{
		if (list[i]->fitness < pivotValue->fitness)
		{
			temp = list[storeIndex];
			list[storeIndex] = list[i];
			list[i] = temp;
			storeIndex = storeIndex + 1;
		}
	}

	temp = list[right];
	list[right] = list[storeIndex];  // Move pivot to its final place
	list[storeIndex] = temp;
	return storeIndex;
}

int partitionGT(vector<Genome*>& list, int left,int right,int pivotIndex)
{
	Genome * pivotValue = list[pivotIndex];
	Genome * temp = list[pivotIndex];
	list[pivotIndex] = list[right];
	list[right] = temp;
	int storeIndex = left;
	
	for (int i=left; i < right; ++i)
	{
		if (list[i]->fitness < pivotValue->fitness)
		{
			temp = list[storeIndex];
			list[storeIndex] = list[i];
			list[i] = temp;
			storeIndex = storeIndex + 1;
		}
	}
	
	temp = list[right];
	list[right] = list[storeIndex];  // Move pivot to its final place
	list[storeIndex] = temp;
	return storeIndex;
}

void findSmallestK(vector<Genome*>& list, int left,int right,int k)
{
	while (right > left)
	{
		int pivotIndex = (int)((right + left)/2);
		int pivotNewIndex = partitionLT(list, left, right, pivotIndex);
		if (pivotNewIndex > k)
			right = pivotNewIndex-1;
		else
		if (pivotNewIndex < k)
		{
			left = pivotNewIndex+1; 
			k = k-pivotNewIndex;
		}
	}
}

void findLargestK(vector<Genome*>& list, int left,int right,int k)
{
	while (right > left)
	{
		int pivotIndex = (int)((right + left)/2);
		int pivotNewIndex = partitionGT(list, left, right, pivotIndex);
		if (pivotNewIndex > k)
			right = pivotNewIndex-1;
		else
		if (pivotNewIndex < k)
		{
			left = pivotNewIndex+1; 
			k = k-pivotNewIndex;
		}
	}
}

/********************
    class Genome
*********************/

Genome::Genome(): 
	model(0), 
	parent1(0), 
	parent2(0),
	fitness(0),
	score(0)
{}

Genome::Genome(const Genome & g) : 
	model(0), 
	params(g.params), 
	parent1(g.parent1),
	parent2(g.parent2),
	fitness(g.fitness),
	score(g.score)
{
	if (g.model)
		model = new SBML_sim(*g.model);
}

Genome::~Genome()
{
	if (model)
		delete model;
}

Genome * Genome::clone() const
{
	return new Genome(*this);
}

Genome * Genome::crossover(Genome * parent2, double pCross) const
{
	Genome * child = this->clone();
	child->parent1 = this;
	child->parent2 = parent2;

	int n = params.size();
	int i = n - (int)(pCross * mtrand() * n);

	for (; i < n; ++i)
		child->params[i] = parent2->params[i];

	return child;
}

void Genome::mutate(double pMut)
{
	int n = params.size();
	int k = (int)(pMut * n);

	for (int i=0; i < k; ++i)
		params[ (int)(mtrand() * n) ] *= (mtrand() * 2.0);
}

double Genome::distance(Genome * g) const
{
	double d = 0.0, sum = 0.0;
	for (int i=0; i < params.size(); ++i)
	{
		d = (params[i] - g->params[i]);
		sum += d*d;
	}
	return sqrt(sum);
}

/********************
   class ModularGA
*********************/

ModularGA::ModularGA(): 
	_calcScore(0), 
	_callback(0), 
	_selectionMethod(RouletteWheel), 
	_popSz(1000), 
	_numGen(100), 
	_neighborhood(0.25), 
	_objectiveType(Maximize),
	_selectionRate(0.9),
	_pMut(0.2),
	_pCross(0.5)
{}

void ModularGA::setPopulationSize(int sz)
{
	if (sz < 1) 
		_popSz = 1;
	else
		_popSz = sz;

	_distMatrix.resize(_popSz);
	_population.resize(_popSz);

	for (int i=0; i < _popSz; ++i)
	{
		_distMatrix[i].resize(_popSz);
		_population[i] = 0;
	}
}


void ModularGA::setParameterSize(int sz, double low, double high)
{
	_population.resize(_popSz);

	for (int i=0; i < _popSz; ++i)
	{
		_population[i] = new Genome;
		_population[i]->params.resize(sz);
		for (int j=0; j < sz; ++j)
			_population[i]->params[j] = low + mtrand() * pow(10.0, high*mtrand());
	}
}

void ModularGA::setGenerations(int n)
{
	if (n > 0)
		_numGen = n;
}

void ModularGA::setSelectionRate(double p)
{
	if (p < 1.0 && p > 0.0)
		_selectionRate = p;
}

void ModularGA::setMutationRate(double p)
{
	if (p < 1.0 && p > 0.0)
		_pMut = p;
}

void ModularGA::setCrossoverRate(double p)
{
	if (p < 1.0 && p > 0.0)
		_pCross = p;
}


void ModularGA::setObjective( double (*f)(Genome*) , ObjectiveType type )
{
	_objectiveType = type;
	_calcScore = f;
}

void ModularGA::setCallback( bool (*f)() )
{
	_callback = f;
}

vector<Genome*> ModularGA::population() const
{
	return _population;
}

Genome* ModularGA::best() const
{
	Genome * g = 0;
	
	for (int i=0; i < _population.size(); ++i)
		if (!g || 
			(_objectiveType == Minimize && _population[i]->score < g->score) ||
			(_objectiveType == Maximize && _population[i]->score > g->score))
			
			g = _population[i];
	
	return g;
}

Genome* ModularGA::worst() const
{
	Genome * g = 0;
	
	for (int i=0; i < _population.size(); ++i)
		if (!g || 
			(_objectiveType == Minimize && _population[i]->score > g->score) ||
			(_objectiveType == Maximize && _population[i]->score < g->score))
			
			g = _population[i];
	
	return g;
}

double ModularGA::evolve()
{
	for (int i=0; i < _numGen; ++i)
	{
		oneStep();
		if (_callback && !_callback())
			break;
	}
	
	double bestScore = 0.0;
	for (int i=0; i < _popSz; ++i)
	{
		if (i==0)
			bestScore = _population[i]->score;
		else
		{
			if ((_objectiveType == Minimize && bestScore > _population[i]->score) ||
				(_objectiveType == Maximize && bestScore < _population[i]->score))
				
				bestScore = _population[i]->score;
		}		
	}
	return bestScore;
}

void ModularGA::oneStep()
{
	//score

	double max_score, min_score, score, max_dist, min_dist, dist, range_score, range_dist;
	int i,j,k,n;

	for (i=0; i < _population.size(); ++i)
	{
		_population[i]->score = score = _calcScore( _population[i] );

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
	
	range_score = max_score - min_score;
	if (range_score <= 0) range_score = 1.0;

	//distances
	
	for (i=0; i < (_popSz-1); ++i)
	{
		for (j=(i+1); j < _popSz; ++j)
		{
			dist = _distMatrix[j][i] = _distMatrix[i][j] = _population[i]->distance(_population[j]);

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
		_distMatrix[i][i] = 0.0;
	}
	_distMatrix[i][i] = 0.0;
	
	range_dist = max_dist - min_dist;
	if (range_dist <= 0) range_dist = 1.0;

	//compute fitness
	
	double sum_dist = 0.0, 
		   sum_fitness = 0.0;
	
	for (i=0; i < _popSz; ++i)
	{
		sum_dist = 1.0;
		
		for (j=0; j < _popSz; ++j)
		{
			dist = (_distMatrix[i][j] - min_dist)/range_dist;
			if (dist < _neighborhood)
			{
				sum_dist += 1.0 - dist/_neighborhood;
			}
		}
	
		_population[i]->fitness =  ((_population[i]->score - min_score)/range_score);///(sum_dist);
		sum_fitness += _population[i]->fitness;
	}
	
	//select parents
	
	vector<Genome*> newPopulation;
	vector<bool> isSelected;
	
	newPopulation.resize(_popSz);
	isSelected.resize(_popSz);
	
	for (i=0; i < _popSz; ++i)
		isSelected[i] = false;
	
	n = (int)(_selectionRate * _popSz);
	
	if (_selectionMethod == RouletteWheel)
	{
		double r, sum;
		
		for (i=0; i < n; ++i)
		{
			r = mtrand();
			sum = 0.0;
			
			for (j=0; j < _popSz; ++j)
			{
				if (_objectiveType == Minimize)
					sum += (1.0 - _population[j]->fitness/sum_fitness);
				else
					sum += _population[j]->fitness/sum_fitness;
				if (r < sum)
					break;
			}
			
			if (isSelected[j])
			{
				newPopulation[i] = _population[j]->clone();
			}	
			else
			{
				newPopulation[i] = _population[j];
				isSelected[j] = true;
			}
		}
	}
	else
	if (_selectionMethod == Tournament)
	{
		double r1, r2;
		
		for (i=0; i < n; ++i)
		{
			r1 = (int)(mtrand() * _popSz);
			r2 = (int)(mtrand() * _popSz);
			
			if (_objectiveType == Minimize && _population[r1]->fitness < _population[r2]->fitness)
			{
				if (isSelected[r1])
				{
					newPopulation[i] = _population[r1]->clone();
				}
				else
				{
					newPopulation[i] = _population[r1];
					isSelected[r1] = true;
				}
			}
			else
			{
				if (isSelected[r2])
				{
					newPopulation[i] = _population[r2]->clone();
				}
				else
				{
					newPopulation[i] = _population[r2];
					isSelected[r2] = true;
				}
			}
		}
	}
	else
	if (_selectionMethod == Elitism)
	{
		if (_objectiveType == Minimize)
			findSmallestK(_population, 0, _popSz, n);
		else
			findLargestK(_population, 0, _popSz, n);
		for (i=0; i < n; ++i)
		{
			newPopulation[i] = _population[i];
			isSelected[i] = true;
		}
	}
	
	//generate children
	for (i=n; i < _popSz; ++i)
	{
		j = (int)(mtrand() * n);
		k = (int)(mtrand() * n);
		newPopulation[i] = newPopulation[j]->crossover(newPopulation[k]);
	}
	
	//mutate population
	for (i=0; i < n; ++i)
		newPopulation[i]->mutate(_pMut);
	
	//remove old population
	for (i=0; i < _popSz; ++i)
		if (!isSelected[i])
			delete _population[i];

	//update
	_population = newPopulation;	
}

