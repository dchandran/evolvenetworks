#include <fstream>
#include <iostream>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>

extern "C"
{
#include "mtrand.h"
}

#include "ParameterMC.h"

static int NUM_THREADS = 8;
static boost::mutex mutex;
std::ofstream fout;

static double ** SUM_XiXj = 0;
static double * SUM_Xi = 0;

void setNumThreads(int n)
{
	if (n > 0)
		NUM_THREADS = n;
}

static void sampleParameterSpaceSingleRun(int numSamples, int numParameters, RandomParameterFunc getRandParams, TargetFunc isGood, int * count)
{
	int i=0,j=0,k=0,n=0,accept=0;
	double * params = (double*)malloc(numParameters * sizeof(double));

	//run n times
	for (k=0; k < numSamples; ++k)
	{
		//get a random set of parameter
		getRandParams(params);

		//check if it is acceptable (the time consuming step, usually)
		accept = isGood(params);

		if (accept)
		{
			//if acceptable, store it
			mutex.lock();
			++(*count);
			for (i=0; i < numParameters; ++i)
			{
				std::cout << params[i] << "\t";
				SUM_Xi[i] += params[i];
				for (j=i; j < numParameters; ++j)
					SUM_XiXj[i][j] += params[i] * params[j];
			}
			std::cout << std::endl;
			mutex.unlock();
		}
	}

	delete params;
}

void sampleParameterSpace(int numSamples, int numParameters, RandomParameterFunc getRandParams, TargetFunc isGood, double * Mu, double ** Sigma)
{
	int i = 0, j = 0, n = 0, n_last = 0, num_threads = 0, count = 0;
	boost::thread ** threads = 0;
	boost::thread * thread = 0;

	/**** global arrays for storing covariance and means ****/
	//SUM_XiXj = new double*[numParameters];
	//SUM_Xi = new double[numParameters];

	SUM_XiXj = Sigma;
	SUM_Xi = Mu;

	for (i=0; i < numParameters; ++i)
	{
		//SUM_XiXj[i] = new double[numParameters];
		SUM_Xi[i] = 0.0;

		for (j=0; j < numParameters; ++j)
			SUM_XiXj[i][j] = 0.0;
	}

	/**** divide the task amongst n threads ****/

	if (numSamples < NUM_THREADS)
		numSamples = NUM_THREADS;

	n = (int)(numSamples/NUM_THREADS); //each thread is allocated n runs
	num_threads = NUM_THREADS;

	n_last = numSamples - (n * NUM_THREADS);  //last thread is given n_last runs
	if (n_last > 0)
		++num_threads;

	threads = new boost::thread*[num_threads];

	/**** make Boost threads ****/

	for (i=0; i < num_threads; ++i)
	{
		if ((n_last > 0) && (i == (num_threads-1)))
			n = n_last;

		thread = new boost::thread(
			boost::bind(
			&sampleParameterSpaceSingleRun,
			n,
			numParameters,
			getRandParams,
			isGood,
			&count));

		threads[i] = thread;
	}

	for (i=0; i < num_threads; ++i)
		threads[i]->join();


	for (i=0; i < num_threads; ++i)
		delete threads[i];

	delete threads;

	--count; //unbiased mean and cov

	for (i = 0; i < numParameters; ++i)
		SUM_Xi[i] /= (double)count;

	for (i = 0; i < numParameters; ++i)
	{
		for (j = i; j < numParameters; ++j)
		{
			SUM_XiXj[i][j] = SUM_XiXj[i][j]/(double)count - SUM_Xi[i]*SUM_Xi[j];
			SUM_XiXj[j][i] = SUM_XiXj[i][j];
		}
	}
}

//testing

int isgood(double * p)
{
	return (int)( ((p[0] - p[1]) > 1.0) || ((p[2] - p[3]) > 1.0) );
}

void randomParams(double * p)
{
	p[0] = 2.0*mtrand();
	p[1] = 2.0*mtrand();
	p[2] = 2.0*mtrand();
	p[3] = 2.0*mtrand();
}

int main(int argc, char* argv[])
{
	double * Mu = new double[4];
	double ** Sigma = new double*[4];

	for (int i=0; i < 4; ++i)
		Sigma[i] = new double[4];

	std::ofstream fout;

	fout.open("params.tab");	

	sampleParameterSpace(1000,4,&randomParams,&isgood,Mu,Sigma);

	printf("\n\nMu\n\n");

	for (int i=0; i < 4; ++i)
	{
		printf("%lf\t",Mu[i]);
	}

	printf("\n\nSigma\n\n");

	for (int i=0; i < 4; ++i)
	{
		for (int j=0; j < 4; ++j)
			printf("%lf\t",Sigma[i][j]);
		printf("\n");
	}

	fout.close();

	for (int i=0; i < 4; ++i)
		delete Sigma[i];

	delete Sigma;
	delete Mu;

	return 0;
}
