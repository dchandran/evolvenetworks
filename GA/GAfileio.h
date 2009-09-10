/********************************************************************************************************

Copyright (C) 2009 Deepak Chandran

This file contains several functions for automating several runs
using the reactionNetwork framework.
The included functions are:

1) repeatedRuns - this function will repeatedly run the GA on the same fitness
function and record relevant information

2) repeatedRunsFromFile - this function will load the information needed for 
repeatedRuns() from a correctly formatted file and call repeatedRuns, producing 
the output file

3) configureRepeatedRuns - used to configure what information should be printed in the output
file produced by repeatedRuns()


*********************************************************************************************************/

#include "reactionNetwork.h"

#ifndef GA_AUTOMATE_FUNCTIONS_H
#define GA_AUTOMATE_FUNCTIONS_H

/*! \brief Runs the same fitness function many times and outputs to a log file.
	The log file will contain the following information: seed for the random number
	generator (to reproduce results), best fitness and corresponding network, 
	mean and st.dev of final fitness values, 
		
	\param char* name of the output log file
	\param GAFitnessFnc fitness function
	\param int population size in each run
	\param int number of generations in each run, i.e each evolution experiment
	\param int number of repeated runs of the evolution experiment
	\param 
*/
void repeatedRuns(char * logfile, GAFitnessFnc fitness, int generations, int runs);

void repeatedRunsFromFile(char * intputFile);

void configureRepeatedRuns();

#endif
