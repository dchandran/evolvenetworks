double** XORtable() //make XOR table
{
	int i,j;
	//XOR logic
	double XOR[][3] =
	{
	   {  0.1,  0.1,  0.0 },  //input = low low,   output = low
	   {  0.1, 10.0, 10.0 },  //input = low high,  output = high
	   { 10.0,  0.1, 10.0 },  //input = high low,  output = high
	   { 10.0, 10.0,  0.0 },  //input = high high, output = low
	};
	
	double ** table = (double**)malloc( 4 * sizeof(double*) );
	
	for (i=0; i < 4; ++i)
	{
		table[i] = (double*)malloc( 3 * sizeof(double) );
		for (j=0; j < 3; ++j)
			table[i][j] = XOR[i][j];
	}
	
	return table;
}

/* 
fitness function that tests how well the steady state outputs match the
given input/output table
*/
double fitness(GAindividual p)
{
	int i, num_rows, num_inputs, num_outputs;
	double score;
	
	double ** table = XORtable();
	
	num_rows = 4;
	num_inputs = 2;
	num_outputs = 1;
	
	score = compareSteadyStates(p, table, num_rows, num_inputs, num_outputs, 1, 0);
	
	//free table
	for (i=0; i < 4; ++i)
	{
		free(table[i]);
	}
	
	free(table);
	
	return score;
}

/* print the number of each generation and the fitness of the best network */
int callback(int iter,int popSz,GApopulation pop,double * fitnessArray, int *** parents)
{
	double f = fitnessArray[0];

	if (f > 0.9)
	{
		fitness(pop[0]);
		return 1;
	}
	return (int)(f > 0.9);
}

