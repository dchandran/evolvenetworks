int FITNESS_NUM = 0;

double** XORtable() //make XOR table
{
	int i,j;
	//XOR logic
	double XOR[][5] =
	{
	   {  0.1,   0.1,   0.1,  0.1 ,   0.0,  0.0 },  
	   {  0.1,   0.1,   0.1, 10.0 ,   0.0, 10.0 }, 
	   {  0.1,   0.1, 10.0, 0.1  ,   0.0,  10.0}, 
	   {  0.1,   0.1, 10.0, 10.0 ,   0.0, 0.0},  
       {  0.1, 10.0,   0.1, 0.1 ,   0.0 , 10.0}, 
	   {  0.1, 10.0,   0.1, 10.0 ,   10.0, 10.0},
	   {  0.1, 10.0,  10.0, 0.1,   10.0, 10.0},
	   {  0.1, 10.0,  10.0, 10.0,   0.0, 10.0}, 
       {  10.0,  0.1,   0.1, 0.1,   0.0 , 10.0}, 
	   {  10.0,  0.1,   0.1, 10.0 ,   10.0 ,  10.0}, 
	   {  10.0,  0.1,  10.0, 0.1 ,   10.0, 10.0}, 
	   {  10.0,  0.1,  10.0, 10.0 ,   0.0, 10.0},  
       {  10.0, 10.0,   0.1, 0.1 ,   0.0, 0.0},  
	   {  10.0, 10.0,   0.1, 10.0 ,   0.0, 10.0}, 
	   {  10.0, 10.0,  10.0, 0.1 ,   0.0, 10.0}, 
	   {  10.0, 10.0,  10.0, 10.0 ,   0.0, 0.0}  
	};
	
	double ** table = (double**)malloc( 16 * sizeof(double*) );
	
	for (i=0; i < 16; ++i)
	{
		table[i] = (double*)malloc( 5 * sizeof(double) );
		for (j=0; j < 4; ++j)
			table[i][j] = XOR[i][j];

		if  (FITNESS_NUM == 0)
			table[i][4] = XOR[i][4];
		else
			table[i][4] = XOR[i][5];
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
	
	num_rows = 16;
	num_inputs = 4;
	num_outputs = 1;
	
	score = compareSteadyStates(p, table, num_rows, num_inputs, num_outputs, 1, 0);
	
	//free table
	for (i=0; i < 16; ++i)
	{
		free(table[i]);
	}
	
	free(table);
	
	return score;
}

int COUNT = 0;
/* print the number of each generation and the fitness of the best network */
int callback(int iter,int popSz,GApopulation pop,double * fitnessArray, int *** parents)
{
	if (COUNT > 10) 
	{
		if (FITNESS_NUM == 0)
			FITNESS_NUM = 1;
		else
			FITNESS_NUM = 0;
		COUNT = 0;
	}
	++COUNT;
	return 0;
}

