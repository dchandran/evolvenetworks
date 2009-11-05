
double fitness(GAindividual net)
{
	int i, N, n;
	double x, * y, time, f, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0, dx = 0.01;

	N = getNumSpecies(net);
	
	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	f = 0;   // Calculate correlation to sine wave
	if (y != 0)
	{
		n = 0;
		mXY = mX = mY = mX2 = mY2 = 0;
		for (i = 0; i < (time); ++i)
		{
			x = sin((double)i/4.0) > 0;
			mX += getValue(y,N+1,i,1);
			mY += x;
			mXY += x * getValue(y,N+1,i,1);
			mX2 += getValue(y,N+1,i,1)*getValue(y,N+1,i,1);
			mY2 += x*x;
			++n;
		}
		
		mX /= (double)(n);
		mY /= (double)(n);
		mXY /= (double)(n);
		mX2 /= (double)(n);
		mY2 /= (double)(n);

		if (((mX2 - mX*mX)) < 0.0001)
		{
			f = 0.0;
		}
		else
		{
			f = ( (mXY - mX*mY)/(sqrt(mX2 - mX*mX)*sqrt(mY2 - mY*mY)) );   // Correlation formula
			if (f < 0) f = -f; // Negative correlation is just as good as positive (for oscillations)
		}
		free(y);
	}

	return (f);
}

int callback(int iter,int popSz, GApopulation P, double * fitnessArray, int *** parents)
{
	return (fitnessArray[0] > 0.6);
}

double maxfitness()
{
	return 1.0;
}
