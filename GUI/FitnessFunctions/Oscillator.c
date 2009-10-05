/* Fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(GAindividual net)
{
	int i, N;
	double * y, time, f, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0;
	
	N = getNumSpecies(net);
	
	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	f = 0;   // Calculate correlation to sine wave
	if (y != 0)
	{
		mXY = mX = mY = mX2 = mY2 = 0;
		
		for (i = 0; i < time; ++i)
		{
			mX += getValue(y,N+1,i,1);
			mY += sin(i/4.0);
			mXY += sin(i/4.0) * getValue(y,N+1,i,1);
			mX2 += getValue(y,N+1,i,1)*getValue(y,N+1,i,1);
			mY2 += sin(i/4.0)*sin(i/4.0);
		}
		mX /= time;
		mY /= time;
		mXY /= time;
		mX2 /= time;
		mY2 /= time;

		if (((mXY - mX*mY) < 0.01) || ((mY2 - mY*mY)) < 0.01)
			f = 0.0;
		else
		{
			f = ( (mXY - mX*mY)/(sqrt(mX2 - mX*mX)*sqrt(mY2 - mY*mY)) );   // Correlation formula
			if (f < 0) f = -f; // Negative correlation is just as good as positive (for oscillations)
		}
		free(y);
	}

	if(getNumSpecies(net) > 30)        // Disallow large networks
	  return (0);
	return (f);
}
