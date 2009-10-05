
/* fitness that calculates the coefficient of variation (CV) */
double fitness(GAindividual p)
{
	int i,r,n,sz;
	double f, sd, dt, time, * y, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0;

	n = getNumSpecies(p);
	r = getNumReactions(p);

	time = 500.0;

	y = simulateNetworkStochastically(p,time,&sz);  //stochastic simulation

	f = 0;
	if (y != 0)         //compute the variance
	{
		mXY = mX = mY = mX2 = mY2 = 0;
		for (i = 0; i < (sz-1); ++i)
		{
			dt = getValue(y,n+1,i+1,0) - getValue(y,n+1,i,0);
			mX += getValue(y,n+1,i,1) * dt;
			mX2 += getValue(y,n+1,i,1)*getValue(y,n+1,i,1)*dt;
		}

		mX /= time;
		mX2 /= time;

		sd = sqrt(mX2 - mX*mX);  //standard deviation

		if (sd <= 0 || mX <= 0 || mX > 5.0)
			f = 0.0;
		else
			f = mX / sd;   // CV = sdev/mean, but the fitness = 1/CV = mean/sdev

		free(y);
	}

	if(getNumSpecies(p) > 5)       //disallow large networks
		f = 0.0;

	return (f);
}



