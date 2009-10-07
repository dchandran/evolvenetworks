
/* Fitness function that tests for oscillations by counting the number of peaks*/
double fitness(GAindividual net)
{
	int i, N, peaks, troughs;
	double * y, time, f, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0, dx = 0.01;

	N = getNumSpecies(net);
	
	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	N = getNumSpecies(net);
	
	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	f = 0;   // Calculate correlation to sine wave
	if (y != 0)
	{
		peaks = 0;
		troughs = 0;
		for (i = 5; i < (time-3); ++i)
		{
			if ( (getValue(y,N+1,i,1) > 0.1) &&
				 (getValue(y,N+1,i,1) < 1000.0) &&
				 (getValue(y,N+1,i,1) > (dx + getValue(y,N+1,i-3,1))) &&
				 (getValue(y,N+1,i,1) > (dx + getValue(y,N+1,i+3,1))) &&
				 (getValue(y,N+1,i-3,1) < getValue(y,N+1,i-2,1)) && 
				 (getValue(y,N+1,i-2,1) < getValue(y,N+1,i-1,1)) && 
				 (getValue(y,N+1,i-1,1) < getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+1,1) < getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+2,1) < getValue(y,N+1,i+1,1)) && 
				 (getValue(y,N+1,i+3,1) < getValue(y,N+1,i+2,1))
				)
			{
				 ++peaks;
				 mX += y[i];
				 mX2 += y[i]*y[i];
			}

			if ( (getValue(y,N+1,i,1) > 0.1) &&
				 (getValue(y,N+1,i,1) < 1000.0) &&
				 (getValue(y,N+1,i,1) < (getValue(y,N+1,i-3,1) - dx)) &&
				 (getValue(y,N+1,i,1) < (getValue(y,N+1,i+3,1) - dx)) &&
				 (getValue(y,N+1,i-3,1) > getValue(y,N+1,i-2,1)) && 
				 (getValue(y,N+1,i-2,1) > getValue(y,N+1,i-1,1)) && 
				 (getValue(y,N+1,i-1,1) > getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+1,1) > getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+2,1) > getValue(y,N+1,i+1,1)) && 
				 (getValue(y,N+1,i+3,1) > getValue(y,N+1,i+2,1))
				)
			{
				 ++troughs;
			}
		}
		
		if ((troughs+peaks) > 30)
		{
			f = (double)30.0 + 1.0/(1.0 + mX2 - mX*mX);
		}
		else
		{
			mXY = mX = mY = mX2 = mY2 = 0;
			
			for (i = 10; i < time; ++i)
			{
				mX += getValue(y,N+1,i,1);
				mY += sin(i/4.0);
				mXY += sin(i/4.0) * getValue(y,N+1,i,1);
				mX2 += getValue(y,N+1,i,1)*getValue(y,N+1,i,1);
				mY2 += sin(i/4.0)*sin(i/4.0);
			}
			mX /= (time-10);
			mY /= (time-10);
			mXY /= (time-10);
			mX2 /= (time-10);
			mY2 /= (time-10);

			if (((mXY - mX*mY) < 0.01) || ((mY2 - mY*mY)) < 0.01)
			{
				f = 0.0;
			}
			else
			{
				f = ( (mXY - mX*mY)/(sqrt(mX2 - mX*mX)*sqrt(mY2 - mY*mY)) );   // Correlation formula
				if (f < 0) f = -f; // Negative correlation is just as good as positive (for oscillations)
			}
			f += (double)troughs + (double)peaks;
		}
		
		free(y);
	}

	return (f);
}
