#include <stdio.h>
#include "ParameterInformationContent.h"

double * parameterMI(OutputFunction F, int N, double * low, double * hi)
{
	int i,j,k,l,iter,iter2, blocks = 100;
	double * params = (double*)malloc(N * sizeof(double)),
			* MI = (double*)malloc(N * sizeof(double)),
			* temp, * temp2;
	double min, max, x, blocksz, u, H, H2, H3;
	
	//get range of x
	iter = 10000;
	temp = (double*)malloc(iter * sizeof(double));
	for (i=0; i < iter; ++i)
	{
		for (j=0; j < N; ++j)
			params[j] = mtrand() * (hi[j] - low[j]) + low[j];
		x = F(params);
		temp[i] = x;
		if (i == 0)
			max = min = x;
		else
			if (x < min) min = x;
			else
				if (x > max) max = x;
	}
	
	blocksz = (max-min)/(double)blocks;
	temp2 = (double*)malloc(blocks * sizeof(double));
	for (i=0; i < blocks; ++i)
		temp2[i] = 1E-5;

	for (i=0; i < iter; ++i)
		temp2[ (int)(temp[i]/blocksz) ] += 1.0;

	for (i=0; i < blocks; ++i)
		temp2[i] /= iter;
	
	free(temp);
	H = 0.0;
	for (i=0; i < blocks; ++i)
		H -= temp2[i] * log(temp2[i]);
	
	printf("H = %lf\n",H);
	free(temp2);
	iter = 1000;
	iter2 = 100;
	
	temp = (double*)malloc(blocks * sizeof(double));
	for (i=0; i < N; ++i)
	{
		MI[i] = 0.0;
		H3 = 0.0;
		
		for (j=0; j < iter2; ++j)
		{
			u = mtrand() * (hi[i] - low[i]) + low[i];
			for (k=0; k < blocks; ++k)
				temp[k] = 1E-5;
			
			for (k=0; k < iter; ++k)
			{
				for (l=0; l < N; ++l)
					params[l] = mtrand() * (hi[l] - low[l]) + low[l];
				params[i] = u;
				x = F(params);
				if (x > max) x = max;
				temp[ (int)(x/blocksz) ] += 1.0;
			}
			
			for (k=0; k < blocks; ++k)
				temp[k] /= iter;
			
			H2 = 0.0;
			for (k=0; k < blocks; ++k)
				H2 -= temp[k] * log(temp[k]);

			H3 += H2/iter2;
		}
		
		printf("H(Y | p%i) = %lf\n",i,H3);
		
		MI[i] = H - H3;
	}
	free(temp);
	
	return MI;
}

void ode(double t, double * u, double * du, void * p)
{
	double * k = (double*)p;
	double s0 = k[0], 
		   vmax1 = k[1], 
		   vmax2 = k[2],
		   ka1 = k[3],
		   ka2 = k[4],
		   ka3 = k[5],
		   h1 = k[6],
		   h2 = k[7],
		   x,
		   y;
	s0 = 1.0;
	vmax2 = 10.0;
	x = ka1*pow(s0,h1);
	y = ka2*pow(u[0],h2);
	du[0] = vmax1 * x/ (1 + x) - 0.5*u[0];
	du[1] = vmax2 * x/ (1 + x + y)  - 0.5*u[1];
}

double ss(double * p)
{
	int i;
	double x0[] = {1.0,1.0};
	double y3 = 0.0;
	double * y = steadyState(2, x0, ode, p, 0.01, 100.0, 2);
	
	if (y)
	{
		y3 = y[1];
	}
	else
	{
		for (i=0; i < 8; ++i) printf("%lf ",p[i]);
		printf("\n");
	}
	free(y);
	return y3;
}

int main()
{
	double low[] = {0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 1.0};
	double hi[] = {20.0, 20.0, 20.0, 100.0, 100.0, 100.0, 8.0, 8.0};
	double * p;
	
	initMTrand();

	int i,n;
	
	FILE * file;
	file = fopen("out2.txt","w");
	
	for (n=0; n < 5; ++n)
	{
		p = parameterMI(&ss,8,low,hi);
		for (i=0; i < 8; ++i)
			fprintf(file,"%lf\t",p[i]);
		fprintf(file,"\n");
		free(p);
	}
	
	fclose(file);
	
	return 0;
}
