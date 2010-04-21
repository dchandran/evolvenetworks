#include "cvodesim.h"

typedef struct 
{
	double a, b;
}
MyStruct;

void myODE(double time, double * y, double * dydt, void * myData)
{
	MyStruct * p = (MyStruct*)myData;
	
	dydt[0] = -p->a * y[1];
	dydt[1] =  p->b * y[0];
}

void myODE2(double time, double * y, double * dydt, void * myData)
{
	MyStruct * p = (MyStruct*)myData;
	
	dydt[0] = p->a * y[1] - p->b * y[0];
	dydt[1] = p->b * y[0] - p->a * y[1];
}

int events(int i, double time, double * y, void * myData)
{
	MyStruct * p = (MyStruct*)myData;

	if (i == 0)
		return (time >= 10.0 && time < 50.0 && p->a < 1.0);
	
	if (i == 1)
		return (time > 80.0 && p->a > 0.0);
}

void responses(int i, double * y, void * myData)
{
	MyStruct * p = (MyStruct*)myData;
	
	if (i==0)
	{
		p->a = 1.0;
		p->b = 1.0;		
	}

	if (i==1)
	{
		p->a = 0.0;
		p->b = 0.0;
	}
}

int main()
{
	int i,j,sz;
	
	int numVars = 2, 
		numEvents = 2;

	double initialValue[] = { 100.0, 0.0 };

	double startTime = 0.0, endTime = 100.0, stepSize = 0.1;
	
	MyStruct p = { 1.0, 1.0 };
	
	/*
	double * y = ODEsim( numVars,  initialValue, myODE, startTime, endTime, stepSize, (void*)(&p), numEvents, events, responses);

	if (y)
	{
		sz = (endTime - startTime)/stepSize;
		for (i=0; i < sz; ++i)
		{
			for (j=0; j < (1+numVars); ++j)
				printf("%lf\t",  y [ i*(numVars+1) + j ]);  //get i,j-th value
			printf("\n");
		}
	}
	else
	{
		printf("integration error\n");
	}*/
	
	double * ss;
	
	for (i=0; i < 10; ++i)
	{
		p.a = i;

		ss = steadyState( numVars,  initialValue,  myODE2, (void*)(&p), 1.0E-4, 1000, 1 , 0 , 0, 0);
		
		if (ss)
		{
			for (j=0; j < numVars; ++j)
				printf("%lf\t",  ss [j]);
			printf("\n");
			free(ss);
		}
	}
	
	return 0;
}

