/*************************
This test code finds the parameters to force two of the ODE functions below to exhibit bistability

Generate the CVode library:
gcc -o cvode/*.*
ar *.o -o libcvode.a

Run this code:
gcc cvodesim.c mat.c neldermead.c ga.c mtrand.c ga_bistable.c test_bistable.c -I./ -L./ -lcvode
./a.out


**************************/
#include "ga_bistable.h"

void ode1(double time,double * u,double * du,void * data)
{
   Parameters * p = (Parameters*)data;
   double * k = (*p).params;
   double * a = (*p).alphas;
   double  r0 = k[2] - k[3]*u[0],  //inflow/outflow of s0
           r1 = k[4] - k[5]*u[1],  //inflow/outflow of s1
           r2 = k[0]*u[0]*u[0]*u[1] - k[1]*u[0]*u[0]*u[0]; //2s0 + s1 <=> 3s0
   du[0] = a[0]*(r0 + r2);
   du[1] = a[1]*(r1 - r2);
}

void ode2(double time,double * u,double * du,void * data)
{
   Parameters * p = (Parameters*)data;
   double * k = (*p).params;
   double * a = (*p).alphas;
   double  r0 = k[0]/(k[1] + pow(u[1],4)) - k[2]*u[0],  //s0
           r1 = k[3]/(k[4] + pow(u[0],4)) - k[5]*u[1];  //s1
   du[0] = a[0]*r0;
   du[1] = a[1]*r1;
}

int main()
{
   int i;
   double iv[] = { 1.8, 2.3 };
   Parameters * p;
   BistablePoint bis;


   bis = makeBistable(2,6,iv,30,2000,&(ode1));

   p = bis.param;
   if (!p) return 0;
   if (bis.unstable)
   {
      printf("\nunstable steady state:   ");
      for (i=0; i < (*p).numVars; ++i)
          printf("%lf ",bis.unstable[i]);

      free(bis.unstable);
   }
   else
   {
      printf("no unstable state\n");
   }
   
   if (bis.stable1 || bis.stable2)
   {
		if (bis.stable1)
		{
			printf("\nstable steady state:   ");
			for (i=0; i < (*p).numVars; ++i)
				printf("%lf ",bis.stable1[i]);

			free(bis.stable1);
		}
		if (bis.stable2)
		{
			printf("\nstable steady state:   ");
			for (i=0; i < (*p).numVars; ++i)
				printf("%lf ",bis.stable2[i]);

			free(bis.stable2);
		}
   }
   else
   {
      printf("no stable states\n");
   }
   
   printf("\nparameters: ");
   for (i=0; i < (*p).numParams; ++i)
   {
       printf("%lf ",(*p).params[i]);
   }
   printf("\n");
   printf("\nalphas: ");
   for (i=0; i < (*p).numVars; ++i)
   {
       printf("%lf ",(*p).alphas[i]);
   }
   printf("\n");
   
   deleteGAindividual(p);
   
   return 0;
}
