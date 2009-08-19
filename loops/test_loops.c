#include <stdlib.h>
#include <stdio.h>
#include "loops.h"

int main()
{
	int i,j;
	
	double J[] = { 
			-0.1 , 	0.5, 	0.0, 0.0,
			-0.56 ,	0.0, 	0.2, 0.0,
			0.0 , 	-0.2, 	0.0, 0.12,
			1.0 , 	0.0, 	0.0, 0.6 };
		
	LoopsInformation info = getLoops(J,2);

	for (i=0; i < info.numLoops; ++i)
	{
		for (j=0; j < info.loopLengths[i]; ++j)
		{
			printf("%i\t",info.nodes[i][j]);
		}
		printf("type = %i\n",info.loopTypes[i]);
	}

	freeLoopsInfo(info);
}
