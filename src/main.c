#include <stdio.h>
#include "cRKSolvers.h"

int main(int argc, char **argv)
{
	int i, N = 100; /* Number of nodes (discretization) */
	double t0 = 0.0, y0 = 0.5;  /* Initial conditions */
	double h = 1.0/N;   /* RK step */
    double time, y; /* y (solution) at time - time */

	FILE *fp; 

    /* Open a file to write in the results */
    /* One can store the results in vectors */
    fp = fopen( "results.dat", "w" );

	for(i = 0; i < N; i++){
		eulerStep(f,t0,y0,h,i,&time,&y); 
        /* rkFehlbergStep(f, fy, t0, y0, h, i, &time, &y); */
        fprintf( fp, "%f  %f\n", time, y );
        t0 = time;
        y0 = y;
	}

    fclose(fp);

	return 0;
}
