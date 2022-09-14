/* Main function - example of how to use the CRKMethods solvers.
* Copyright (C) 2014  Georgios Is. Detorakis (gdetor@protonmail.com)
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>
*/

#include <stdio.h>
#include "solvers.h"


int main(int argc, char **argv)
{
	int i, N = 100; /* Number of nodes (discretization) */
	float t0 = 0.0, y0 = 0.5;  /* Initial conditions */
	float h = 1.0/N;   /* RK step */
    float time, y; /* y (solution) at time - time */

	FILE *fp; 

    /* Open a file to write in the results */
    /* One can store the results in vectors */
    fp = fopen("results.dat", "w");

	for(i = 0; i < N; i++){
		eulerStep(f, t0, y0, h, i, &time, &y); 
        /* rkFehlbergStep(f, fy, t0, y0, h, i, &time, &y); */
        fprintf(fp, "%f  %f\n", time, y);
        t0 = time;
        y0 = y;
	}

    fclose(fp);

	return 0;
}
