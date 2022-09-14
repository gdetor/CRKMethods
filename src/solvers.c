/* RK Solvers - This file contains one-step solvers of some of the most basic
 * Runge-Kutta ODE solvers.
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
#include "solvers.h"


/* RHS of the ODE. Here you just have to declare your function f. */
REAL f(REAL t, REAL y) {
	return -y;
}

/* First derivative of function f. */ 
REAL fy(REAL t, REAL y) {
	return 1.0;
}


/*
 * This function implements a single step of Euler's Method. f is the 
 * right-hand side (RHS) of the ODE, t0 is the initial time, y0 is the initial 
 * value of the ODE, h is the step size of the method and i is the step 
 * counter. t is the time nodes on which the computation is performed and 
 * y is the solution. 
*/
void eulerStep(REAL (*f)(REAL,REAL),
               REAL t0,
               REAL y0,
               REAL h,
               int i,
               REAL *t,
               REAL *y) {
	*t = t0 + h * i;
	*y = y0 + h * (*f)(t0,y0);
}


/*
 * This function implements a single step of Runge-Kutta Method of order 4. 
 * f is the RHS side of the ODE, t0 is the initial time, y0 is the initial 
 * value of the ODE, h is the step size of the method and i is the step 
 * counter. t is the time nodes on which the computation takes place and y is
 * the solution.
*/
void rk4Step(REAL (*f)(REAL,REAL),
             REAL t0,
             REAL y0,
             REAL h,
             int i,
             REAL *t,
             REAL *y) {
	REAL k1, k2, k3, k4;

	k1 = y0 + h * (*f)(t0,y0);
	k2 = y0 + h * (*f)(t0 + h*.5, y0 + k1*.5);
	k3 = y0 + h * (*f)(t0 + h*.5, y0 + k2*.5);
	k4 = y0 + h * (*f)(t0 + h, y0 + k3 * h);

	*t = t0 + h * i;
	*y = y0 + ( k1 + 2. * k2 + 2. * k3 + k4 ) * 1.0/6.0;
}


/* 
 * This function implements a single step of Runge-Kutta Method of order 4, 
 * with different coefficients from the classic RK4. f is the RHS of the ODE, 
 * fy is the first derivative of the RHS of the ODE, t0 is the initial time, 
 * y0 is the initial value of the ODE, h is the step size of the method and i 
 * is the step counter. t is the time on which the computation is performed 
 * and y is the solution. 
*/
void rk4RefinedStep(REAL (*f)(REAL,REAL),
                    REAL (*fy)(REAL,REAL),
                    REAL t0,
                    REAL y0,
                    REAL h,
                    int i,
                    REAL *t,
                    REAL *y) {
	REAL k1, k2, k3, k4;

	k1 = h * f( t0, y0 );
	k2 = h * f( t0 + h/2, y0 + 1.0/3.0 * k1 + 1.0/18.0 * h * fy(t0,y0) * k1 );
	k3 = h * f( t0 + h/2, y0 - 152.0/125.0 * k1 + 252.0/125.0 * k2 - 44.0/125.0
            * h * fy(t0,y0) * k1 );
	k4 = h * f( t0 + h, y0 + 19.0/2.0 * k1 - 72.0/7.0 * k2 + 25.0/14.0 * k3 + 
            5.0/2.0 * h * fy(t0,y0) * k1 );
	
	*t = t0 + h * i;
	*y = y0 + 5.0/48.0 * k1 + 27.0/56.0 * k2 + 125.0/336.0 * k3 + 1.0/24.0 * k4;
}


/*
 * This function implements a single step of Runge-Kutta Method - Fehlberg. 
 * f is the RHS of the ODE, fy is the first derivative of f, t0 is the 
 * initial time, y0 is the initial value of the ODE, h is the step size of the
 * method, i is the step counter, t is the time points and y is the obtained 
 * solution on the corresponding time points. 
*/
void rkFehlbergStep(REAL (*f)(REAL,REAL),
                    REAL (*fy)(REAL,REAL),
                    REAL t0,
                    REAL y0,
                    REAL h,
                    int i,
                    REAL *t,
                    REAL *y) {
	REAL k1, k2, k3, k4, k5, k6;

	k1 = h * (*f)(t0,y0);
	k2 = h * (*f)(t0 + 1.0/4.0 * h, y0 + 1.0/4.0 * h * (*fy)(t0,y0) * k1);
	k3 = h * (*f)(t0 + 3.0/8.0 * h, y0 + 3.0/32.0 * h * k1 + 
            9.0/32.0 * h * (*fy)(t0,y0) * k2);
	k4 = h * (*f)( t0 + 12.0/13.0 * h, y0 + 1932.0/2197.0 * h * k1 
            - 7200.0/2197.0 * h * k2 + 7296.0/2197.0 * h * (*fy)(t0,y0) * k3);
	k5 = h * (*f)(t0 + 1.0 * h, y0 + 439.0/216.0 * h * k1 - 8.0 * h * k2 + 
            3680.0/513.0 * h * k3 - 845/4104 * h * (*fy)(t0,y0) * k4 );
	k6 = h * (*f)(t0 - 1.0/2.0 * h, y0 - 8.0/27.0 * h * k1 + 2.0 * h * k2 
            - 3544.0/2565.0 * h * k3 + 1859.0/4104.0 * h * k4 - 
            11.0/40.0 * h * (*fy)(t0,y0) * k5);
		
	*t = t0 + h * i;
	/* Fehlberg of order 4. */
	*y = y0 + 25.0/216.0 * k1 + 0.0 * k2 + 1408.0/2565.0 * k3 + 
        2197.0/4104.0 * k4 - 1.0/5.0 * k5 + 0.0 * k6; 
	/* Fehlberg of order 5. */
	/* *y = y0 + 16.0/135.0 * k1 + 0.0 * k2 + 6656.0/12825.0 * k3 + 
        28561.0/56430.0 * k4 - 9.0/50.0 * k5 + 2.0/55.0 * k6; //5th order */
}
