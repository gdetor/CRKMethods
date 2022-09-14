/* Header file for CRKMethods solvers.
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
#ifndef CRKSOLVERS_H
#define CRKSOLVERS_H
#include <stdlib.h>

#define REAL float


REAL f(REAL , REAL);
REAL fy(REAL , REAL);

void eulerStep(REAL (*)(REAL,REAL), REAL, REAL, REAL, int, REAL *, REAL *);
    
void rk4Step(REAL (*)(REAL,REAL), REAL, REAL, REAL, int, REAL *, REAL *);

void rk4RefinedStep(REAL (*)(REAL,REAL), REAL (*)(REAL,REAL),
                    REAL, REAL, REAL, int, REAL *, REAL *);

void rkFehlbergStep(REAL (*)(REAL,REAL), REAL (*)(REAL,REAL), REAL, REAL,
                    REAL, int, REAL *, REAL *);

#endif /* CRKSOLVERS_H */
