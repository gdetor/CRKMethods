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
