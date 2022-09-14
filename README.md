## CRKMethods

CRKMEthods provides some of the most used Runge-Kutta methods in the C
Programming Language. The algorithms that have been implemented as one-step
are:
  - Forward Euler
  - RK4
  - RK4 refined
  - RK - Fehlberg

More information about those solvers can be found on
[Runge-Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)'s Wikipedia
web-page and [Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method)'s
wikipedia webpage.

### Numerical precision

The user can define the numerical precision from the *solver.h* file. They 
can define the macro `REAL` to either `float` or `double`.


### Compile an example
 
The given example can be compiled from the parent directory by executing the
following command:
```unix
$ gcc -I./include src/solvers.c src/main.c
```
