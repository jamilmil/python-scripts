# python-scripts
The 4 Python scripts in this repository were driven from lectures that a I took in the TCM3: Programming in Chemistry Lecture, originally written in Fortran,
at the University of Tuebingen held by Pr.Dr.Reinhold.Fink from the Institute of Physical and Theoretical Chemistry.

The monte-carlo.py uses the monte carlo itegration method to estimate the Integrale I(x) in an interval [a,b] through random procedure. The script generates
N random points in the same interval and evalutes wether it lay below the cure of I(x). It counts the points that fullfil this conditions and divides it by
the total number of random points N, providing a rough estimation of I(x).

The zero_of_functions.py uses 3 numerical  methods which are: Regula-Falsi, Newton's and Interval methods to estimate the value of pi
using a Sin(x) function.

The debye_solid.py uses uses 3 numerical integration methods which are: Simpson, Trapezoidal and Gauss-Legender to evaluet the value of an I(x) in an inetravl
[a,b] by sampling it into N subintervals with a specific width. At the end, the heat capacity of some metals is evaluated from this Integrale.

The kinetic.py uses numerical methods for solving 3 differential equations that describe a kinetic reaction (A <=> B <=> C). The objectif is to access the
concentration of all compounds at any time using 3 different methods: Euler, RK2 and RK4 methods.
