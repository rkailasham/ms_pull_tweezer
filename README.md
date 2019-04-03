# ms_pull_tweezer
FORTRAN77 code for performing Brownian dynamics simulations on a dumbbell with fluctuating internal friction and hydrodynamic interactions

Developed at the Molecular Rheology Group (MRG, Supervisor: Ravi Jagadeeshan), Monash University, Australia.
Please read the LICENSE file carefully.


1) ms_tweezer.f is the driver routine that performs Brownian dynamics simulations

2) inp.dat contains the necessary input values for the simulation, in the following order.
The parameters are to be entered in the same line, from left to right. For the sake of clarity,
the descriptions of the parameters are provided in a top-down fashion here.

[Excluded volume interaction parameter (set=0)]
[distance cut-off for Narrow Gaussian excluded volume potential (set any non-zero value)]
[Finite extensibility parameter, b]
[Length scale, l_H (nm)]
[First trap strength, in units of dumbbell spring constant (c_1)]
[x-coordinate of first trap position]
[y-coordinate of first trap position]
[z-coordinate of first trap position]
[Second trap strength, in units of dumbbell spring constant (c_2)]
[x-coordinate of initial position of second trap]
[y-coordinate of initial position of second trap]
[z-coordinate of initial position of second trap]
[x-coordinate of final position of second trap]
[y-coordinate of final position of second trap]
[z-coordinate of final position of second trap]
[Type of hydrodynamic interaction tensor to be used: 2 for Rotne-Prager-Yamakawa, 3 for no HI]
[Default parameter, pick initial configurations from Gaussian distribution, set=1]
[Flow parameter, use (1), (2), or (3). No flow condition implemented in code itself]

3) param.dat contains additional input parameters, in the following order:

[Temperature (K)] [Viscosity of solvent (Pa.s)] [Bead radius (nm)] [Internal friction coefficient (kg/s)]
