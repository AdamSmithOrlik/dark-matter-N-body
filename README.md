# dark-matter-N-body

Building a dark matter N-body simulation

UPDATES:

26/06/24

- Rewrote the sampling code to correctly sample velocities from a Maxwell Boltzmann distribution.
- Added capital constants to control key aspects of the initial conditions.
- Added dimensionful constant, giving time steps a concrete meaning -> - Calculate what this relationship is.
- Added a new softening length according to pg 125 Binney & Termaine -> Not sure if it should be positive or negative (POSITIVE).
- Made a rotational velocity test.

28/06/24

- Made a sun earth orbital system to test
- Fixed velocities by relating them to the radial distance by adding circular velocity and approximating
- Debugged dimensionful quantities by relating G to rhos and rs
- Added loading bar for the simulation
-

TO DO:

1. Figure out how to model orbits as a check that the model is working DONE.
2. Fix the velocities to represent their distance from the center DONE.
3. Rewrite the c++ code with the changes
4. Think about parallelization or schemes to speed up the simulations.

ULTIMATE GOALS:
Full 3D N-body simulation generated in c++ with analysis code in python.

IDEAS:
Add more complicated scattering probability (cross-section).
Trace collisions for each particle in a dataframe.
Figure out how to add asymmetric velocity dispersions, or, nonspherical halo shapes to start.
