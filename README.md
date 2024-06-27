# dark-matter-N-body

Building a dark matter N-body simulation

UPDATES:
Rewrote the sampling code to correctly sample velocities from a Maxwell Boltzmann distribution.
Added capital constants to control key aspects of the initial conditions
Added dimensionful constant, giving time steps a concrete meaning -> Calculate what this relationship is
Added a new softening length according to pg 125 Binney & Termaine -> Not sure if it should be positive or negative
Made a rotational velocity test

TO DO:
c++ functions and constants code now needs to be updated to reflect the above changes
Test starting positions and velocities in c++
Finish a working version in python and then translate into c++ code

ULTIMATE GOALS:
Full 3D N-body simulation generated in c++ with analysis code in python

IDEAS:
Add more complicated scattering probability (cross-section)
Trace collisions for each particle in a dataframe
Figure out how to add asymmetric velocity dispersions, or, nonspherical halo shapes to start
