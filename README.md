# HPCFS
High-Pressure Compressible Flow Solver

This is written and developed by Marc Bernades during the PhD (2021-2024). The solver is open source but its distribution is forbidden unless previously accepted, each usage must be cited accordingly.

The solver is MATLAB-based DNS with the main target to serve a flexible tool to develop methods, numerical schemes and design new cases and approaches which includes both low-pressure (ideal-gas) and high-pressure (real-gas) frameworks. It is not designed to run large 3D DNS cases although 3D cases can also be computed. Nevertheless, it is recommended to use HPC-based solvers to run large scales problems. Our research group provides RHEA, an open source C++ DNS solver to run on extra-scale.

To this extent, this DNS code supports numerical methods to simulate transcritical turbulent flows based on Peng-Robinson equation of state and high-pressure coefficients.

The DNS can compute the tyical canonical geometries on a cartasian domain with potential to stretch on either direction. In this regard, several tests are available at this current distribution:
- 1D advective
- 2D TGV
- 2D Mixing layer
- 2D Channel flow
- 3D TGV inviscid and Viscous
- 3D Channel flow

Each of the "Main" matlab files will need to be adjusted based on test objectives and the functions to run will be common for all. In each test file the user needs to set:
- Solver: Ideal or Real
- Substance: defined on substance library function, currently implemented for N2 and CO2 or ideal-gas
- High-pressure coefficients (mu, kappa): highpressure, power, sutherland or constant
- Boundary conditions
- Domain: origin (x,y,z) and size
- Convective scheme
- 
