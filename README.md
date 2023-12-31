# HPCFS
High-Pressure Compressible Flow Solver

This is written and developed by Marc Bernades during the PhD (2021-2024). The solver is open source but its distribution is forbidden unless previously accepted at request, each usage must be cited accordingly.

1. Introduction
   
The solver is MATLAB-based high-fidelity DNS with the main target to serve a flexible tool to develop methods, numerical schemes and design new cases and approaches which includes both low-pressure (ideal-gas) and high-pressure (real-gas) frameworks. The solver is well-commented, but in addition to its comments there is a "Instructions" guide with the main setup steps and definitions to avoid confusion.
The solver is not designed to run large 3D DNS cases although 3D (and large) cases can also be computed, in other words, it is not parallelized and its computation speed is not optimized and limited to the single-core machine capacity. Nevertheless, it is recommended to use HPC-based solvers to run large scales problems. Our research group provides RHEA, an open source C++ DNS solver to run on extra-scale.

To this extent, this DNS code supports numerical methods to simulate transcritical turbulent flows based on Peng-Robinson equation of state and high-pressure coefficients based on finite difference second order discretization stencil.

The DNS can compute the tyical canonical geometries on a cartasian domain with potential to stretch on either direction. The current tests are:
- 1D High-pressure sweep
- 1D Advective
- 2D TGV
- 3D TGV (viscous and inviscid)
- 2D Mixing layer (also periodic version)
- 2D Vortex advection
- 2D Lid-driven cavity
- 2D Channel flow (and 3D)

Nevertheless, the main purpose of this distribution is to replicate the results obtained on the numerical assessment benchmark (Kinetic-energy- and pressure-equilibrium-preserving schemes for real-gas turbulence in the transcritical regime, Bernades et al 2023, JCP). To this extent the test supplied with extensive documentation to mimic the results are based on the following, although the repositary contains the complete wide range of tests designed:
- 1D Advective
- 2D Mixing layer

Furthermore, the user can implement any further test, providing that the new test dependant (name, initial conditions and any additional feature are correctly modelled).

2. Discretization
   
The solver is discretized on a finite-difference second-order scheme on a cartesian grid (X,Y,Z) with uniform mesh with stretching functionality along each of the direction seperatly.
The discretization has N inner grid points + 2 additional outer grid points which are used to set and adjust the boundary conditions (reference to the instruction guide).
The main objective of the PhD is working on high-pressure turbulence, hence, within this environment the simulations are strongly susceptible to numerical instabilities due to the presence of nonlinear thermodynamic phenomena and large density gradients, which can trigger spurious pressure oscillations that may contaminate the solution and even lead to its divergence.
Consequently, it is highly beneficial that the numerical schemes utilized, in addition to being kinetic-energy preserving (KEP), attain the so-called pressure-equilibrium-preservation (PEP) property. The numerical scheme utilized in this work has been developed specifically to be simultaneously KEP and PEP. The latter property is achieved by solving a pressure evolution equation.
In brief, the transport equations are numerically solved by adopting a standard semi-discretization procedure; viz. they are first discretized in space and then integrated in time.
In particular, spatial operators are treated using second-order central-differencing schemes, and time-advancement is performed by means of a third-order (or four-order) strong-stability preserving (SSP) Runge-Kutta explicit approach with Courant–Friedrichs–Lewy (CFL) stability criterion.
The temporal errors that arise due to the time-integration scheme are assumed to be kept under control by using sufficiently small time steps.
The convective terms are expanded according to the Kennedy-Gruber-Pirozzoli (KGP) splitting, which has been recently assessed for high-pressure supercritical fluids turbulence.
As a result, the method utilized (i) preserves kinetic energy by convection, (ii) is locally conservative for mass and momentum, (iii) preserves pressure equilibrium and (iv) yields stable and robust numerical simulations without adding any numerical diffusion to the solution or stabilization procedures.

3. Thermodynamics
   
The thermodynamic space of solutions for the state variables pressure $P$, temperature $T$, and density $\rho$ of a single substance is described by an equation of state.
One popular choice for systems at high pressures, which is used in this study, is the Peng-Robinson. In addition, the user will have the ideal-gas equation of state available
The high pressures involved in the analyses conducted in this work prevent the use of simple relations for the calculation of the dynamic viscosity $\mu$ and thermal conductivity $\kappa$.
In this regard, standard methods for computing these coefficients for Newtonian fluids are based on the correlation expressions proposed by Chung et al.
These correlation expressions are mainly function of critical temperature $T_c$ and density $\rho_c$, molecular weight $W$, acentric factor $\omega$, association factor $\kappa_a$ and dipole moment $\mathcal{M}$, and the NASA 7-coefficient polynomial (Burcat, 2005); further details can be found in dedicated works, like for example Poling.

4. Validation

The main validation test results are also enclose on the repositary, but they can be summarized as follows. The roboustness of the code has ensure sucessfull results in all the analysis and investigations performed.

- The former order of accuracy of the spatial discretization (second order) and the time intergator (Runge-kutta 4th order - or the order selected by user).
- Real gas Thermodynamics: Peng-Robinson equation of state implementation and high-pressure coefficients based on Chung et al. validated against NIST in different sub-, trans- and supercritical conditions (for example using the 1D High-pressure test to evaluate the thermodynamic model at different input conditions).
- Viscous scheme validated using lid-driven cavity (Ghia et al. 1982)
- 2D Vortex advection based on Pirozzoli JCP 2011 Vortex test
- 3D TGV ideal-gas inviscid to validate the invariants of transported and induced quantities (Coppola AMR 2019)
- Enthalpy split comparison discretizing rhoE or independently the internal and kineatic energy counterparts.

5. Data output

The discretization is computed based on meshgrid MATLAB functions, as a consequence the grid is transposed so that x-direction is horizontal on the physical space and y-direction vertical. Therefore, the first position of the matrix is for y-direction, second for x-direction and third for z-direction, it is important to consider for plotting and post-processing.
Based on the user pre-defined timestamp, different snapshots will be saved. There will be a pair of files, (i) a file containing the instantaneous fields and (ii) a file containing time-data for each time step up to this instant, these will be the first and second order invariants and kinetic energy.
The structure of data output will be the name defined plus the timestamp in csv format for instanteneous snapshot and the time will be distinguished by "_time" in the end of the file name for the invariants file.

6. Postprocess

The repository will also contain postprocess functions to evaluate the snapshots saved from the simulation outputs and additional plot functionalities.
In any case, the results can be postprocessed by the user with any other tool using the output data in comma-separated values format.
