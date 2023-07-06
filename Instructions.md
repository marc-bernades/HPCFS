
HPCFS instructions, setup and procedure

1. Set Problem Parameters

Each of the "Main (+ test ID)" matlab files will need to be adjusted based on test objectives and the functions to run will be common for all. In each test file the user needs to set:
- Solver: Ideal or Real.
- Substance: defined on substance library function, currently implemented for N2 and CO2 or ideal-gas
- Fluid properties: specially for ideal-gas case the R, gamma, reference properties should be listed
- HP model - High-pressure coefficients (mu, kappa): highpressure, power, sutherland or constant+
- Problem parameters:
    > x_0, y_0, z_0: domain origin
    > L_x, L_y, L_z: length domain
    > t_c, t_end: characteristic time and final time (simulation to abort)
    > Test ID (important to set the label correctly so that it can correctly switch cases)
- Computational parameters
    > grid points in each direction. The un-used dimensions will be set to 1 single node i.e., 1D adv num_y = num_z = 1, and 2D case num_z = 1 with periodic boundaries
    > CFL
    > Max number of iterations
    > Snapshot for iterations: it will produce a file at each of the instant positions of vector
    > Runge-kutta order, i.e, RK = 4
- Stretching: under Fluix.A_x,y,z defined according to the Calculate_stretching_factor function. Currently only implemented for channel flow
- Source terms: at each transported variable
- Controller: only applies for channel flow to control de body force to obtain correct target (Re_tau)
- Numerical scheme: attention is needed here as multiple options are avalable
    > bPressureModel: enabling to 1 transports the pressure in-lieu of E.
    > Filter: enable by setting .Gate = 1 and select type either Implicit or Gaussian
    > Convective matrix form: this should be bEnergySplit = 1 and bEnergySplitEnthalpy = 0 as the scheme is linearly defined below.
    > Flux form: enable .gate = 1 only if any of the flux-based schemes is to be assessed (HLLC, HLLC_plus, etc) or Hybrid
    > Scheme definition: This is the definition for linear combination of several split forms. Set .def = 1. This is the main functionality of the flexibility of schemes:
      > Structure for each of the transported variables with the enthalpy split for E, hence linear combination needs to be defined for kinetic energy (k), internal energy (e) and pressure gradient (p)
      > Convection scheme definition weights are as commented on the code following the order of C_D (divergence), C_fi, C_u (advective), C_rho and C_L
- Viscous scheme: this needs to be set as it is by default, i.e. bViscous_5PointsStencil = 0 (hence 3-point stencil based), bkappaVariable = 1 and bViscosityVariable = 1
- Output file name: name_file_out
- Boundary conditions: 3 options available, i.e., Periodic (1), Neumann (2) and Direchlet (3) to be defined for each border (west, east, north, south, back and front)
The boundary conditions are to be defined based on the self-intuitive comments on the code BC for BC_type & [u,v,w,P,T]
    > Periodic: {BC_types{1},0,0,0,0,0}
    > Neumann:  {BC_types{1},0,0,0,0,0}
    > Dirichlet: Set P<0 for impermeable boundary, Set T<=0 for impermeable boundary and define the desired u,v,w,P or T
      Feature for summetry by setting u=v=w=-1 is defined


2. Compressible flow solver (HPCFS)

Based on "compressibe_solver.m" function where the inputs are previously defined by the problem parameters. The structure of this function is the same regardless the type of test
This will save on the workspace the output functions. However, data will be saved for each snapshot as well which can be further postprocessed

3. Embedded setup functions

In case user desires to adjust further flexibility to create new tests or add functionalities, here a brief of test-dependant functions is listed:
- Initialize_Fields: sets the u,v,w,P,T fields depending on the selected test. Hence, user can modify this function to obtain any additional initial conditions.
For analytical solution-known tests, the initial condition is based on "Reference_Solution_X.m" function (i.e., 1D advective test)
- Filter_Matrix: one of the filter functions if implicit F2 or F4 options. Hence, inside these functions the user cna adjust the alpha_f governing parameter, while the type of filter (F2 or F4) for an implicit selection will be defined within the main function "Compressible_solver.m"
