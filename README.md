
# DIUH MODEL
River networks systematically expand and contract as a function of climatic input and catchment wetness. As a consequence, 
systematic variations in water drainage path length and hydrologic response time arises. This code implements a generaized 
unit IUH theory that explicitly accounts for time-varying networks, implementing a dynamic formulation of the convolution 
relationship linking effective rainfall and discharge.

The main features are:
- a specified exponential IUH for a given network L
- the impact coefficient specifies how the IUH evolves as a function of L
- rainfall can be generated synthetically as a poisson process
- a length-discharge relation is used to define the temporal evolution of L
- the generalized convolution is used to transform effective rainfall into discharge

## FOLDER STRUCTURE
- `./`						(this folder) contains the scripts to run to build the figures reported in the paper
- `./modules/model`			contains all the functions related to model behaviour
- `./modules/parameters`	contains all the functions related to the definition of the parameters
- `./modults/plotters`		contains some functions with additional plot possibilities.
- `./results`				used to save & load stocastically generated simulations

## PARAMETERS
All parameters and default values are defined in `/modules/parameters/definitions.m`. See comments reported in that file
for more details. Some parameters are automatically calculated based on others (e.g., the number of timesteps is 
calculated as a function of the simulation time step and total duration. Some parameters can be selected from a list of 
available catchments. See definitions.m file for details.

## SCRIPTS
Note: some scripts require the parallel toolbox in order to work. This happens when many simulations need to be performed.

- `example_alpha_iuh.m`
	Uses default parameter values to produce an example of how the impact coefficient and the IUH varies with network length.

- `sensitivity_IUH.m`
	Tests a number of values of uh, k and beta to produce the sensitivity analyisis of the IUH = f(tau, L)

- `sensitivity_overallAffectedFraction.m`
	Shows the relative fraction of catchment affected by network dynamics as a function of the parameters.

- `sim_dynamicVsStatic.m`
	Produces a long-term stochastic simulation comparing flowrate obtained with a static and a dynamic river network.

- `sim_pulseSensitivityAnalysis.m`
	Runs a number of simulations with a single rain pulse of increasing intensity, and compares the results.

## TODO
- [ ] fix instability for small dt
- [ ] Counter-correction for G during re-expansion

### Lower priority
- [ ] link to paper when available
- [ ] explicit integration resolution for dFdL (dL and dtau, can be important if uh is small)
- [ ] integration error check for dFdL
- [ ] additional default catchments in makeDefaultParameters
