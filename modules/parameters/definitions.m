% THIS FILE PROVIDES DEFINITONS AND DEFAULT VALUES / OPTIONS FOR ALL PARAMETERS
function [par, catchments, derived] = definitions(rainType)

% For catchment-defined parameters, here we provide default values. Values will be overwritten if a different catchment
%	name is specified. See the second part of this file for a table of catchments and parameter values.

% For derived parameters, no need to specify a value here. The proper value is calculated based on the other paramters
%	using the function "makeDependentParamaters" and the definitions in the last section of this file.










%% ================================================================= SECTION 1: PARAMETERS DEFINITION AND DEFAULT VALUES

% ==== TIME DOMAINS
% chronological time
par.domainT.D = 90*24;			% [h]	 total duration of simulation
par.domainT.dt = 1/8;			% [h]	 fixed time step
par.domainT.Nt = NaN;			% [-]	 (derived parameter) number of time steps

% tau time (time after rainfall falling)
par.domainTau.tauMax = 365*24;	% [h]	 make sure it's long enough such that IUH is always complete
par.domainTau.dtau = NaN;		% [h]	 (derived parameter) resolution in tau = resolution in t to simplify convolution 
par.domainTau.Ntau = NaN;		% [-]	 (derived parameter) number of steps in tau

% ==== LENGTH DOMAIN
par.domainL.dL = 50;			% [m]	 resolution of length domain for f(tau, L) (dF/dL is solved at 1m resolution)
par.domainL.A = 5;				% [km2]  (catchment parameter) catchment area 
par.domainL.phiC = 0.75;		% [-]	 (catchment parameter) non-perennial fraction 
par.domainL.Dg = 3.68;			% [km-1] (catchment parameter) geomorphic drainage density 
par.domainL.Lmax = NaN;			% [m]	 (derived parameter) maximum network length 
par.domainL.Lmin = NaN;			% [m]	 (derived parameter) minimum network length 
par.domainL.Lrange = NaN;		% [m]	 (derived parameter) range of variability of network length (max-min)

% ==== IUH AS FUNCTION OF LENGTH, f(tau, L)
% alpha function (impact coefficient)
par.alpha.beta = 120;			% [m/h]	 parameter of the hyperbolic function for alpha
par.alpha.fun = NaN;			%		 (derived parameter) hyperbolic functional form of alpha 

% boundary condition for f
par.f0.k = 1/24;				% [1/h]  inverse of mean response time for exponential iuh as boundary cond of f
par.f0.uh = 1e-02 * 3600;		% [m/h]  water celerity in hyporheic zone of dry network;
par.f0.L0 = NaN;				% [m]    (derived parameter) where the initial condition is applied 
par.f0.fun = NaN;				% [1/h]  (derived parameter) boundary condition (exponential iuh at L = L0)

% ==== LENGTH DISCHARGE RELATION
% see SI of this paper for empirical values: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024WR038732
par.LQ.k = 0.767;				% [-]	  (catchment parameter) shape parameter of gamma function
par.LQ.theta = 4.45/24;			% [mm/h]  (catchment parameter) scale parameter of the distribution
par.LQ.fun = NaN;				%		  (derived parameter) functional form of L(Q)

% ==== rainfall parameters
par.rain.baseflow = NaN;		% [mm/h]  (derived parameter) constant baseflow to be added to flowrate 

if ~exist("rainType", "var")
	rainType = "missing";
end
switch lower(rainType)
	case "pulse" 
		% define rainfall (single constant event)
		par.rain.type = @generateRainfallPulse;
		par.rain.start = 24;				 % [h]  starting hour of rainfall
		par.rain.duration = 1;				 % [h]  duration of rainfall
		par.rain.volume = 40;				 % [mm] total rain volume (100/24*rain.duration = 100 mm/d)
		par.rain.shape = "rectangular";		 %		"rectangular", "triangular"
	case "poisson"
		% define rainfall (poisson)
		par.rain.type = @generateRainfallPoisson;
		par.rain.alpha = 15;				 % [mm/d] avg daily rain depth, exp distributed
		par.rain.lambda = 0.3;				 % [1/d]  avg daily rain frequency (i.e., probability of rainfall for any given day)
		par.rain.downscale = "conservative"; %        downscaling type: "linear", "triangular", "conservative"
		par.rain.downscaleLimit = 80/24;	 % [mm/h] rain intensity limit after downscaling
		par.rain.downscaleBeta = 1;			 %		 beta distrib param to extract downscaling factors.
		%											 1 = uniform distribution; >1 = peaked around middle
	otherwise
		error("rainType must be either 'pulse' or 'poisson'");
end






%% ============================================================================= SECTION 2: CATCHMENT-DEFINED PARAMETERS
% These paramters define catchment size, length variability and L(Q) shape for a number of catchments.
% See SI of this paper for empirical values: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024WR038732

catchments = table('Size', [0 6], 'VariableTypes', ["string", repmat("doublenan", 1, 5)]);
catchments.Properties.VariableNames = ["name", "A", "Dg",	"phiC", "LQ_k", "LQ_theta"];
catchments(end+1, :) = {"HubbardBrook25",	0.25,	8.56,	0.79,	0.767,	4.43/24};
catchments(end+1, :) = {"Valfredda",		5.21,	3.68,	0.75,	0.175,	61.9/24};
catchments(end+1, :) = {"PovertyCreek",		0.35,	5.13,	0.99,	1.43,	19.2/24};
catchments(end+1, :) = {"DevonA",			0.35,	6.91,	0.84,	0.321,	79.3/24};






%% ======================================================================================= SECTION 3: DERIVED PARAMETERS
% This section defines how derived parameters are calculated based on the pre-defined parameters.

derived = {
	"domainL.Lmax",		@(par) ceil(par.domainL.A * par.domainL.Dg * 1000);
	"domainL.Lmin",		@(par) floor((1 - par.domainL.phiC) * par.domainL.Lmax);
	"domainL.Lrange",	@(par) par.domainL.Lmax - par.domainL.Lmin;
	"domainL.NL",		@(par) ceil(par.domainL.Lrange / par.domainL.dL) + 1;
	"domainL.L",		@(par) linspace(par.domainL.Lmin, par.domainL.Lmax, par.domainL.NL);
	"domainL.dL",		@(par) par.domainL.L(2) - par.domainL.L(1);

	"alpha.fun",		@(par) @(L, tau) par.alpha.beta * tau ./ (par.alpha.beta * tau + L);
	
	"f0.L0",			@(par) par.domainL.Lmax;
	"f0.fun",			@(par) @(tau) par.f0.k * exp(-par.f0.k * tau);
	
	"LQ.fun",			@(par) @(Q) par.domainL.Lmax * gammainc(Q / par.LQ.theta, par.LQ.k);
	
	"domainT.Nt",		@(par) ceil(par.domainT.D / par.domainT.dt) + 1;
	"domainT.t",		@(par) linspace(0, par.domainT.D, par.domainT.Nt);
	"domainT.dt",		@(par) par.domainT.t(2) - par.domainT.t(1);
	
	"domainTau.Ntau",	@(par) ceil(par.domainTau.tauMax / par.domainT.dt) + 1;
	"domainTau.tau",	@(par) linspace(0, par.domainTau.tauMax, par.domainTau.Ntau);
	"domainTau.dtau",	@(par) par.domainTau.tau(2) - par.domainTau.tau(1);

	"rain.baseflow",	@(par) fzero(@(q) par.LQ.fun(q) - par.domainL.Lmin, [0 1000]);
};