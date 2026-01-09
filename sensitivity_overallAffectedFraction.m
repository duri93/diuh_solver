clc; 
clear; 
close all;
addpath("modules/parameters", "modules/model", "modules/plotters");

%% PARAMETERS
% define parameters
par = makeDefaultParameters("pulse", "Valfredda");
par.domainL.dL = 1000;
par.domainTau.tauMax = 180*24;
% note: potentially update parameters here if you want to change the default configuration
par = makeDependentParameters(par);

% define sensitivity parameters
sensitivityPar.parameter = "alpha.beta";
sensitivityPar.values = 10:10:500;

%% RUN SIMULATIONS
F_L_tau = integratedFdL(par.domainTau, par.domainL, par.alpha, par.f0);
[oaf, fmin_tau_beta]  = runOAFSensitivityAnalysis(par, sensitivityPar);

%% plots
plotOAF(sensitivityPar, oaf, F_L_tau);
% plotfLmin(sensitivityPar, fmin_tau_beta, F_L_tau)

%% FUNCTIONS
function [oaf, fmin_tau_beta] = runOAFSensitivityAnalysis(par, sensitivityPar)
	% allocate results
	Nvalues = numel(sensitivityPar.values);
	oaf = NaN(par.domainL.NL, Nvalues);
	fmin_tau_beta = NaN(par.domainTau.Ntau, Nvalues);

	parameter = sensitivityPar.parameter;
	values = sensitivityPar.values;

	% try all param values
	parfor v = 1:Nvalues
		% update param
		testPar = updateParameter(par, parameter, values(v));
		testPar = makeDependentParameters(testPar);
	
		% simulate dynamic
		F_L_tau = integratedFdL(testPar.domainTau, testPar.domainL, testPar.alpha, testPar.f0);

		fmin_tau_beta(:, v) = F_L_tau.f_L_tau(1, :);
		oaf(:, v) = checkOverallAffectedFraction(F_L_tau, testPar.alpha.fun);
	end
end
function par = updateParameter(par, parName, parValue)
	parName = strsplit(parName, ".");

	switch numel(parName)
		case 1
			par.(parName) = parValue;
		case 2
			par.(parName(1)).(parName(2)) = parValue;
		case 3
			par.(parName(1)).(parName(2)).(parName(3)) = parValue;
		case 4
			par.(parName(1)).(parName(2)).(parName(3)).(parName(4)) = parValue;
	end
end

%% PLOTS 
function plotOAF(sensitivityPar, oaf, F_L_tau)
	NL = F_L_tau.NL;
	
	fig = figure();
	colors = colormap("copper");
	colors = interp1(linspace(0,1,height(colors)), colors, linspace(0,1, NL));
	close(fig);
	
	displayNames = arrayfun(@(x) sprintf("L = %.1f km", x/1000), F_L_tau.L);
	

	figure();
	hold on;
	for v = 1:NL
		plot(sensitivityPar.values, oaf(v, :), '-', 'Color', colors(v, :), 'DisplayName', displayNames(v));
	end
	xlabel("Beta");
	ylabel(["Overal affected fraction," "OAF(L) = \int_0^\infty f \alpha d\tau"]);
	legend("Location", "best");
	ylim([0 1]);
end
