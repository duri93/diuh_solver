clc; 
clear; 
close all;
addpath("modules/parameters", "modules/model", "modules/plotters");

%% PARAMETERS
% define parameters
par = makeDefaultParameters("pulse", "Valfredda");
par.domainL.dL = 5000;
par.domainTau.tauMax = 20*24;
% note: potentially update parameters here if you want to change the default configuration
par = makeDependentParameters(par);

%% RUN SIMULATIONS
F_L_tau = integratedFdL(par.domainTau, par.domainL, par.alpha, par.f0);
F_L_tau.alpha = par.alpha.fun(repmat(F_L_tau.L, 1, F_L_tau.Ntau), repmat(F_L_tau.tau', F_L_tau.NL, 1));

%% PLOT
plot_alpha_iuh(F_L_tau);

%% FUNCTIONS
function plot_alpha_iuh(F_L_tau)
	% prepare values
	NL = F_L_tau.NL;
	
	% prepare colors
	fig = figure("Position", [200 100 450 800]);
	colors = colormap("copper");
	colors = interp1(linspace(0,1,height(colors)), colors, linspace(0,1,NL));
	close(fig);
	
	% prepare figure
	fig = figure("Position", [200 100 450 800]);
	tiledlayout(3, 1, "TileSpacing", "tight", "Padding", "tight");
	
	% panel a) maps
	nexttile();
	title("a)");
	
	% panel b) alpha
	nexttile();
	hold on;
	box on;
	title("b)");
	xlabel("\tau [days]");
	ylabel("\alpha(\tau, L)");
	legend("Location", "best", "box", "off");
	for i = 1:NL
		str = sprintf("L = %.1f km", F_L_tau.L(i)/1000);
		plot(F_L_tau.tau/24, F_L_tau.alpha(i, :), '-', 'Color', colors(i, :), "DisplayName", str);
	end
	
	% panel c) f
	nexttile();
	hold on;
	box on;
	title("c)");
	xlabel("\tau [days]");
	ylabel("f(\tau, L) [days^{-1}]");
	legend("Location", "best", "box", "off");
	for i = 1:NL
		str = sprintf("L = %.1f km", F_L_tau.L(i)/1000);
		plot(F_L_tau.tau/24, F_L_tau.f_L_tau(i, :), '-', 'Color', colors(i, :), "DisplayName", str);
	end
end