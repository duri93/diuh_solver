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

% define sensitivity parameters
sensitivityPar.parameter(1) = "f0.uh";
sensitivityPar.label(1) = "uh = %u m/h";
sensitivityPar.values{1} = [36 72]; % [1 0.8] * 1e-2 * 3600;

sensitivityPar.parameter(3) = "f0.k";
sensitivityPar.label(3) = "k = %.2f h^{-1}";
sensitivityPar.values{3} = 1 ./ [1 0.5] / 24;

sensitivityPar.parameter(2) = "alpha.beta";
sensitivityPar.label(2) = "beta = %u m/h";
sensitivityPar.values{2} = [120 240];

%% RUN SIMULATIONS
result = runIUHSensitivityAnalysis(par, sensitivityPar);

%% plots
plot_IUHsensitivity(par, sensitivityPar, result);


%% FUNCTIONS
function result = runIUHSensitivityAnalysis(par, sensitivityPar)
	% allocate results
	Nvalues = cellfun(@numel, sensitivityPar.values);
	Nparameters = numel(Nvalues);
	Ivalues = ones(Nparameters, 1);

	result = cell(Nvalues);

	% try all param value combinations
	for v = 1:prod(Nvalues)
		% update param
		testPar = par;
		for p = 1:Nparameters
			testPar = updateParameter(testPar, sensitivityPar.parameter(p), sensitivityPar.values{p}(Ivalues(p)));
		end
		testPar = makeDependentParameters(testPar, ["domainL.L", "domainL.NL", "domainL.dL"]);
	
		% simulate dynamic
		result{v} = integratedFdL(testPar.domainTau, testPar.domainL, testPar.alpha, testPar.f0);

		% update indices
		Ivalues = updateIndices(Ivalues, Nvalues);
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
function Ivalues = updateIndices(Ivalues, Nvalues)
	Ivalues(1) = Ivalues(1) + 1;

	for i = 1:numel(Ivalues)-1
		if Ivalues(i) > Nvalues(i)
			Ivalues(i) = Ivalues(i) - Nvalues(i);
			Ivalues(i + 1) = Ivalues(i + 1) + 1;
		end
	end
end

%% PLOT FUNCTIONS
function plot_IUHsensitivity(par, sensitivityPar, result)
	% prepare constants
	NL = par.domainL.NL;
	displayNames = arrayfun(@(x) sprintf("L = %.1f km", x/1000), par.domainL.L);
	
	% prepare colors
	fig = figure();
	colors = colormap("copper");
	colors = interp1(linspace(0,1,height(colors)), colors, linspace(0,1, NL));
	close(fig);
	
	% prepare sensitivity parameter values
	Nvalues = cellfun(@numel, sensitivityPar.values);
	Nparameters = numel(Nvalues);
	Ivalues = ones(Nparameters, 1);
	
	ax = gobjects(prod(Nvalues));
	
	for f = 1:prod(Nvalues(1:end-2))
		figure("Position", [(f-1)*700 100 700 600]);
		tiledlayout(Nvalues(end), Nvalues(end-1), "TileSpacing", "tight", "Padding", "tight");
	
		for c = 1:Nvalues(end-1)
			Ivalues(end-1) = c;
	
			for r = 1:Nvalues(end)
				% get combination id and result
				Ivalues(end) = r;
				% id = ind2sub(Nvalues, Ivalues);
				id = f + prod(Nvalues(1:end-2))*(c-1) + prod(Nvalues(1:end-1))*(r-1);
				F = result{id};
	
				% prepare plot
				ax(id) = nexttile();
				hold on;
				box on;
				title(sprintf("%s)", char(c + (r-1)*Nvalues(end-1) + 96)));
				xlabel("\tau [days]");
				ylabel("f(\tau, L) [d^{-1}]");
				leg = legend("Location", "best", "box", "off");
				
				% plot results
				for i = 1:NL
					plot(F.tau / 24, F.f_L_tau(i, :)*24, '-', 'Color', colors(i, :), 'DisplayName', displayNames(i));
				end
	
				% prepare text properties
				position = [leg.Position(1), leg.Position(2) - 0.275, leg.Position(3), 0.25];
	
				values = NaN(Nparameters, 1);
				for p = 1:Nparameters
					values(p) = sensitivityPar.values{p}(Ivalues(p));
				end
			
				string = strjoin(sensitivityPar.label, "\n");
				string = sprintf(string, values);
				
				annotation('textbox', position, 'String', string, 'FitBoxToText','off', 'EdgeColor', 'none');
	
			end
		end
	
		Ivalues = updateIndices(Ivalues, Nvalues);
	end
	
	linkaxes(ax, 'xy');
	ylim([0 1]);
	xlim([0 15]);
end
