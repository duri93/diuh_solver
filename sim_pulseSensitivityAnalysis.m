clc; 
clear; 
close all;
addpath("modules/parameters", "modules/model", "modules/plotters");

%% PARAMETERS
% define parameters
par = makeDefaultParameters("pulse", "Valfredda"); % HubbardBrook25, PovertyCreek, Valfredda, DevonA
par.domainT.D = 10*24;
% note: potentially update parameters here if you want to change the default configuration
par = makeDependentParameters(par);

% define sensitivity parameters
sensitivityPar.parameter = "rain.volume";
sensitivityPar.values = [5 10 20 40]; % [mm]

%% RUN SIMULATIONS

% allocate results
results = table('Size', [numel(sensitivityPar.values), 6], ...
				'VariableTypes', ["string", "double", "cell", "cell", "doublenan", "doublenan"], ...
				'VariableNames', ["parameter", "value", "dynamic", "static", "peakQ", "peakTime"]);
results.parameter = repmat(sensitivityPar.parameter, height(results), 1);
results.value = sensitivityPar.values';

% solve dF/dL
F_L_tau = integratedFdL(par.domainTau, par.domainL, par.alpha, par.f0);
OAF = checkOverallAffectedFraction(F_L_tau, par.alpha.fun);

% run all sims
results = runSensitivityAnalysis(results, F_L_tau, par);

%% plots
% L(t), Q(t), Q/Qmax(t)
plotComparisonLQ(results, par, "I = %.0f mm/h", [0 par.domainT.D]/24);

% timeseries comparison static-dynamic
plotComparisonDynamicStatic(results, par, [0 par.domainT.D]/24, 0.5);

%% FUNCTIONS
function results = runSensitivityAnalysis(results, F_L_tau, par)
	% extract table vars
	parameter = results.parameter;
	value = results.value;

	% allocate
	dynamic = cell(height(results), 1);
	static = cell(height(results), 1);
	peakQ = NaN(height(results), 1);
	peakTime = NaN(height(results), 1);

	% to the parallel tests
	for r = 1:height(results)
		% update param
		parTest = updateParameters(par, parameter(r), value(r));
	
		% simulate dynamic
		result_t = generateRainfallPulse(parTest.domainT, parTest.rain);
		result_t = simulateDynamic(result_t, parTest.alpha.fun, parTest.f0.uh, parTest.LQ, F_L_tau, parTest.rain.baseflow);
		checkMassBalance(result_t, F_L_tau, parTest);

		% simulate static
		static_result_t = simulateStatic(result_t, F_L_tau, parTest.rain.baseflow, max(result_t.L));

		% get statistics
		[peakQ(r), pT] = max(result_t.Q);
		peakTime(r) = result_t.t(pT);
	
		% save
		dynamic{r} = result_t;
		static{r} = static_result_t;
	end

	% save back
	results.dynamic = dynamic;
	results.static = static;
	results.peakQ = peakQ;
	results.peakTime = peakTime;
end
function par = updateParameters(par, parName, parValue)
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
function plotComparisonLQ(results, par, legendFormat, xlimits)
	% prepare values
	fig = figure("Position", [200 100 450 800]);
	colors = colormap("copper");
	colors = interp1(linspace(0,1,height(colors)), colors, linspace(0,1,height(results)));
	close(fig);
	
	displayNames = arrayfun(@(x) sprintf(legendFormat, x), results.value);

	% prepare figure
	figure("Position", [300 75 450 800]);
	tiledlayout(3, 1, "TileSpacing", "tight", "Padding", "tight");
	
	% panel a) L(t)
	ax(1) = nexttile();
	box("on");
	hold("on");
	plot(xlimits, par.domainL.Lmax*[1 1] / 1000, '--', 'Color', 0.4*[1 1 1]);
	plot(xlimits, par.domainL.Lmin*[1 1] / 1000, '--', 'Color', 0.4*[1 1 1]);
	for r = height(results):-1:1
		plot(results.dynamic{r}.t / 24, results.dynamic{r}.L / 1000, ...
			'LineWidth', 1.5, "Color", colors(r, :), "DisplayName", displayNames(r));
	end
	xlabel("t [d]");
	ylabel("L [km]");
	title("a)");
	
	% panel b) Q(t)
	ax(2) = nexttile();
	box("on");
	hold("on");
	for r = height(results):-1:1
		plot(results.dynamic{r}.t / 24, results.dynamic{r}.Q*24, ...
			'LineWidth', 1.5, "Color", colors(r, :), "DisplayName", displayNames(r));
	end
	xlabel("t [d]");
	ylabel("Q [mm/d]");
	title("b)");
	legend("location", "northeast", "box", "off");
	
	% panel c) Q(t) / Qmax
	ax(3) = nexttile();
	box("on");
	hold("on");
	for r = height(results):-1:1
		plot(results.dynamic{r}.t / 24, results.dynamic{r}.Q / results.peakQ(r), ...
			'LineWidth', 1.5, "Color", colors(r, :), "DisplayName", displayNames(r));
	end
	xlabel("t [d]");
	ylabel("Q/Q_{max}");
	title("c)");
	
	% end stuff
	linkaxes(ax(1:3), 'x');
	xlim(ax(1), xlimits);
	set(ax, 'FontSize', 11);
end
function plotComparisonDynamicStatic(results, par, xlimits, peakZoomDuration)
	colorDynamic = [0.65 0.80 0.95]; % azzurro pastello
    colorStatic = [0.98 0.75 0.70]; % pesca pastello

	% prepare figure
	figure("Position", [800 100 600 800]);
	tiledlayout(height(results), 3, "TileSpacing", "tight", "Padding", "tight");
	ax = gobjects(height(results), 2);

	% plot each test
	for r = 1:height(results)
		% whole sim
		ax(r, 1) = nexttile([1 2]);
		box("on");
		hold("on");
		plot(results.dynamic{r}.t / 24, results.dynamic{r}.Q * 24, ...
			"Color", colorDynamic, "LineWidth", 2, "DisplayName", "Dynamic network");
		plot(results.dynamic{r}.t / 24, results.dynamic{r}.Qs * 24, '-', ...
			"Color", [0 0 1], "LineWidth", 0.5, "DisplayName", "Static Q for dynamic network");
		plot(results.static{r}.t / 24, results.static{r}.Q * 24, '-.', ...
			"Color", colorStatic, "LineWidth", 1, "DisplayName", "Static network");
		title(sprintf("%s)", char(r+96)));
		xlabel("t [d]");
		ylabel("Q [mm/d]");
		xlim(xlimits);
		legend("Location", "best");

		% zoom
		ax(r, 2) = nexttile();
		box("on");
		hold("on");
		plot(results.dynamic{r}.t / 24, results.dynamic{r}.Q * 24, ...
			"Color", colorDynamic, "LineWidth", 2, "DisplayName", "Dynamic network");
		plot(results.static{r}.t / 24, results.static{r}.Q * 24, '-.', ...
			"Color", colorStatic, "LineWidth", 1, "DisplayName", "Static network");
		xlabel("t [d]");
		ylabel("Q [mm/d]");

		xlim(par.rain.start/24  + peakZoomDuration*[-0.25 0.75])
		% xlim(results.peakTime(r)/24 + peakZoomDuration*[-0.25 0.75]);

		linkaxes(ax(r, :), 'y');
	end
end
