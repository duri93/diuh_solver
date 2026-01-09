clc; 
clear; 
% close all;
addpath("modules/parameters", "modules/model", "modules/plotters");

%% define parameters
par = makeDefaultParameters("poisson", "Valfredda");
% note: potentially update parameters here if you want to change the default configuration
par.domainT.dt = 1; % h
par = makeDependentParameters(par);

%% simulation

% solve dF/dL
F_L_tau = integratedFdL(par.domainTau, par.domainL, par.alpha, par.f0);
checkOverallAffectedFraction(F_L_tau, par.alpha.fun);

% generate rainfall (or load it if it exists already)
if exist("results\sim_dynamicVsStatic.mat", "file")
	load("results\sim_dynamicVsStatic.mat", "result_t");
else
	result_t = par.rain.type(par.domainT, par.rain);
end

% dynamic simulation
[result_t, result_t_tau] = simulateDynamic(result_t, par.alpha.fun, par.f0.uh, par.LQ, F_L_tau, par.rain.baseflow);
massBalance = checkMassBalance(result_t, F_L_tau, par);

% static simulation
result_t_Lmax = simulateStatic(result_t, F_L_tau, par.rain.baseflow, max(result_t.L));
result_t_Lmed = simulateStatic(result_t, F_L_tau, par.rain.baseflow, mean(result_t.L));
result_t_Lcal = calibrateLStatic(par.domainL, result_t, F_L_tau, par.rain.baseflow);

% save: since it's stochastic, save this exact simulation if you need to update the figure while keeping same rainfall
save("results\sim_dynamicVsStatic.mat");

%% plots
plot_dynamicVsStaticTimeseries(result_t, {result_t_Lmed}, ["mean(L)"], par.domainL);
% plot_dynamicVsStaticTimeseries(result_t, {result_t_Lmax, result_t_Lmed, result_t_Lcal}, ["max(L)", "mean(L)", "calibrated"], par.domainL);
plot_dynamicVsStaticFDC(result_t, {result_t_Lmax}, "Static network, max(L)");

%% FUNCTIONS
function result_t_Lcal = calibrateLStatic(domainL, result_t, F_L_tau, baseflow)
	
	% try many L
	NL = 100;
	L = linspace(domainL.Lmin, domainL.Lmax, NL);
	err = NaN(NL, 1);

	for i = 1:NL
		res_t   = simulateStatic(result_t, F_L_tau, baseflow, L(i));
		err(i) = sum( (res_t.Q - result_t.Q).^2, "omitmissing");
	end

	% show error
	% figure();
	% plot(L, err, '.-');

	% select best one
	[~, i] = min(err);
	result_t_Lcal = simulateStatic(result_t, F_L_tau, baseflow, L(i));
end

function plot_dynamicVsStaticTimeseries(result_t, static_res, static_labels, domainL)
	dt = result_t.t(2) - result_t.t(1);

	% colors and style
	colors = colororder();
	linewidth = 1;
	facealpha = 0.1;
	fontsize = 10;

	% figure layout
	figure("Position", [60 50 800 900]);
	tiledlayout('vertical', 'TileSpacing','tight', 'Padding','tight');

	% === first panel: rainfall
	ax(1) = nexttile();
	hold on;
	box on;
	title("a)");
	ylabel('Rainfall h(t) [mm/d]');

	% bar(result_t.t/24, result_t.hDaily/dt*24, 1, "FaceColor", 0.9*[1, 1, 1], "EdgeColor", "none", "DisplayName", "Daily average");
	bar(result_t.t/24, result_t.h/dt*24, 1, "FaceColor", colors(1,:), "EdgeColor", "none", "DisplayName", "dt-scale rainfall");
	% legend("Location", "best", "Box", "off");

	% === second panel: length L(t)
	ax(2) = nexttile();
	hold on;
	box on;
	title("b)");
	ylabel('Length L [m]');
	ylim([domainL.Lmin domainL.Lmax] / 1000);
	legend("Location", "best", "Orientation", "horizontal", "Box", "off");

	area(result_t.t/24, result_t.L/1000, 0, "FaceColor", colors(1,:), "FaceAlpha", facealpha, "EdgeColor", colors(1, :), "LineWidth", linewidth, "DisplayName", "L(t)");
	for r = 1:numel(static_res)
		res_t = static_res{r};
		plot(res_t.t/24, res_t.L/1000, '-.', "Color", colors(r+1,:), "LineWidth", linewidth, "DisplayName", static_labels(r));
	end
	
	% === third panel: flowrate Q(t) ===
	ax(3) = nexttile();
	hold on;
	box on;
	title("c)");
	xlabel('Time [days]');
	ylabel('Flowrate Q [mm/d]');
	legend("Location", "best", "Orientation", "horizontal", "Box", "off");

	% plot(result_t.t([1 end])/24, baseflow*[1 1]*24, '--', 'color', blue, 'DisplayName', 'Constant baseflow');
	% plot(result_t.t/24, result_t.Qs*24, '-', 'color', blue, 'DisplayName', 'Static component of Q with dinamic network');
	area(result_t.t/24, result_t.Q*24, 0, "FaceColor", colors(1,:), "FaceAlpha", facealpha, "EdgeColor", colors(1, :), "LineWidth", linewidth, "DisplayName", "Dynamic network");
	for r = 1:numel(static_res)
		res_t = static_res{r};
		plot(res_t.t/24, res_t.Q*24, '-.', "Color", colors(r+1,:), "LineWidth", linewidth, "DisplayName", strcat("Static network, ", static_labels(r)));
	end
	
	% --- finishing touches
	linkaxes(ax,"x");
	xlim([30 90]);

	for a = 1:3
		ax(a).FontSize = fontsize;
	end
end

function plot_dynamicVsStaticFDC(result_t, static_result_t, static_labels)
	% constant params
	Nclasses = 200;
	Qrange = [0 15];
	Nstatic = numel(static_result_t);
	
	% prepare pdfs
	dynamic_pdf = pdfMaker(24*result_t.Q, Nclasses, Qrange);

	static_pdf = cell(Nstatic, 1);
	for n = 1:Nstatic
		static_pdf{n} = pdfMaker(24*static_result_t{n}.Q, Nclasses, Qrange);
	end

	% prepare figure
	figure("Position", [100 100 400 670]);
	tiledlayout(2, 1);

	% panel a) flow pdf
	nexttile();
	hold on;
	box on;
	xlabel("Q [mm/d]");
	ylabel("P(Q) [d/mm]");
	legend("location", "best", "box", "off");
	xlim(Qrange);
	ylim([0 0.3]);
	title("a)");
	plot(dynamic_pdf.x, dynamic_pdf.pdf, '-', 'DisplayName', "Dynamic network");
	for n = 1:Nstatic
		plot(static_pdf{n}.x, static_pdf{n}.pdf, '-.', 'DisplayName', static_labels(n));
	end

	% panel b) flow duration curve (semilog scale)
	x = linspace(1, 0, height(result_t));

	nexttile();
	hold on;
	box on;
	xlabel("Duration [-]");
	ylabel("Q [mm/d]");
	legend("location", "best", "box", "off");
	xscale('log');
	title("b)");
	% ylim([0 20]);
	% xlim([0 1]);
	plot(x, sort(24*result_t.Q), '-', 'DisplayName', "Dynamic network");
	for n = 1:Nstatic
		plot(x, sort(24*static_result_t{n}.Q), '-.', 'DisplayName', static_labels(n));
	end
end