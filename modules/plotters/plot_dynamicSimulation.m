function plot_dynamicSimulation(result, resultF, Ldomain)

	% extract vectors
	t = result.t / 24; % days
	dt = t(2) - t(1); % days
	% trange = [min(t) max(t)];

	% figure
    figure();
	tiledlayout(3, 2);

	% rainfall
	ax(1) = nexttile(1);
	hold on;
	grid on;
	bar(t-dt/2, result.h/dt, 1);
	xlabel("t [d]");
	ylabel("Rainfall [mm/d]");

	% length
	ax(3) = nexttile(3);
	hold on;
	grid on;
	plot(t, result.L/1000, '-', 'DisplayName', "L(t)");
	plot([min(t) max(t)], Ldomain.Lmin*[1 1]/1000, ':k', 'DisplayName', 'Lmin');
	plot([min(t) max(t)], Ldomain.Lmax*[1 1]/1000, ':k', 'DisplayName', 'Lmax');
	bar(t-dt/2, double(result.c == 2)*max(result.L)/1000, 1, ...
		"EdgeColor", "none", "FaceColor", [0 0 0], "FaceAlpha", 0.15, "DisplayName", "Contraction correction");
	bar(t-dt/2, double(result.c == 3)*max(result.L)/1000, ...
		1, "EdgeColor", "none", "FaceColor", [0 0 0], "FaceAlpha", 0.3, "DisplayName", "Countercorrection");
	xlabel("t [d]");
	ylabel("L [km]");
    %ylim([Ldomain.Lmin/1000, Ldomain.Lmax/1000])
	legend('Location', 'best');

	% flowrate
	ax(5) = nexttile(5);
	hold on;
	grid on;
	plot(t, result.Q*24, '-');
	ylabel("Q [mm/d]");
	xlabel("t [d]");

	% F_t_tau
	ax(2) = nexttile(2);
	hold on;
	% waterfall(resultF.tau/24, resultF.t/24, resultF.F_t_tau);
	imagesc(resultF.t/24, resultF.tau/24, resultF.F_t_tau');
    title('F(\tau, t)');
    ylabel('\tau [d]');
    xlabel('t [d]');
	clim([0 1]);
	ylim([0 20]);
	colorbar();

	% dFdtau
	ax(4) = nexttile(4);
	hold on;
	% waterfall(resultF.tau/24, resultF.t/24, resultF.dFdtau_t_tau);
	imagesc(resultF.t/24, resultF.tau/24, resultF.dFdtau_t_tau');
    title('dF/d\tau(\tau, t)');
    ylabel('\tau [d]');
    xlabel('t [d]');
	ylim([0 20]);
	clim([0 0.006]);
	colorbar();

	% dFdt
	ax(6) = nexttile(6);
	hold on;
	% waterfall(resultF.tau/24, resultF.t/24, resultF.dFdt_t_tau);
	imagesc(resultF.t/24, resultF.tau/24, resultF.dFdt_t_tau');
    title('dF/dt(\tau, t)'); 
    xlabel('\tau [d]');
    ylabel('t [d]');
	ylim([0 20]);
	clim([0 0.006]);
	colorbar();

	% finishing touches
	linkaxes(ax, 'x');
end