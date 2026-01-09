function plot_posteriorLQ(result, LQfun)
	Q = linspace(min(result.Q), max(result.Q));
	L = LQfun(Q);

	figure();
	hold on;
	plot(result.Q*24, result.L/1000, '.', 'DisplayName', 'Posterior L(Q)');
	plot(Q*24, L/1000, '-k', 'DisplayName', 'Prior L(Q)');
	xlabel('Q [mm/d]');
	ylabel('L [km]');
	legend('location', 'best');
end