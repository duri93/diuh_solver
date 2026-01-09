function oaf = checkOverallAffectedFraction(F_L_tau, alphaFun)

	L = repmat(F_L_tau.L, 1, F_L_tau.Ntau);
	tau = repmat(F_L_tau.tau', F_L_tau.NL, 1);
	oaf = trapz(F_L_tau.dtau, F_L_tau.f_L_tau .* alphaFun(L, tau), 2);

	fprintf("Overall affected fraction range: %.2f - %.2f\n", min(oaf), max(oaf));

	if min(oaf) < 0.1
		warning("Warning: overall affected fraction is less than 10%");
	end
end