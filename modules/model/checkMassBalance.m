function balance = checkMassBalance(result_t, F_L_tau, par)
	
	% calc volumes during simulation
	Vin = sum(result_t.h); % [mm]
	Vout = trapz(result_t.t(1:end-1), result_t.Q(1:end-1) - par.rain.baseflow); % [mm]

	% calc volume not yet drained by the end of the simulation
	tmentau = result_t.t(end) - F_L_tau.tau;
	idL = find(F_L_tau.L >= result_t.L(end), 1);
	dt = par.domainT.dt;
	dtau = par.domainTau.dtau;

	J_tau = interp1(result_t.t, result_t.h, tmentau, 'previous', 0);
	Vafter = (1 - F_L_tau.F_L_tau(idL, :)) * J_tau * dtau / dt;
	
	% calc errors
	massError = Vin - Vout - Vafter;
	percentError = massError / Vin * 100; % [%]
	
	% display error
	% fprintf('Mass balance error [mm, %%]: %.4g, %6.2g\n', massError, percentError);

	if percentError > 0.5
		fprintf("WARNING: mass balance error is %.1f mm = %.1f%%, which is bigger than 0.5%%.\n", massError, percentError);
		% warning("Relative error is larger than 0.5%");
	end

	% store values
	balance.Vin = Vin;
	balance.Vout = Vout;
	balance.Vafter = Vafter;
	balance.massError = massError;
	balance.percentError = percentError;
end