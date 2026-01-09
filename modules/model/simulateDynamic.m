function [result_t, result_t_tau] = simulateDynamic(result_t, alphaFun, uh, LQ, F, baseflow)
	fprintf('Simulating... ');

	% define domains
	time  = result_t.t;
	dt = time(2) - time(1);
	Nt = height(result_t);

	tau = F.tau;
	dtau = tau(2) - tau(1);
	Ntau = numel(tau);

	L = F.L;
	Lmin = L(1);
	dL = L(2) - L(1);
	NL = numel(L);

	% previous timestep status
	Ltm1 = LQ.fun(baseflow);
	
	% extract vectors from table for speed
	h_t = result_t.h;
	L_t = NaN(Nt, 1); L_t(1) = Ltm1;
	Qs_t = NaN(Nt, 1); Qs_t(1) = 0;
	Qd_t = NaN(Nt, 1); Qd_t(1) = 0;
	Qt_t = NaN(Nt, 1); Qt_t(1) = 0;
	c_t = NaN(Nt, 1);

	F_L_tau = F.F_L_tau;
	f_L_tau = F.f_L_tau;

	% prepare for saving
	saveTauStep = 6*ceil(1/dt); % save tau at 6h resolution
	saveTStep = ceil(2/dt); % save t at 2h resolution
	NSaveT = ceil(Nt / saveTStep); % how many timesteps in saved matrices
	NSaveTau = ceil(Ntau / saveTauStep); % how many tausteps in saved matrices
	saveCount = 0; % count timesteps; every saveTStep we save a new row
	saveId = 0; % matrices are saved one row at a time; this is the last saved row

	F_t_tau = NaN(NSaveT, NSaveTau);
	dFdt_t_tau = NaN(NSaveT, NSaveTau);
	dFdtau_t_tau = NaN(NSaveT, NSaveTau);

	

	% prepare to show iteration progress
	progressFormat = '% 6u/% 6u';
	progressBack = repmat('\b', 1, 13);
	fprintf(progressFormat, Nt, Nt);

	% solve each timestep
	for t = 1:Nt-1
		if saveCount == saveTStep % if we save F ecc
			% update save counter
			saveId = saveId + 1;
			saveCount = 1;

			% calc Q at time t
			[Qt_t(t), Qs_t(t), Qd_t(t), c_t(t), saveF_tau, savedFdt_tau, savedFdtau_tau] = ...
				convolutionQ(time, dt, t, tau, dtau, Ntau, L_t, dL, NL, Lmin, f_L_tau, F_L_tau, h_t, uh, alphaFun, baseflow, Ltm1, L_t(t));

			% also save matrices
			F_t_tau(saveId, :) = saveF_tau(1:saveTauStep:end);
			dFdt_t_tau(saveId, :) = savedFdt_tau(1:saveTauStep:end);
			dFdtau_t_tau(saveId, :) = savedFdtau_tau(1:saveTauStep:end);

			% update progress
			fprintf([progressBack progressFormat], Nt-t, Nt);

		else % if we don't save F
			% update save counter
			saveCount = saveCount + 1;

			% calc Q
			[Qt_t(t), Qs_t(t), Qd_t(t), c_t(t)] = ...
				convolutionQ(time, dt, t, tau, dtau, Ntau, L_t, dL, NL, Lmin, f_L_tau, F_L_tau, h_t, uh, alphaFun, baseflow, Ltm1, L_t(t));
		end

		% update length variables
		L_t(t+1) = LQ.fun(Qs_t(t));
		Ltm1 = L_t(t);
	end
	fprintf(progressBack);

	% save results back to timeseries
	result_t.L = L_t;
	result_t.Qs = Qs_t;
	result_t.Qd = Qd_t;
	result_t.Q = Qt_t;
	result_t.c = c_t;

	result_t_tau = struct();
	result_t_tau.t = time(1:saveTStep:end);
	result_t_tau.tau = tau(1:saveTauStep:end);
	result_t_tau.F_t_tau = F_t_tau;
	result_t_tau.dFdt_t_tau = dFdt_t_tau;
	result_t_tau.dFdtau_t_tau = dFdtau_t_tau;

	% fprintf('\n');
end

function [Qt_Lt, Qs_Lt, Qd_Lt, corr_Lt, F_Lt_tau, dFdt_Lt_tau, dFdtau_Lt_tau] = convolutionQ(time, dt, t, tau, dtau, Ntau, L_t, dL, NL, Lmin, f_L_tau, F_L_tau, h_t, uh, alphaFun, baseflow, Ltm1, Lt)
	
	% calcolo f(tau, t) ed F(tau, t)
	iL = max((Lt - Lmin) / dL + 1, 1);
	iLsx = floor(iL);
	iLdx = min(iLsx + 1, NL);

	dFdtau_Lt_tau = f_L_tau(iLsx, :) .* (iLdx - iL) + f_L_tau(iLdx, :) .* (iL - iLsx);
	F_Lt_tau = F_L_tau(iLsx, :) .* (iLdx - iL) + F_L_tau(iLdx, :) .* (iL - iLsx);

	% calcolo F(tau, t-1) e dFdt
	iL = (Ltm1 - Lmin + 1e-4) / dL + 1; % 1e-4 acts as tolerance when we're close to min length
	iLsx = floor(iL);
	iLdx = min(iLsx + 1, NL);

	F_tm1_tau = F_L_tau(iLsx, :) .* (iLdx - iL) + F_L_tau(iLdx, :) .* (iL - iLsx);
	dFdt_Lt_tau = (F_Lt_tau - F_tm1_tau) / dt;

	% definiamo Lstar e tstar (per ogni Ltest)
	Lstar_Lt_t = Lt - uh * ( time(1:t-1)' - time(t) );

	NLt = numel(Lt);
	tstar_Lt = NaN(NLt, 1);
	corr_Lt = NaN(NLt, 1);

	for i = 1:NLt
		id = find(L_t(1:t-1)' - Lstar_Lt_t(i,:) >= 0, 1, 'last');
		
		if isempty(id)
			tstar_Lt(i) = NaN;
			corr_Lt(i) = 1;
		elseif ( Lt(i) - Ltm1) / dt < -uh
			tstar_Lt(i) = time(id);
			corr_Lt(i) = 2;
		else
			tstar_Lt(i) = time(id);
			corr_Lt(i) = 3;
		end
	end
	
	% define mask (which Ltest and which tau need which correction)
	mask_Lt_tau = ones(NLt, Ntau); % base case: no correction
	mask_Lt_tau(corr_Lt == 2, :) = 2; % fast contraction correction
	mask_Lt_tau(corr_Lt == 3, :) = 3; % counter-correction
	% TODO smaller taus dont have correction 3
	
	% define G
	Gs_Lt_tau = dFdtau_Lt_tau;

	alpha = alphaFun(Lt, tau');
	Gd_Lt_tau = NaN(NLt, Ntau);
	Gd_Lt_tau(mask_Lt_tau == 1) = dFdt_Lt_tau(mask_Lt_tau == 1);
	Gd_Lt_tau(mask_Lt_tau == 2) = - alpha(mask_Lt_tau == 2) .* dFdtau_Lt_tau(mask_Lt_tau == 2);
	Gd_Lt_tau(mask_Lt_tau == 3) = - alpha(mask_Lt_tau == 3) .* dFdtau_Lt_tau(mask_Lt_tau == 3); 
	% TODO: add third correction term

	% convolute G and rain
	tmentau = time(t) - tau;
	J_tau = interp1(time, h_t, tmentau, 'previous', 0);

	Qs_Lt = Gs_Lt_tau * J_tau * dtau / dt + baseflow; % mm/h
	Qd_Lt = Gd_Lt_tau * J_tau * dtau / dt;
	Qt_Lt = Qs_Lt + Qd_Lt; % mm/h
end