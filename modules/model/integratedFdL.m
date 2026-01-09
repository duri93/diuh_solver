% TODO
% - make integration resolution (dL_i and dtau_i) either explicit parameter or automatically set
% - OR BETTER: use matlab functions to solve the equation!

function F = integratedFdL(domainTau, domainL, alpha, f0)
	fprintf('Integrating df/dL...\n');

	% get params
	uh = f0.uh; % [m/h] hyporheic velocity in dry network
	alpha_fun = alpha.fun; % impact coefficient



	% integration domain over L
	Lmin = domainL.Lmin; % [m] minimum network length
	Lmax = domainL.Lmax; % [m] maximum network length
	dL_i = 1; % [m] resolution for integration 
	L_i = Lmin:dL_i:Lmax; % [m] integration length domain
	NL_i = numel(L_i); % number of length steps

	% integration domain over tau
	tauMin = 0; % [h] min tau
	tauMax = domainTau.tauMax; % [h] max tau
	dtau_i = 0.2; % [h] resolution for integration
	tau_i = tauMin:dtau_i:tauMax; % [h] integration tau domain
	Ntau_i = numel(tau_i); % number of tau steps
	


	% allocate f matrix
	f_matrix_i = zeros(NL_i, Ntau_i);

	% set boundary condition
	if f0.L0 > Lmax || f0.L0 < Lmin
		error("Boundary condition of f must be in the length domain (i.e., par.domainL.Lmin <= par.f0.L0 <= par.domainL.Lmax)");
	end

	L0 = f0.L0;
	id_L0 = find(L_i >= L0, 1);
	f_matrix_i(id_L0, :) = f0.fun(tau_i);

	% integrate forwards (from L0 to Lmax)
	f_matrix_i = integrate_f(f_matrix_i, L_i, dL_i, id_L0,  1, NL_i, tau_i, dtau_i, Ntau_i, alpha_fun, uh);

	% integrate backwards (from L0 to Lmin)
	f_matrix_i = integrate_f(f_matrix_i, L_i, dL_i, id_L0, -1, 1,    tau_i, dtau_i, Ntau_i, alpha_fun, uh);
	


	% define tau saving domain
	dtau = domainTau.dtau;
	tau = domainTau.tau;
	Ntau = domainTau.Ntau;

	% define L saving domain
	NL = domainL.NL;
	L = domainL.L;
	dL = domainL.dL;

	% resample f_matrix_i from integration resolution to saving resolution
	f_matrix = interp2(tau_i', L_i, f_matrix_i, tau', L, 'linear', 0);
	


	% accumulate f over tau to get F
	F_matrix = cumtrapz(tau, f_matrix, 2);

	% check tauMax is enough to deplete all the IUH
	totalF = min(F_matrix(:, end));
	if totalF < 0.995
		fprintf("WARNING: F(tau = tauMax, L) = %.4f is less than 0.995 for some network lengths.\n", totalF)
		fprintf("Consider increasing domainTau.tauMax to avoid mass balance problems.\n");
	end
	
	

	% save results
	F = struct();
	F.Ntau = Ntau;
	F.tau = tau';
	F.dtau = dtau;
	F.NL = NL;
	F.L = L';
	F.dL = dL;
	F.F_L_tau = F_matrix;
	F.f_L_tau = f_matrix;
end

function f_matrix_i = integrate_f(f_matrix_i, L_i, dL_i, iStart, iStep, iStop, tau_i, dtau_i, Ntau_i, alpha_fun, uh)
		f_prev = f_matrix_i(iStart, :);

		for i = iStart:iStep:iStop-iStep
			F = alpha_fun(L_i(i), tau_i) .* f_prev;

			f_next = f_prev;

			for j = 2:Ntau_i-1
				F_plus = F(j); % F_{j+1/2}
				F_minus = F(j-1); % F_{j-1/2}

				RHS = (F_plus - F_minus) / dtau_i;

				f_next(j) = f_prev(j) + iStep * dL_i / uh * RHS;
			end

			f_next(1) = 2 * f_next(2) - f_next(3); % linear extrapolation for first value, starting from values 2 and 3
			f_next(end) = f_prev(end); % constant value at end

			f_matrix_i(i + iStep, :) = f_next;
			f_prev = f_next;
		end
end

% integrate backwards (from max L)
	% for i = NL_i:-1:2
	% 
	% 	Li = L_i(i);  
	% 	g_vec = alpha_fun(Li, tau_i);  
	% 	F = g_vec .* f_prev;  % flusso
	% 
	% 	f_next = f_prev;  
	% 
	% 	% scelta dei valori a j+1/2 e j-1/2: cella precedente a
	% 	for j = 2:Ntau_i-1
    % 		F_plus  =  F(j);    % F_{j+1/2}
    % 		F_minus = F(j-1);   % F_{j-1/2}
    % 		RHS = (F_plus - F_minus) / dtau_i;
    % 		f_next(j) = f_prev(j) + (-dL_i / uh) * RHS;
	% 	end
	% 
	% 	f_next(1) = 2 * f_next(2) - f_next(3); 
	% 	%estrapolazione lineare in avanti, la soluzione iniziale viene trovata
	% 	%tramite una retta che ha la medesima pendenza che la curva ha appena
	% 	%dentro al dominio
	% 	f_next(end) = f_prev(end);
	% 
	% 	f_matrix_i(i-1, :) = f_next;
	% 	f_prev = f_next;
	% end
