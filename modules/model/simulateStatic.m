function static_result_t = simulateStatic(result_t, F_L_tau, baseflow, Lstatic)
	Nt = height(result_t);

	static_result_t = result_t;
	static_result_t.L = repmat(Lstatic, Nt, 1);
	static_result_t.Qd = zeros(Nt, 1);
	static_result_t.c = ones(Nt, 1);

	% get iuh for Lstatic
	iuh = interp2(F_L_tau.tau, F_L_tau.L, F_L_tau.f_L_tau, F_L_tau.tau, Lstatic); % [1/h]

	% convolution for Q
	q = conv(static_result_t.h, iuh', "full");

	static_result_t.Qs = q(1:Nt) + baseflow;
	static_result_t.Q = static_result_t.Qs;
end