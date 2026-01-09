function timeseries = generateRainfallPoisson(domainT, rain)
	% fprintf('Generating poisson rainfall...\n');

	% get parameters
	Nt = domainT.Nt; % number of timesteps
	dt = domainT.dt; % [h] time step
	t = domainT.t; % [h] time domain

	% allocate table
	timeseries   = table();
	timeseries.t = t'; % [h] time domain
	timeseries.h = zeros(Nt, 1); % [mm] rainfall depth

	% generate daily rain
	NtDaily = ceil(domainT.D/24) + 1; % number of days in simulation
	hDaily = rand(NtDaily, 1) < rain.lambda; % rainy days
	hDaily = double(hDaily);
	hDaily(hDaily == 1) = exprnd(rain.alpha, sum(hDaily == 1), 1); % [mm] random exponentially distributed daily rain

	timeseries.hDaily = hDaily(floor(t/24) + 1) * dt/24; % [mm] daily rainfall save at dt step

	% downscale
	switch rain.downscale
		case "linear"
			timeseries.h = movmean(timeseries.hDaily, ceil(24/dt)); % mm
		case "triangular"
			tau = -0.5:dt/24:1.5;
			shift = ceil(24/dt);
			kernel = (1 - abs(tau - 0.5)) *dt/24;

			% tau = -0.25:dt/24:0.25;
			% shift = ceil(24/dt/4);
			% kernel = (0.5 - 2*abs(tau))*8 * dt/24;

			h = conv(timeseries.hDaily, kernel, "full");
			timeseries.h = h(shift + (1:Nt));
		case "conservative"
			% span will start as 1 day (rainfall resolution) and be divided over and over until it becomes the desired dt
			if rem(24, dt) ~= 0
				error("Timestep must be a finite fraction of a day (e.g. whole minutes)");
			end

			span = 24/dt; % how many timesteps in one span (initially one day)
			nspan = NtDaily; % how many spans in simulation
			factors = factor(span);

			h = hDaily;
			for s = 1:numel(factors)
				% generate random factor from beta distribution (note: if a=1, then it's uniform in binary downscaling)
				alpha = rain.downscaleBeta * ones(factors(s), nspan);
				randFactors = gamrnd(alpha, alpha, factors(s), nspan);
				randFactors = randFactors ./ sum(randFactors);

				% calc max factor (for each timestep) to ensure max downscaled rain is 80 mm/h
				limitRandFact = rain.downscaleLimit * span*dt/factors(s)./h';
				limitRandFact(limitRandFact > 1) = 1;

				% find factors exceeding the limit
				maxRandFact = max(randFactors);
				exceedingCol = maxRandFact > limitRandFact;

				% in those cases, make the downscaling more equal by a weight omega
				% (omega = 0 -> use random factors; omega = 1 -> split rain equally; omega is calculated such that the
				% max downscaled rain is exactly equal to the limit)
				omega = (1 - factors(s) * limitRandFact) ./ (1 - factors(s) * maxRandFact);
				omega = omega .* rand(1, nspan);
				randFactors(:, exceedingCol) = omega(exceedingCol)		.* randFactors(:, exceedingCol) + ...
											   (1 - omega(exceedingCol)) / factors(s);

				span = span / factors(s);
				nspan = nspan * factors(s);

				h = repelem(h, factors(s), 1) .* randFactors(:); % mm

				if max(h) > rain.downscaleLimit*dt*span + 1e4 || min(h) < 0
					error("Error while limiting rain intensity during downscaling. Try again.");
				end
			end

			timeseries.h = h(1:Nt);

		otherwise
			warning("Unrecognized rainfall downscaling type. Using daily rainfall");
			timeseries.h = timeseries.hDaily;
	end
end