function timeseries = generateRainfallPulse(domainT, rain)
	% fprintf('Generating pulse rainfall...\n');

	% get parameters
	Nt = domainT.Nt; % number of timesteps
	dt = domainT.dt; % [h] time step
	t = domainT.t';  % [h] time domain

	start = rain.start; % [h] start time of rainfall pulse
	duration = rain.duration; % [h] duration of rain pulse
	volume = rain.volume; % [h] total rain volume

	% allocate table
	timeseries   = table();
	timeseries.t = t; % [h] time domain
	timeseries.h = zeros(Nt, 1); % [mm] rainfall depth

	% define rainfall pulse
	id = timeseries.t >= start & timeseries.t < start + duration;
	timeseries.h(id) = volume / sum(id); % [mm] rainfall only during pulse

	% fix shape
	switch rain.shape
		case "rectangular"
		case "triangular"
			window = ceil(duration / dt);
			movmean(timeseries.h, window);
		otherwise
			error("Parameter rain.shape must be either ""rectangular"" or ""triangular"".");
end
