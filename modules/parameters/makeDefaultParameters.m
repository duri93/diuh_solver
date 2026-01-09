function par = makeDefaultParameters(rainType, baseCatchment)
	% rainType: type of rainfall. "pulse" or "poisson"
	% baseCatchment: specify a catchment name to use its empirical parameters

	% define catchment combination if not specified
	if ~exist("baseCatchment", "var")
		baseCatchment = "HubbardBrook25";
	end
	if ~exist("rainType", "var")
		rainType = "poisson";
	end

	% define parameters
	par = definitions(rainType);
	par = updateByCatchment(par, baseCatchment);
	par = makeDependentParameters(par);
end

function par = updateByCatchment(par, baseCatchment)
	[~, catchments] = definitions();

	% find catchment
	id = find(catchments.name == baseCatchment, 1);

	% check specified one exists
	if isempty(id)
		id = 1;
		warning("Specified catchment is not defined. Falling back to default parameters (%s)", catchments.name(id));
	end

	% update params
	par.baseCatchment = catchments.name(id);
	par.domainL.A = catchments.A(id);
	par.domainL.Dg = catchments.Dg(id);
	par.domainL.phiC = catchments.phiC(id);
	par.LQ.k = catchments.LQ_k(id);
	par.LQ.theta = catchments.LQ_theta(id);
end