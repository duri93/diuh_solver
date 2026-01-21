function par = makeDependentParameters(par, exclude)
	% Updates parameter struct to ensure dependent parameters are congruent with the associated independent parameters.
	% par: the parameters struct to be updated
	% exlude: an array of strings of parameters you don't want to update. Use with caution, this may introduce fatal errors.

	% define default exclude
	if ~exist("exclude", "var")
		exclude = "";
	end

	% define dependencies
	[~, ~, dependencies] = definitions('pulse');

	% filter out the ones that are excluded
	id = ismember([dependencies{:,1}]', exclude);
	dependencies = dependencies(~id, :);

	% update parameters using dependencies
	for d = 1:height(dependencies)
		names = split(dependencies{d,1}, ".");
		par.(names(1)).(names(2)) = dependencies{d,2}(par);
	end
end