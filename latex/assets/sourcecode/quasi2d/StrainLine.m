% calculate the number of windows along this scanline
M = maxX / deltaAx;

% these variables will return the alpha and tau coefficients
% along this scanline
alpha = [];
tau = [];

% loop through every window along this scanline
for m = 1:M
	% calculate the locations of the two windows to start the search from
	pr1 = m .* deltaAx;
	pr2 = sum(deltaAx ./ alpha);

	% extract the window in the uncompressed image
	x1 = linspace(pr1, pr1 + L, numPoints);
	r1 = interp1(1 : length(I1), I1, x1, 'linear', 'extrap');

	% initialize search variables
	alphas = maxAlpha:-0.001:1;
	taus = zeros(size(alphas));
	correlations = zeros(size(alphas));
	
	% loop through the possible values of alpha
	for i = 1 : length(alphas)
		% keep track of which scanline gives us the best correlation
		columnCorrelations = [];

		% loop through the adjacent scanlines
		for c = (column - columnRadius) : (column + columnRadius)
			% make sure we have a valid scanline
			if (c < 1) || (c > columns)
				continue;
			end
			
			% extract the window in the compressed image
			x2 = linspace(round(pr2), round(pr2) + round(L), numPoints);
			r2 = interp1(1:length(I2(:, c)), I2(:, c), x2, 'linear', 'extrap');
			
			% get the correlation for the two windows
			columnCorrelations = [columnCorrelations, Correlate(r1, r2)];
		end
		
		% pick the scanline that had the best correlation for this alpha
		c = -columnRadius : column + columnRadius;
		[correlations(i), mindex] = min(columnCorrelations);
		taus(i) = c(mindex);
	end

	% employ a-priori smoothing to the data
	% in case of errant outliers
	[~, mindex] = min(correlations);
	if abs(alphas(mindex) - alpha(length(alpha))) > 0.02
		alpha = [alpha; alpha(length(alpha))];
		tau = [tau; tau(length(tau))];
	else
		alpha = [alpha; alphas(mindex)];
		tau = [tau; taus(mindex)];
	end
end
