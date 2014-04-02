% extract the displacement along a lateral line traversing the focal point
% throughout the complete simulation time
[x, t, lateralCut] = getLateralFocalCutDisplacement(model);

% use the mean value of the displacement through time and location
% as the isoline value to track
targetValue = mean(mean(lateralCut));

% loop along the x-coordinates of the dataset
points = NaN(length(x), 2);
% ignore the first 4 data points which lie within the width of the
% initial radiation force
for xi = 5 : length(x)
	% find the point in time where the displacement at the given x-coordinate
	% crosses the target value previously established
	for ti = 2 : length(t)
		if (lateralCut(ti, xi) < targetValue) && (lateralCut(ti - 1, xi) >= targetValue)
			% if the cross-over point was found, store it and continue with
			% the next x-coordinate
			points(xi, :) = [t(ti), x(xi)];
			ti = length(t);
			break;
		end
	end
end

% remove NaNs from the dataset
points = points(~any(isnan(points), 2), :);

% differentiate to get shear velocity
% (use a center-weighted moving window average filter first
% otherwise the data will become unusable)
Ct = diff(smooth(points(:, 2))) ./ diff(smooth(points(:, 1)));
x = linspace(min(points(:, 2)), max(points(:, 2)), length(Ct));