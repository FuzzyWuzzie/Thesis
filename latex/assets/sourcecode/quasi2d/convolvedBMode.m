% generate normally distributed noise across the domain
rng(domainSeed);
backgroundMap = 1 + noiseMagnitude * randn([Ny, Nx]);

% transform the background map into the domain -1 -> 1
backgroundMap = (backgroundMap - min(min(backgroundMap))) / (max(max(backgroundMap)) - min(min(backgroundMap)));
backgroundMap = 2 * backgroundMap - 1;

% if a finite-element model of tissue compression is being used,
% compress the background map
if model ~= 0
	% extract the resultant degrees of freedom from the model
	[x, y, u, v] = extractUV([transducerWidth, depth], size(backgroundMap), model);

	% generate the domain over which interpolation will occur
	[xx, yy] = meshgrid(x, y);

	% interpolate the data to deform it
	backgroundMap = interp2(xx, yy, backgroundMap, xx - u, yy - v, 'spline', mean(mean(backgroundMap)));
end

% use a cosine function to create the point-spread function shape in the axial direction
xpsf = linspace(-windowWidth / 2, windowWidth / 2, 1 * pointsPerWaveLength);
ypsf = linspace(0, 4 * waveLength, 2 * pointsPerWaveLength);
[xmpsf, ympsf] = meshgrid(xpsf, ypsf);
psf = cos(2 * pi * frequency * ympsf / waveSpeed);

% apply a lateral gaussian "filter" to the signal
mu = 0;
sigma = windowWidth / 4;
gauss = (1 / (sigma * sqrt(2 * pi))) * exp(-(xmpsf - mu) .^ 2 / (2 * sigma ^ 2));
psf = psf .* gauss;

% apply an axial gaussian "filter" to the signal
mu = 2 * waveLength;
sigma = waveLength * 2;
gauss = (1 / (sigma * sqrt(2 * pi))) * exp(-(ympsf - mu) .^ 2 / (2 * sigma ^ 2));
gauss = (gauss - min(min(gauss))) / (max(max(gauss)) - min(min(gauss)));
psf = psf .* gauss;

% normalize it to -1 -> 1
psf = (psf - min(min(psf)))/(max(max(psf)) - min(min(psf)));
psf = 2 * psf - 1;

% convolve the scattering map with the point spread function
bmode = conv2(backgroundMap, psf);

% crop the image to the appropriate size
sz = size(bmode);
bmode = bmode(int32(floor((sz(1) - pointDepth) / 2)) : int32(floor((sz(1) - pointDepth) / 2) + pointDepth - 1), int32(floor((sz(2) - numElements) / 2)) : int32(floor((sz(2) - numElements) / 2) + numElements - 1));

% apply classical ultrasound post-processing
% the signal is currently oscillating a high frequency
%  - extract the envelope of the signal
bmode = envelopeDetection(bmode')';
% apply log compression to allow to be readily viewed
bmode = logCompression(bmode, 3, true);

% normalize it to 0 -> 1
bmode = (bmode - min(min(bmode))) ./ (max(max(bmode)) - min(min(bmode)));