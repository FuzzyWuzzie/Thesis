function bmode = convolvedBMode(dbConnection, caseIndex, model)

	% build our point spread function
	query = [...
		'select waveSpeed, '...
		'domainWidth, '...
		'domainDepth, '...
		'samplingFrequency, '...
		'probingFrequency, '...
		'numElements, '...
		'numActiveElements, '...
		'domainSeed, '...
		'noiseMagnitude, '...
		'Nx, '...
		'Ny '...
		'from convBModes where id=''' sprintf('%d', caseIndex) ''' limit 1;'];
	results = fetch(dbConnection, query);

	waveSpeed = results.waveSpeed;
	samplingFrequency = results.samplingFrequency;
	depth = results.domainDepth;
	pointsPerM = samplingFrequency / waveSpeed;
	pointDepth = pointsPerM * depth;
	frequency = results.probingFrequency;
	waveLength = waveSpeed / frequency;
	pointsPerWaveLength = round(pointsPerM * waveLength);
	numElements = results.numElements;
	numActiveElements = results.numActiveElements;
	transducerWidth = results.domainWidth;
	windowWidth = numActiveElements * transducerWidth / numElements;

	% interpolate our data for higher fidelity
	[xx, yy] = meshgrid(linspace(-results.domainWidth/2, results.domainWidth/2, results.Nx), linspace(0, depth, results.Ny));
	[x, y] = meshgrid(linspace(-results.domainWidth/2, results.domainWidth/2, numElements), linspace(0, depth, pointDepth));

	% get our scatterer map
	rng(results.domainSeed);
	backgroundMap = 1 + results.noiseMagnitude * randn([results.Ny, results.Nx]);

	backgroundMap = interp2(xx, yy, backgroundMap, x, y, 'cubic');
	backgroundMap = (backgroundMap - min(min(backgroundMap)))/(max(max(backgroundMap)) - min(min(backgroundMap)));
	backgroundMap = 2*backgroundMap-1;

	% if we have a model input, deform our background map
	if model ~= 0
		% extract our U and V from the model
		[x, y, u, v] = extractUV([transducerWidth, depth], size(backgroundMap), model);

		% and store it in the db!
		saveMatrixInDB(dbConnection, 'convBModes', caseIndex, u, 'uMap');
		saveMatrixInDB(dbConnection, 'convBModes', caseIndex, v, 'vMap');

		backgroundMap = compressBMode(x, y, backgroundMap, u, v);
	end

	% use a cos function to create the shape in the axial direction
	xpsf = linspace(-windowWidth/2, windowWidth/2, 1*pointsPerWaveLength);
	ypsf = linspace(0, 4*waveLength, 2*pointsPerWaveLength);
	[xmpsf, ympsf] = meshgrid(xpsf, ypsf);
	psf = cos(2 * pi * frequency * ympsf / waveSpeed);
	
	% lateral gaussian
	mu = 0; sigma = windowWidth/4;
	gauss = (1 / (sigma * sqrt(2*pi))) * exp(-(xmpsf - mu).^2 / (2*sigma^2));
	psf = psf .* gauss;

	% axial gaussian
	mu = 2*waveLength; sigma = waveLength*2;
	gauss = (1 / (sigma * sqrt(2*pi))) * exp(-(ympsf - mu).^2 / (2*sigma^2));
	gauss = (gauss - min(min(gauss)))/(max(max(gauss)) - min(min(gauss)));
	psf = psf .* gauss;

	% normalize it
	psf = (psf - min(min(psf)))/(max(max(psf)) - min(min(psf)));
	psf = 2*psf-1;

	% now convolve!
	bmode = conv2(backgroundMap, psf);

	% and window it out
	sz = size(bmode);
	bmode = bmode(int32(floor((sz(1) - pointDepth) / 2)) : int32(floor((sz(1) - pointDepth) / 2) + pointDepth - 1), int32(floor((sz(2) - numElements) / 2)) : int32(floor((sz(2) - numElements) / 2) + numElements - 1));

	% and post-process it
	bmode = envelopeDetection(bmode')';
	bmode = logCompression(bmode, 3, true);

	% normalize it
	bmode = (bmode - min(min(bmode))) ./ (max(max(bmode)) - min(min(bmode)));

end