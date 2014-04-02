% calculate the best grid size for computational efficiency
% taking into account Nyquist so that appropriate sampling is used
dx = medium.sound_speed / (2 * probingFrequency);
dy = dx;
Nx = bestFactor(domainDepth / dx + (2 * PML_X_SIZE)) - (2 * PML_X_SIZE);
Ny = bestFactor(domainWidth / dy + (2 * PML_Y_SIZE)) - (2 * PML_X_SIZE);
dx = domainDepth / Nx;
dy = domainWidth / Ny;

% make the grid!
kGrid = makeGrid(Nx, dx, Ny, dy);

% create the time array
kGrid.t_array = makeTime(kGrid, medium.sound_speed, [], pulseCycles / probingFrequency);

% beam-form the data
% create indices for our elements
numActiveElements = Ny;
elementIndices = -(numActiveElements - 1) / 2 : (numActiveElements - 1) / 2;

% calculate time delays for a steered and focussed beam
elementWidth = dy;
delayTimes = focalDepth / medium.sound_speed * (1 - sqrt(1 + (elementIndices * elementWidth ./ focalDepth) .^ 2 )); % [s]

% convert the delays to be in units of time points
delayTimes = delayTimes - min(delayTimes);
delayTimes = delayTimes ./ kGrid.dt;

% use the k-wave toolbox's function "toneBurst" to create the signal
inputSignal = toneBurst(1 / kGrid.dt, probingFrequency, pulseCycles, 'SignalOffset', delayTimes);

% scale the signal by the source pressure
source.p = sourcePressure .* inputSignal;

% truncate the input signal to the appropriate length
source.p = source.p(:, 1:length(kGrid.t_array));

% make only the nodes along the top boundary apply pressure to the domain
source.p_mask = zeros([Nx, Ny]);
source.p_mask(1, 1 : (numActiveElements)) = 1;

% tell the k-wave toolbox to record the pressure and intensity for
% the entire domain, continuously
sensor.mask = ones(Nx, Ny);
sensor.record = {'I', 'p'};

% set up simulation settings
inputArgs = {'PlotSim', false, 'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE], 'DataCast', 'single', 'DataRecast', true, 'DisplayMask', 'off'};

% run the simulation
sensorData = kspaceFirstOrder2D(kGrid, medium, source, sensor, inputArgs{:});

% reshape the output data to make sense
Ix = reshape(sensorData.Iy, [Nx, Ny, kGrid.Nt]);
Iy = reshape(sensorData.Ix, [Nx, Ny, kGrid.Nt]);
P = reshape(sensorData.p, [Nx, Ny, kGrid.Nt]);

% calculate body forces
alpha = (attenuationCoefficient * 100 * probingFrequency / 1e6) / (20 / log(10));
Fx = 2 * alpha .* Ix ./ soundSpeed;
Fy = -2 * alpha .* Iy ./ soundSpeed;