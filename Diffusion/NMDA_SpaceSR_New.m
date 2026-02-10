%% NMDA_SpaceSR.m
%   Spatiotemporal glutamate simulation + NMDA state occupancy maps
%   (refactored – 30-Jun-2025)

%% ---------------------------- User parameters ---------------------------
TimeGlu = 8;                        % Duration of glutamate diffusion (ms)
BinGlu = 100;                        % Number of bins for distance
ScaleParameter = 1;                  % Scaling factor for NMDA computation
timeCalculus = 200 * ScaleParameter; % Basic 200 ms of NMDA ODE integration horizon (ms)
inputPattern = 'DistanceFree*.txt';  % Input file pattern

% Optimized distance calculation
DistanceReal = linspace(1/BinGlu, 1 + 1/BinGlu, BinGlu + 1);
BinTotal = round(BinGlu*timeCalculus/TimeGlu);
NMDATimeBin = linspace(0, timeCalculus, BinTotal); % Generate t for f

% Basic upload files inputPattern
InitialDistribution = processDistanceFreeFiles();

%NormParameter = 20*4*10^19/(6.023*10^23);
NormParameter = 10^21/(6.023*10^23);
Volume = ((4/3) * pi * ((2 * DistanceReal(2:end)').^3 - (2 * DistanceReal(1:end-1)').^3));
%InitialDistribution(:, 2:end) = NormParameter * InitialDistribution(:, 2:end) ./ Volume;
ConcentrationOneMolecules = 1.66 * 10^-6; % mM 1 molecules creat such concentration in mM
Con= (1000)*ConcentrationOneMolecules; % uM
InitialDistribution(:, 2:end) = Con*InitialDistribution(:, 2:end) ./ Volume;
%%
% Plot glutamate concentration
figure(1)
data = InitialDistribution(:, 2:end);
[nRows, nCols] = size(data);
contourf(log(data), 10, 'LineColor', 'none')
colorbar
clim([-8, 4.5])

set(gca, 'XTick', linspace(1, nCols, TimeGlu), ...
         'XTickLabel', 0:TimeGlu-1, ...
         'YTick', linspace(1, nRows, 5), ...
         'YTickLabel', linspace(0, 2, 5))

xlabel('Time (ms)')
ylabel('Distance (um)')
title('Glu concentration (uM)')

% Save glutamate concentration data
save('glutamate_concentration_data.mat', 'data');
fprintf('Glutamate concentration data saved to glutamate_concentration_data.mat\n');

%% *********************Computation NMDA open probability  ****************
GluTotal=zeros(1,BinTotal);
TwoBount=zeros(BinGlu, BinGlu);
TwoBountRaw=1:1:BinGlu;
% Maintain original solver but enhance numerical handling
% Pre-compute constants outside loop for efficiency
options = odeset('RelTol', 1e-12, ...      % Balanced precision for solver
                'AbsTol', 1e-20, ...       % Match concentration scale
                'NonNegative', 1:5, ...    % Physical constraint (all states non-negative)
                'Refine', 10, ...          % Increase output resolution
                'MaxStep', 0.1);           % Prevent over-aggressive time stepping

% Pre-allocate interpolation target times
interp_times = round(ScaleParameter * TwoBountRaw);
eps_threshold = eps;                       % Cache eps value
realmin_threshold = realmin('double');     % Cache realmin value
zero_result = zeros(BinGlu, 1);           % Pre-allocate zero array

% Main processing loop
for Jitter = 1:BinGlu
    fprintf('Processing iteration %d of %d\n', Jitter, BinGlu);
    
    % Extract and prepare glutamate test data
    TestGlu = InitialDistribution(Jitter, 2:end);
    TestGlu(numel(GluTotal)) = 0;          % Reset final element
    
    % Define ODE system and solve
    ode = @(t,y) NMDA(t, y, TestGlu', NMDATimeBin, timeCalculus);
    [T, Y] = ode45(ode, [0 timeCalculus], [1 0 0 0 0], options);
    
    % Extract state variable of interest (4th component)
    TempY = Y(:, 4);
    valid_mask = TempY > eps_threshold;    % Identify numerically valid data points
    
    if any(valid_mask)
        % Log-space interpolation for numerical stability
        log_data = log10(TempY(valid_mask));
        TempY_interp = interp1(T(valid_mask), log_data, interp_times, ...
                              'pchip', -inf);  % PCHIP with -inf extrapolation
        
        % Convert back from log space
        TempY_result = 10.^TempY_interp;
        
        % Apply physical floor constraint
        TempY_result(TempY_result < realmin_threshold) = 0;
    else
        % Handle case with no valid data points
        TempY_result = zero_result;
    end
    
    % Store results
    TwoBount(Jitter, :) = TempY_result;
end

% Save NMDA open state data
save('NMDA_open_state_TwoBount.mat', 'TwoBount');
fprintf('NMDA open state data saved to NMDA_open_state_TwoBount.mat\n');

axis manual
figure(2);                  % display image
%axis([0 timeCalculus 0 1])

subplot(2, 4, 1); plot(T,Y(:,1))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 2); plot(T,Y(:,2))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 3); plot(T,Y(:,3))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 4); plot(T,Y(:,4))
%axis([0 timeCalculus 0 0.02])

subplot(2, 4, 5); plot(T,Y(:,5))
%axis([0 timeCalculus 0 0.01])
MinPlot=min(min(TwoBount(:,1:end)));
MaxPlot=max(max(TwoBount(:,1:end)));
figure(3)
%TwoBount(1,15)=0.49; SCALING MAX
contourf((TwoBount(:,1:end)), 20, 'LineColor','none')
colorbar
clim([0,0.4]);
ax = gca;
ax.XTick = linspace(1, size(TwoBount, 2), 5);
ax.XTickLabel = linspace(0, timeCalculus, 5);
ax.YTick = linspace(1, size(TwoBount, 1), 5);
ax.YTickLabel = linspace(0, 2, 5);
xlabel('Time (ms)')
ylabel('Distance (um)')
title(' Open state')



% contourf((TwoBount(:,1:end)), 20, 'LineColor','none')
% colorbar;
% clim([0,0.4]);
%caxis([MinPlot MaxPlot]);

function InitialDistribution = processDistanceFreeFiles(myFolder)
% Process all DistanceFree*.txt files and compute average
% Input: myFolder - path to folder (optional, defaults to current directory)
% Output: InitialDistribution - averaged data from all files

% Set default folder if not provided
if nargin < 1
    myFolder = pwd;
end

% Validate folder existence
if ~isfolder(myFolder)
    error('Folder does not exist: %s', myFolder);
end

% Get file list efficiently
filePattern = fullfile(myFolder, 'DistanceFree*.txt');
theFiles = dir(filePattern);

% Check if files exist
numFiles = length(theFiles);
if numFiles == 0
    warning('No files matching pattern found: %s', filePattern);
    InitialDistribution = [];
    return;
end

% Pre-allocate and process files
fprintf('Processing %d files...\n', numFiles);
InitialDistribution = 0;

try
    for k = 1:numFiles
        fullFileName = fullfile(myFolder, theFiles(k).name);
        
        % Load and accumulate data
        thisStructure = load(fullFileName);
        InitialDistribution = InitialDistribution + thisStructure;
        
        % Progress indicator for large datasets
        if mod(k, 10) == 0 || k == numFiles
            fprintf('Processed %d/%d files\n', k, numFiles);
        end
    end
    
    % Compute average
    InitialDistribution = InitialDistribution / numFiles;
    fprintf('Successfully processed and averaged %d files\n', numFiles);
    
catch ME
    error('Error processing files: %s', ME.message);
end

end

% Alternative vectorized approach for better performance
function InitialDistribution = processDistanceFreeFilesVectorized(myFolder)
% Vectorized version for better performance with large datasets

if nargin < 1
    myFolder = pwd;
end

if ~isfolder(myFolder)
    error('Folder does not exist: %s', myFolder);
end

% Get file list
filePattern = fullfile(myFolder, 'DistanceFree*.txt');
theFiles = dir(filePattern);
numFiles = length(theFiles);

if numFiles == 0
    warning('No files found matching pattern: %s', filePattern);
    InitialDistribution = [];
    return;
end

% Get full file paths
fileNames = fullfile({theFiles.folder}, {theFiles.name});

% Load first file to determine dimensions
firstData = load(fileNames{1});
dataSize = size(firstData);

% Pre-allocate 3D array for all data
allData = zeros([dataSize, numFiles]);

% Load all files
fprintf('Loading %d files...\n', numFiles);
parfor k = 1:numFiles  % Use parallel processing if available
    allData(:,:,k) = load(fileNames{k});
end

% Compute mean across third dimension
InitialDistribution = mean(allData, 3);
fprintf('Successfully processed %d files\n', numFiles);

end