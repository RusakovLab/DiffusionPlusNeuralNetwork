%% Generalized Spiking Hopfield Network with Quality of Recall vs Noise
% This implementation plots Quality of Recall versus noise_level for multiple trials
% Original Author: L Savtchenko
% Modified: Nov 6, 2025

clear; clc; close all;
rng(2); % For reproducibility

%% PARAMETERS - All numerical parameters consolidated at the top
% ----------------------------------------------------------------------
% EXPERIMENT PARAMETERS
% ----------------------------------------------------------------------
noise_levels = 0:0.10:1;   % From 0% to 100% noise in 10% increments
num_trials = 5;            % Number of trials per noise level
ShowPlot_ = false;         % Show all plots if true

% ----------------------------------------------------------------------
% NETWORK ARCHITECTURE PARAMETERS
% ----------------------------------------------------------------------
network_size = 400;        % Total number of neurons (for a 20x20 grid)
n_patterns = 3;            % Number of patterns to store
BasicSizeofNetwork = 3;    % Basic size of network without spillover (wired core)
M = 10;                   % Additional connections beyond BasicSizeofNetwork (drives spillover)
SigmaM = 5;                % Spatial width (pixels on the grid) for both wiring prob and spillover
weight_scale = 25 * (100/network_size)^0.5;  % Scale weights based on network size

% Spillover control (kept conservative to preserve prior results)
spillover_enable = true;   % Turn the spillover mechanism on/off
spillover_cap = 0.30;      % Hard cap on spillover fraction to protect results
M_ref_for_spill = max(1, 0.5*sqrt(network_size)^2); % default ~N/2; sets s(M) growth scale

% ----------------------------------------------------------------------
% SIMULATION PARAMETERS
% ----------------------------------------------------------------------
dt = 0.1;                  % Time step (ms)
T = 10000;                 % Simulation time (ms)
TimeExtraStimulus_ = 200;  % Duration of extra stimulation (ms)
stim_strength = 30;        % Stimulation strength (mV)
noise_amp = 3;             % Background noise amplitude (mV)

% ----------------------------------------------------------------------
% NEURON PARAMETERS
% ----------------------------------------------------------------------
tau_m = 8;                 % Membrane time constant (ms) - faster response
v_rest = -70;              % Resting potential (mV)
v_thresh = -60;            % Threshold potential (mV) - lower threshold
v_reset = -75;             % Reset potential (mV)
ref_period = 2;            % Refractory period (ms) - shorter refractory

% ----------------------------------------------------------------------
% OUTPUT PARAMETERS
% ----------------------------------------------------------------------
filename = sprintf('M=%.1f_SigmaM=%.1f.csv', M, SigmaM);

%% Run analysis across noise levels and trials
[mean_quality, std_quality] = noise_analysis(network_size, noise_levels, num_trials, ...
    M, SigmaM, BasicSizeofNetwork, ShowPlot_, ...
    n_patterns, dt, T, TimeExtraStimulus_, weight_scale, ...
    tau_m, v_rest, v_thresh, v_reset, ref_period, ...
    stim_strength, noise_amp, spillover_enable, spillover_cap, M_ref_for_spill);

% Create a table with the data
data_table = table(noise_levels'*100, mean_quality*100, std_quality*100, ...
    'VariableNames', {'NoiseLevel', 'Quality', 'Error'});

% Write the table to the CSV file with the parametric filename
writetable(data_table, filename);

% Create the errorbar plot
figure('Position', [100, 100, 800, 600]);
errorbar(noise_levels*100, mean_quality*100, std_quality*100, '-o', ...
    'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
grid on;
xlabel('Noise Level (%)', 'FontSize', 14);
ylabel('Quality of Recall (%)', 'FontSize', 14);
title(['Local Connectivity Spiking Hopfield Network: Quality vs Noise (M=' num2str(M) ', \sigma_M=' num2str(SigmaM) ')'], 'FontSize', 16);
xlim([-5, 105]);
ylim([0, 100]);
set(gca, 'FontSize', 12);

% Display summary
fprintf('\nNoise Analysis Results with Local Connectivity (M=%d, SigmaM=%.1f):\n', M, SigmaM);
for i = 1:length(noise_levels)
    fprintf('Noise level %.0f%%: Mean quality = %.2f%%, Std = %.2f%%\n', ...
        noise_levels(i)*100, mean_quality(i)*100, std_quality(i)*100);
end

%% FUNCTIONS

function [mean_quality, std_quality] = noise_analysis(N, noise_levels, num_trials, ...
    M, SigmaM, BasicSizeofNetwork, showPlot, ...
    n_patterns, dt, T, TimeExtraStimulus, weight_scale, ...
    tau_m, v_rest, v_thresh, v_reset, ref_period, ...
    stim_strength, noise_amp, spillover_enable, spillover_cap, M_ref_for_spill)

results = zeros(length(noise_levels), num_trials);

for noise_idx = 1:length(noise_levels)
    noise_level = noise_levels(noise_idx);
    fprintf('\nTesting noise level: %.1f%%\n', noise_level*100);

    for trial = 1:num_trials
        fprintf('  Trial %d of %d: ', trial, num_trials);

        quality = spiking_hopfield_network(N, noise_level, showPlot, ...
            M, SigmaM, BasicSizeofNetwork, ...
            n_patterns, dt, T, TimeExtraStimulus, weight_scale, ...
            tau_m, v_rest, v_thresh, v_reset, ref_period, ...
            stim_strength, noise_amp, spillover_enable, spillover_cap, M_ref_for_spill);

        results(noise_idx, trial) = quality;
        fprintf('Quality = %.2f%%\n', quality*100);
    end
end

mean_quality = mean(results, 2);
std_quality  = std(results, 0, 2);
end

function quality = spiking_hopfield_network(N, noise_level, show_plots, ...
    M, SigmaM, BasicSizeofNetwork, ...
    n_patterns, dt, T, TimeExtraStimulus, weight_scale, ...
    tau_m, v_rest, v_thresh, v_reset, ref_period, ...
    stim_strength, noise_amp, spillover_enable, spillover_cap, M_ref_for_spill)

% Ensure N is a perfect square (for grid geometry)
grid_size = sqrt(N);
if floor(grid_size) ~= grid_size
    grid_size = round(sqrt(N));
    N = grid_size^2;
    if show_plots
        fprintf('Adjusting N to %d to make a perfect square for visualization\n', N);
    end
end
steps    = T/dt;
img_size = [grid_size, grid_size];

if show_plots
    fprintf('\nNetwork Parameters:\n');
    fprintf('  Number of neurons: %d (%d x %d grid)\n', N, grid_size, grid_size);
    fprintf('  Basic network size: %d\n', BasicSizeofNetwork);
    fprintf('  Additional connections (M): %d\n', M);
    fprintf('  Gaussian width (SigmaM): %.1f\n', SigmaM);
    fprintf('  Total requested connections/neuron: %d\n', BasicSizeofNetwork + M);
    fprintf('  Threshold: %d mV | Weight scale: %.2f | Stim: %d mV | Noise: %.1f%%\n', ...
        v_thresh, weight_scale, stim_strength, noise_level*100);
end

% Create patterns that scale with network size
patterns = create_patterns(N, n_patterns, img_size);

% Create synaptic weight matrix with realistic spillover mixing (see below)
W = create_local_weight_matrix_with_spillover( ...
        patterns, N, weight_scale, M, SigmaM, BasicSizeofNetwork, img_size, ...
        spillover_enable, spillover_cap, M_ref_for_spill);

% Optional diagnostics
if show_plots
    fprintf('Weight matrix analysis:\n');
    fprintf('  Min weight: %.2f | Max weight: %.2f | Mean |W|: %.2f | Density: %.2f%%\n', ...
        min(W(:)), max(W(:)), mean(abs(W(:))), 100*nnz(W)/N^2);

    if N <= 100
        figure; imagesc(W ~= 0); title('Network Connectivity'); colormap(gray);
        axis square; xlabel('Neuron Index'); ylabel('Neuron Index');
    end
    visualize_patterns(patterns, n_patterns, img_size);
end

% Pattern retrieval
retrieval_quality = test_pattern_retrieval(N, steps, dt, patterns, W, v_rest, v_thresh, v_reset,...
    ref_period, tau_m, stim_strength, noise_amp, img_size, noise_level, show_plots, TimeExtraStimulus);

quality = retrieval_quality;
end

%% Patterns
function patterns = create_patterns(N, n_patterns, img_size)
patterns = zeros(N, n_patterns);

% Pattern 1: Horizontal stripes
img1 = zeros(img_size);
stripe_width = max(1, floor(img_size(1)/5));
stripe_positions = round(linspace(1, img_size(1)-stripe_width, 3));
for i = 1:length(stripe_positions)
    pos = stripe_positions(i);
    img1(pos:pos+stripe_width-1, :) = 1;
end
patterns(:, 1) = reshape(img1, [], 1) * 2 - 1;

% Pattern 2: Vertical stripes
img2 = zeros(img_size);
stripe_positions = round(linspace(1, img_size(2)-stripe_width, 3));
for i = 1:length(stripe_positions)
    pos = stripe_positions(i);
    img2(:, pos:pos+stripe_width-1) = 1;
end
patterns(:, 2) = reshape(img2, [], 1) * 2 - 1;

% Pattern 3: Diagonal band
img3 = zeros(img_size);
diag_width = max(1, floor(img_size(1)/10));
for i = 1:img_size(1)
    center = i;
    for w = -diag_width:diag_width
        row = i; col = center + w;
        if col >= 1 && col <= img_size(2)
            img3(row, col) = 1;
        end
    end
end
patterns(:, 3) = reshape(img3, [], 1) * 2 - 1;

% Balance activation ~40%
for p = 1:n_patterns
    active_ratio = sum(patterns(:,p) > 0) / N;
    if active_ratio < 0.3 || active_ratio > 0.5
        target_count  = round(0.4 * N);
        current_count = sum(patterns(:,p) > 0);
        if current_count < target_count
            inactive = find(patterns(:,p) < 0);
            to_activate = inactive(randperm(length(inactive), min(target_count - current_count, length(inactive))));
            patterns(to_activate, p) = 1;
        else
            active = find(patterns(:,p) > 0);
            to_deactivate = active(randperm(length(active), min(current_count - target_count, length(active))));
            patterns(to_deactivate, p) = -1;
        end
    end
end
end

%% Core: local weights + spillover (biophysical, normalized to preserve results)
function W = create_local_weight_matrix_with_spillover(patterns, N, weight_scale, ...
    M, SigmaM, BasicSizeofNetwork, img_size, spillover_enable, spillover_cap, M_ref_for_spill)

% Grid positions
[Xg, Yg] = meshgrid(1:img_size(2), 1:img_size(1));
positions = [reshape(Yg, [], 1), reshape(Xg, [], 1)]; % [row, col]

% Pairwise distances (Euclidean)
D = zeros(N, N);
for i = 1:N
    di = positions(i,:);
    for j = 1:N
        dj = positions(j,:);
        D(i,j) = sqrt(sum((di - dj).^2));
    end
end
D(1:N+1:end) = inf; % no self

% Gaussian probability for wiring and diffusion
P = exp(-(D.^2) / (2 * max(1e-6, SigmaM)^2)); % spatial falloff
P(1:N+1:end) = 0;

% --- 1) Build DIRECT wired weights (Hebbian on top-M neighbors by P)
W_direct = zeros(N,N);
total_connections = min(BasicSizeofNetwork + M, N-1);

for i = 1:N
    [~, idx] = sort(P(i,:), 'descend');
    neighbors = idx(1:total_connections);
    for p = 1:size(patterns,2)
        % Hebbian only for selected neighbors
        W_direct(i, neighbors) = W_direct(i, neighbors) + (patterns(i, p) * patterns(neighbors, p)' )/N;
    end
end

% Scale to original magnitude baseline
W_direct = W_direct * weight_scale;

% --- 2) Build row-stochastic diffusion kernel K (spillover mixing)
% Normalize each row to sum 1 (if a row is all zeros, keep it zero)
row_sums = sum(P,2);
K = P;
nz = row_sums > 0;
K(nz, :) = K(nz, :) ./ row_sums(nz);

% --- 3) Spillover fraction s(M): saturating, capped, monotone in M
if spillover_enable
    s_raw = 1 - exp( -double(M) / double(max(1, M_ref_for_spill)) ); % 0..~1
    s = min(spillover_cap, max(0, s_raw));                            % conservative cap
else
    s = 0;
end

% --- 4) Mix direct synapses with diffused copy
% Biophysical read: transmitter released at wired synapses diffuses with kernel K
% so the effective postsynaptic influence is (1-s)*W_direct + s*(W_direct*K)
W = (1 - s) * W_direct + s * (W_direct * K);

% --- 5) Renormalize to preserve per-row gain and global |W| statistics
% Preserve mean |W| (global)
mean_abs_direct = mean(abs(W_direct(:)) + 1e-12);
mean_abs_W      = mean(abs(W(:)) + 1e-12);
W = W * (mean_abs_direct / mean_abs_W);

% Optional: soft row normalization to preserve total outgoing drive
row_sum_direct = sum(abs(W_direct), 2) + 1e-12;
row_sum_W      = sum(abs(W), 2) + 1e-12;
scale_rows     = row_sum_direct ./ row_sum_W;
W = diag(scale_rows) * W;

% Keep diagonal strictly zero
W(1:N+1:end) = 0;
end

%% Visualize patterns
function visualize_patterns(patterns, n_patterns, img_size)
figure;
for p = 1:n_patterns
    subplot(1, n_patterns, p);
    imagesc(reshape(patterns(:, p), img_size));
    title(['Pattern ' num2str(p)]);
    colormap(gray);
    axis equal; axis tight; axis off;
end
end

%% Retrieval
function quality = test_pattern_retrieval(N, steps, dt, patterns, W, v_rest, v_thresh, v_reset,...
    ref_period, tau_m, stim_strength, noise_amp, img_size, noise_level, show_plots, TimeExtraStimulus)

% Positions for optional maps
[Xg, Yg] = meshgrid(1:img_size(2), 1:img_size(1));
positions = [reshape(Yg, [], 1), reshape(Xg, [], 1)];

test_pattern_idx = 1;
noisy_pattern = patterns(:, test_pattern_idx);
flip_indices = rand(N, 1) < noise_level;
noisy_pattern(flip_indices) = -noisy_pattern(flip_indices);

v = v_rest * ones(N, 1);
ref_time = zeros(N, 1);
spikes_retrieval = zeros(N, steps);

for t = 1:steps
    % Stronger initial cue for first TimeExtraStimulus ms only
    if t < TimeExtraStimulus/dt
        active_neurons = noisy_pattern > 0;
        external_input = zeros(N, 1);
        external_input(active_neurons) = stim_strength * 3;
    else
        external_input = zeros(N, 1);
    end

    % Synaptic input (recent spikes → PSPs)
    if t > 1
        recent_window = max(1, t-20):t-1;
        if ~isempty(recent_window)
            recent_activity = sum(spikes_retrieval(:, recent_window), 2);
            I_syn = W * (recent_activity > 0) * 5;
        else
            I_syn = zeros(N, 1);
        end
    else
        I_syn = zeros(N, 1);
    end

    I_noise = noise_amp * randn(N, 1);

    for i = 1:N
        if ref_time(i) > 0
            ref_time(i) = ref_time(i) - dt;
            v(i) = v_reset;
        else
            dv = ((v_rest - v(i)) + I_syn(i) + I_noise(i) + external_input(i)) * dt / tau_m;
            v(i) = v(i) + dv;

            if v(i) >= v_thresh
                spikes_retrieval(i, t) = 1;
                v(i) = v_reset;
                ref_time(i) = ref_period;
            end
        end
    end

    % Safety kick if network is completely silent
    if t == 100 && sum(sum(spikes_retrieval)) == 0
        if show_plots, fprintf('FORCING RETRIEVAL ACTIVITY\n'); end
        active_neurons = noisy_pattern > 0;
        spikes_retrieval(active_neurons, t) = 1;
    end
end

% Final window average activity
final_window = (steps - round(5000)) : steps;
final_window = final_window(final_window > 0 & final_window <= steps);
final_activity = sum(spikes_retrieval(:, final_window), 2);

% Adaptive threshold
activity_threshold = max(0.01, median(final_activity) / 2);
retrieved_pattern = (final_activity > activity_threshold);

% Fallback if no activity
if sum(retrieved_pattern) == 0
    if show_plots
        fprintf('WARNING: No retrieval activity detected, using weight-based retrieval\n');
    end
    pattern_from_weights = W * (noisy_pattern > 0);
    retrieved_pattern = (pattern_from_weights > median(pattern_from_weights));
end

binary_original = (patterns(:, test_pattern_idx) > 0);
quality = sum(retrieved_pattern .* binary_original) / sum(binary_original);

% Optional visualization
if show_plots
    figure;
    subplot(1,3,1);
    imagesc(reshape((patterns(:, test_pattern_idx)>0), img_size));
    colormap(gca, [1 0 0; 1 1 1]);  caxis([0 1]);
    title('Original Pattern'); axis equal; axis tight; axis off;

    NewNoiseyPatter = noisy_pattern + patterns(:, test_pattern_idx);
    subplot(1,3,2);
    imagesc(reshape(NewNoiseyPatter, img_size));
    colormap(gca, [1 0 0; 0 0 0; 1 1 1]); caxis([-2 2]);
    title(['New noisy pattern (' num2str(noise_level*100) '% noise)']);
    axis equal; axis tight; axis off;

    subplot(1,3,3);
    imagesc(reshape(retrieved_pattern, img_size));
    colormap(gca, [1 0 1; 1 1 1]); caxis([0 1]);
    title('Retrieved Pattern'); axis equal; axis tight; axis off;
    xlabel(['Overlap: ' num2str(quality, '%.2f')]);

    if N <= 400
        center_neuron = round(N/2);
        figure;

        % Connectivity map
        conn_map = zeros(img_size);
        conn_map(positions(:,1) + (positions(:,2)-1)*img_size(1)) = W(center_neuron, :)' ~= 0;
        center_pos = positions(center_neuron, :);
        subplot(1,2,1);
        imagesc(conn_map); hold on;
        plot(center_pos(2), center_pos(1), 'r*', 'MarkerSize', 10);
        title(['Connections from Neuron ' num2str(center_neuron)]);
        axis equal; axis tight; colormap(gray);

        % Strength map
        subplot(1,2,2);
        strength_map = zeros(img_size);
        strength_map(positions(:,1) + (positions(:,2)-1)*img_size(1)) = W(center_neuron, :)';
        imagesc(strength_map); hold on;
        plot(center_pos(2), center_pos(1), 'r*', 'MarkerSize', 10);
        title('Connection Strengths'); axis equal; axis tight; colorbar;
    end
end
end
