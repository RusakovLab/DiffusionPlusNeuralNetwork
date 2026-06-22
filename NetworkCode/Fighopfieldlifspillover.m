%% Generalized Spiking Hopfield Network with Quality of Recall vs Noise
% This implementation plots Quality of Recall versus noise_level for multiple trials
% Original Author: L Savtchenko
% Modified: Nov 6, 2025
% MODIFICATION: Random patterns + Paper figure generation
%
% =====================================================================
% BIOPHYSICAL REFORMULATION (this version):
% ---------------------------------------------------------------------
% In the original implementation, M controlled the Hebbian fan-in
% (number of wired neighbours per neuron). Here, M is reinterpreted
% as a NEUROTRANSMITTER SPILLOVER COEFFICIENT acting during retrieval.
%
% Two postsynaptic pathways coexist, both with FIXED anatomy:
%
%   (1) Wired core W_core: each neuron has Hebbian connections to its
%       BasicSizeofNetwork (=3) nearest spatial neighbours. This is the
%       intra-cleft synaptic component. Fan-in is the same for every M.
%
%   (2) Spillover matrix W_spill: full Hebbian (Hopfield) matrix
%       multiplied element-wise by a row-normalised Gaussian footprint
%       of spatial width SigmaM. Biophysical reading: transmitter that
%       escapes the cleft activates EXTRASYNAPTIC receptors at
%       functionally co-tuned neighbours within an extracellular
%       diffusion radius set by SigmaM (a property of tissue geometry
%       and astrocytic uptake).
%
% Retrieval-phase synaptic drive at neuron i:
%       I_syn(i) = (1 - alpha) * (W_core  * r)_i
%                +     alpha   * (W_spill * r)_i
% with the spillover amplitude saturating in M:
%       alpha(M) = alpha_max * (1 - exp(-M / M_tau)).
%
% At low M the wired core dominates -> sparse local drive -> poor
% completion. At high M the extrasynaptic content-addressable drive
% dominates -> full Hopfield-style recall. M sets the contribution
% of spillover to retrieval, reproducing the monotone Q(M) trend of
% the original sparsity-driven model.
% =====================================================================

clear; clc; close all;
rng('shuffle'); % Random seed for trial variability
pattern_seed = 2; % Fixed seed for pattern generation only

%% ============== PAPER FIGURE CONTROL ==============
GENERATE_PAPER_FIGURES = true;  % SET TO true TO GENERATE 7-PATTERN FIGURES
%% ==================================================

%% PARAMETERS - All numerical parameters consolidated at the top
% ----------------------------------------------------------------------
% EXPERIMENT PARAMETERS
% ----------------------------------------------------------------------
noise_levels = 0:0.10:1;   % From 0% to 100% noise in 10% increments
num_trials = 100;          % Number of trials per noise level
ShowPlot_ = true;         % Show individual trial plots if true

% ----------------------------------------------------------------------
% NETWORK ARCHITECTURE PARAMETERS
% ----------------------------------------------------------------------
network_size = 400;        % Total number of neurons (for a 20x20 grid)
n_patterns = 3;            % Number of patterns to store
BasicSizeofNetwork = 3;    % Basic size of network without spillover (wired core)
M_list = [1, 10, 25, 50, 100, 200];  % M values for the final plot
SigmaM = 5;                % Spatial width (pixels on the grid) for both wiring prob and spillover
weight_scale = 25 * (100/network_size)^0.5;  % Scale weights based on network size

% ----------------------------------------------------------------------
% PAPER FIGURE M INDICES
% ----------------------------------------------------------------------
% Pattern 3: Recall with M (for Noise=0% figure)
M_indices_noise0 = [1];    % Corresponds to M=[1]

% Patterns 5-7: Recalls with M (for Noise=20% figure)  
M_indices_noise20 = [1, 2, 6]; % Corresponds to M=[1, 10, 200]
% ----------------------------------------------------------------------

% Spillover model parameters (M now acts as a spillover coefficient,
% not a sparsity / fan-in parameter)
spillover_enable = true;   % Turn the spillover mechanism on/off
alpha_max        = 1.0;    % Asymptotic spillover amplitude as M -> infinity
M_tau            = 30;     % M-scale at which alpha saturates (alpha(M_tau) = alpha_max*(1 - 1/e))

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

%% Run analysis across M values and noise levels
fprintf('\n========================================\n');
fprintf('Running Hopfield Network Analysis\n');
fprintf('M values: %s\n', mat2str(M_list));
fprintf('Noise levels: 0%% to 100%% (10%% steps)\n');
fprintf('Trials per condition: %d\n', num_trials);
fprintf('========================================\n\n');

% Storage for results across all M values
all_mean_quality = zeros(length(M_list), length(noise_levels));
all_std_quality = zeros(length(M_list), length(noise_levels));

% Storage for paper figures
if GENERATE_PAPER_FIGURES
    paper_data = struct();
    paper_data.stored_pattern = [];
    paper_data.noise0_input = [];
    paper_data.noise0_recalls = cell(length(M_list), 1);
    paper_data.noise20_input = [];
    paper_data.noise20_recalls = cell(length(M_list), 1);
    paper_data.qualities_noise0 = zeros(length(M_list), 1);
    paper_data.qualities_noise20 = zeros(length(M_list), 1);
    paper_data.std_noise0 = zeros(length(M_list), 1);
    paper_data.std_noise20 = zeros(length(M_list), 1);
end

for m_idx = 1:length(M_list)
    M = M_list(m_idx);
    fprintf('\n--- Processing M = %d ---\n', M);
    
    [mean_quality, std_quality, fig_data] = noise_analysis(network_size, noise_levels, num_trials, ...
        M, SigmaM, BasicSizeofNetwork, ShowPlot_, ...
        n_patterns, dt, T, TimeExtraStimulus_, weight_scale, ...
        tau_m, v_rest, v_thresh, v_reset, ref_period, ...
        stim_strength, noise_amp, spillover_enable, alpha_max, M_tau, ...
        GENERATE_PAPER_FIGURES, pattern_seed);
    
    all_mean_quality(m_idx, :) = mean_quality;
    all_std_quality(m_idx, :) = std_quality;
    
    % Store paper figure data
    if GENERATE_PAPER_FIGURES && ~isempty(fig_data)
        if m_idx == 1
            % Store patterns (same for all M)
            paper_data.stored_pattern = fig_data.stored_pattern;
            paper_data.noise0_input = fig_data.noise0_input;
            paper_data.noise20_input = fig_data.noise20_input;
        end
        % Store recalls for each M
        paper_data.noise0_recalls{m_idx} = fig_data.noise0_recall;
        paper_data.noise20_recalls{m_idx} = fig_data.noise20_recall;
        
        % Store quality values
        noise0_idx = find(noise_levels == 0);
        noise20_idx = find(abs(noise_levels - 0.2) < 0.01);
        paper_data.qualities_noise0(m_idx) = mean_quality(noise0_idx);
        paper_data.qualities_noise20(m_idx) = mean_quality(noise20_idx);
        paper_data.std_noise0(m_idx) = std_quality(noise0_idx);
        paper_data.std_noise20(m_idx) = std_quality(noise20_idx);
    end
end

%% Generate paper figures
if GENERATE_PAPER_FIGURES
    generate_paper_figures(paper_data, M_list, network_size, M_indices_noise0, M_indices_noise20);
end

%% Generate final quality vs noise plot for all M values
generate_quality_plot(all_mean_quality, all_std_quality, noise_levels, M_list);

fprintf('\n========================================\n');
fprintf('Analysis Complete!\n');
fprintf('========================================\n\n');

%% ========================= FUNCTIONS =========================

function [mean_quality, std_quality, fig_data] = noise_analysis(N, noise_levels, num_trials, ...
    M, SigmaM, BasicSizeofNetwork, showPlot, ...
    n_patterns, dt, T, TimeExtraStimulus, weight_scale, ...
    tau_m, v_rest, v_thresh, v_reset, ref_period, ...
    stim_strength, noise_amp, spillover_enable, alpha_max, M_tau, ...
    generate_fig_data, pattern_seed)

results = zeros(length(noise_levels), num_trials);

% Storage for paper figure data
fig_data = [];
if generate_fig_data
    fig_data = struct();
    fig_data.stored_pattern = [];
    fig_data.noise0_input = [];
    fig_data.noise0_recall = [];
    fig_data.noise20_input = [];
    fig_data.noise20_recall = [];
end

for noise_idx = 1:length(noise_levels)
    noise_level = noise_levels(noise_idx);
    fprintf('  Noise level: %.0f%% ... ', noise_level*100);
    
    for trial = 1:num_trials
        [quality, trial_data] = spiking_hopfield_network(N, noise_level, showPlot, ...
            M, SigmaM, BasicSizeofNetwork, ...
            n_patterns, dt, T, TimeExtraStimulus, weight_scale, ...
            tau_m, v_rest, v_thresh, v_reset, ref_period, ...
            stim_strength, noise_amp, spillover_enable, alpha_max, M_tau, ...
            pattern_seed);
        
        results(noise_idx, trial) = quality;
        
        % Store data for paper figures (first trial only)
        if generate_fig_data && trial == 1
            if noise_level == 0
                fig_data.stored_pattern = trial_data.stored_pattern;
                fig_data.noise0_input = trial_data.noisy_input;
                fig_data.noise0_recall = trial_data.retrieved_pattern;
            elseif abs(noise_level - 0.2) < 0.01  % noise = 20%
                fig_data.noise20_input = trial_data.noisy_input;
                fig_data.noise20_recall = trial_data.retrieved_pattern;
            end
        end
    end
    
    fprintf('Mean Q = %.1f%%, Std = %.1f%%\n', ...
        mean(results(noise_idx, :))*100, std(results(noise_idx, :))*100);
end

mean_quality = mean(results, 2);
std_quality  = std(results, 0, 2);
end

function [quality, trial_data] = spiking_hopfield_network(N, noise_level, show_plots, ...
    M, SigmaM, BasicSizeofNetwork, ...
    n_patterns, dt, T, TimeExtraStimulus, weight_scale, ...
    tau_m, v_rest, v_thresh, v_reset, ref_period, ...
    stim_strength, noise_amp, spillover_enable, alpha_max, M_tau, ...
    pattern_seed)

% Ensure N is a perfect square (for grid geometry)
grid_size = sqrt(N);
if floor(grid_size) ~= grid_size
    grid_size = round(sqrt(N));
    N = grid_size^2;
end
steps    = T/dt;
img_size = [grid_size, grid_size];

% Create patterns with fixed seed (same patterns for all trials)
patterns = create_patterns(N, n_patterns, img_size, pattern_seed);

% --- Build the FIXED wired Hebbian core W_core (independent of M) ---
% Fan-in is fixed at BasicSizeofNetwork; W_core is the same for every M.
% Also returns W_spill: the spatially-gated full Hebbian matrix
% representing extrasynaptic spillover coupling (also M-independent).
[W_core, W_spill] = build_core_and_kernel(patterns, N, weight_scale, ...
                                          BasicSizeofNetwork, SigmaM, img_size);

% --- Spillover amplitude (saturating in M) ---
% alpha = 0  -> no transmitter escape  (memory is purely the wired core)
% alpha -> 1 -> all wired drive is laterally diffused
if spillover_enable
    alpha_M = alpha_max * (1 - exp(-double(M) / double(max(1, M_tau))));
    alpha_M = min(1.0, max(0, alpha_M));
else
    alpha_M = 0;
end

% Pattern retrieval with activity-dependent spillover
[retrieval_quality, retrieved_pattern, noisy_input] = test_pattern_retrieval( ...
    N, steps, dt, patterns, W_core, W_spill, alpha_M, ...
    v_rest, v_thresh, v_reset, ref_period, tau_m, ...
    stim_strength, noise_amp, img_size, noise_level, show_plots, TimeExtraStimulus);

quality = retrieval_quality;

% Return trial data for paper figures
trial_data = struct();
trial_data.stored_pattern = patterns(:, 1);  % Use first pattern
trial_data.noisy_input = noisy_input;
trial_data.retrieved_pattern = retrieved_pattern;
end

%% Patterns - Random noise patterns
function patterns = create_patterns(N, n_patterns, img_size, pattern_seed)
patterns = zeros(N, n_patterns);

% Save current RNG state
current_rng = rng;

% Pattern 1: Random noise pattern (with fixed seed)
rng(pattern_seed + 101);
img1 = rand(img_size) > 0.5;
patterns(:, 1) = reshape(img1, [], 1) * 2 - 1;

% Pattern 2: Random noise pattern
rng(pattern_seed + 202);
img2 = rand(img_size) > 0.5;
patterns(:, 2) = reshape(img2, [], 1) * 2 - 1;

% Pattern 3: Random noise pattern
rng(pattern_seed + 303);
img3 = rand(img_size) > 0.5;
patterns(:, 3) = reshape(img3, [], 1) * 2 - 1;

% Restore RNG state (so trials remain independent)
rng(current_rng);

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

%% Core: fixed wired Hebbian matrix (M-independent) + spillover matrix
% W_core:  each neuron has Hebbian connections to its BasicSizeofNetwork
%          nearest spatial neighbours. This is the static synaptic anatomy
%          (the wired, intra-cleft component). Identical for every M.
%
% W_spill: spatially-gated full Hebbian matrix
%              W_spill(i,j) = G(i,j) * (1/N) sum_p xi_i^p * xi_j^p
%          where G is a row-normalised Gaussian footprint of width SigmaM.
%          Biophysical reading: transmitter that escapes the cleft at site j
%          activates extrasynaptic receptors at neurons i within a diffusion
%          radius set by SigmaM (controlled by astrocytic uptake / tortuosity);
%          those receptors are functionally co-tuned because Hopfield memory
%          embeds pattern-correlation in pairwise interactions.
%          W_spill is also M-independent; the M-dependence enters via the
%          mixing amplitude alpha(M) at retrieval time.
function [W_core, W_spill] = build_core_and_kernel(patterns, N, weight_scale, ...
    BasicSizeofNetwork, SigmaM, img_size)

% Grid positions
[Xg, Yg] = meshgrid(1:img_size(2), 1:img_size(1));
positions = [reshape(Yg, [], 1), reshape(Xg, [], 1)];

% Pairwise distances
D = zeros(N, N);
for i = 1:N
    di = positions(i,:);
    for j = 1:N
        dj = positions(j,:);
        D(i,j) = sqrt(sum((di - dj).^2));
    end
end

% Gaussian spatial profile (no self-mass: D=0 on diagonal -> P=1, zeroed below)
P = exp(-(D.^2) / (2 * max(1e-6, SigmaM)^2));
P(1:N+1:end) = 0;   % no self-spillover

% Sort neighbours by spatial proximity (used to pick wired partners)
P_for_sort = P;
P_for_sort(1:N+1:end) = -inf;

% --- Wired core: Hebbian on top-Basic spatial neighbours -----------------
W_core = zeros(N, N);
total_connections = min(BasicSizeofNetwork, N - 1);

for i = 1:N
    [~, idx] = sort(P_for_sort(i,:), 'descend');
    neighbors = idx(1:total_connections);
    for p = 1:size(patterns, 2)
        W_core(i, neighbors) = W_core(i, neighbors) + ...
            (patterns(i, p) * patterns(neighbors, p)') / N;
    end
end
W_core = W_core * weight_scale;
W_core(1:N+1:end) = 0;

% --- Spillover: spatially-gated full Hebbian matrix ----------------------
% Step 1: full Hebbian (Hopfield) matrix, same scaling as wired core
W_full = (patterns * patterns') / N;
W_full(1:N+1:end) = 0;
W_full = W_full * weight_scale;

% Step 2: row-normalised Gaussian gate G — restricts spillover to the
% extracellular diffusion radius (~SigmaM)
row_sums = sum(P, 2);
G = P;
nz = row_sums > 0;
G(nz, :) = G(nz, :) ./ row_sums(nz);

% Step 3: gated Hebbian — only co-tuned neighbours within diffusion reach
W_spill = G .* W_full;

% Match overall scale of W_spill to W_core so alpha=0..1 is a meaningful
% interpolation of effective drive magnitude (rather than a quiet -> loud
% transition that confounds gain with selectivity)
mean_abs_core  = mean(abs(W_core(:)))  + 1e-12;
mean_abs_spill = mean(abs(W_spill(:))) + 1e-12;
W_spill = W_spill * (mean_abs_core / mean_abs_spill);

W_spill(1:N+1:end) = 0;
end

%% Retrieval with activity-dependent spillover
% At each time step, the wired Hebbian drive is mixed with the
% extrasynaptic spillover drive, with mixing amplitude alpha_M (the
% spillover coefficient set by M):
%
%   I_wired = W_core  * r        (sparse, local, intra-cleft)
%   I_spill = W_spill * r        (full Hebbian gated by Gaussian footprint;
%                                 extrasynaptic / extracellular-diffusion
%                                 coupling to functionally co-tuned
%                                 neighbours within radius ~SigmaM)
%   I_syn   = (1 - alpha_M) * I_wired + alpha_M * I_spill
%
% This is the biophysical reformulation: M no longer changes the
% synaptic anatomy (W_core is fixed), it changes how much of the
% retrieval drive comes from extrasynaptic spillover-mediated coupling.
function [quality, retrieved_pattern, noisy_pattern] = test_pattern_retrieval( ...
    N, steps, dt, patterns, W_core, W_spill, alpha_M, ...
    v_rest, v_thresh, v_reset, ref_period, tau_m, ...
    stim_strength, noise_amp, img_size, noise_level, show_plots, TimeExtraStimulus)

test_pattern_idx = 1;
noisy_pattern = patterns(:, test_pattern_idx);
flip_indices = rand(N, 1) < noise_level;
noisy_pattern(flip_indices) = -noisy_pattern(flip_indices);

% Effective postsynaptic matrix (precomputed once per trial)
W_eff = (1 - alpha_M) * W_core + alpha_M * W_spill;

v = v_rest * ones(N, 1);
ref_time = zeros(N, 1);
spikes_retrieval = zeros(N, steps);

for t = 1:steps
    if t < TimeExtraStimulus/dt
        active_neurons = noisy_pattern > 0;
        external_input = zeros(N, 1);
        external_input(active_neurons) = stim_strength * 3;
    else
        external_input = zeros(N, 1);
    end

    if t > 1
        recent_window = max(1, t-20):t-1;
        if ~isempty(recent_window)
            recent_activity = sum(spikes_retrieval(:, recent_window), 2);
            r = double(recent_activity > 0);

            % Combined wired + spillover synaptic drive
            I_syn = W_eff * r * 5;   % factor 5 matches original gain
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

    if t == 100 && sum(sum(spikes_retrieval)) == 0
        active_neurons = noisy_pattern > 0;
        spikes_retrieval(active_neurons, t) = 1;
    end
end

% Final window average activity
final_window = (steps - round(5000)) : steps;
final_window = final_window(final_window > 0 & final_window <= steps);
final_activity = sum(spikes_retrieval(:, final_window), 2);

activity_threshold = max(0.01, median(final_activity) / 2);
retrieved_pattern = (final_activity > activity_threshold);

if sum(retrieved_pattern) == 0
    % Fallback: decode directly from W_eff acting on the noisy probe
    pattern_from_weights = W_eff * (noisy_pattern > 0);
    retrieved_pattern = (pattern_from_weights > median(pattern_from_weights));
end

binary_original = (patterns(:, test_pattern_idx) > 0);
quality = sum(retrieved_pattern .* binary_original) / sum(binary_original);
end

%% Generate paper figures (7 patterns in 2 figures)
function generate_paper_figures(paper_data, M_list, N, M_indices_noise0, M_indices_noise20)

grid_size = sqrt(N);

% Convert patterns to 2D images
stored_img = reshape((paper_data.stored_pattern > 0), grid_size, grid_size);
noise0_input_img = reshape((paper_data.noise0_input > 0), grid_size, grid_size);

% For noise=20% input: show original (red) and flipped bits (black)
stored_binary = (paper_data.stored_pattern > 0);
noisy_binary = (paper_data.noise20_input > 0);
flipped_bits = (stored_binary ~= noisy_binary); % Where bits were flipped

% Create visualization: 0=white, 1=red (original active), 2=black (flipped)
noise20_visualization = zeros(grid_size, grid_size);
noise20_visualization(reshape(stored_binary, grid_size, grid_size)) = 1; % Original pattern in red
noise20_visualization(reshape(flipped_bits, grid_size, grid_size)) = 2;  % Flipped bits in black

% Figure 1: Noise = 0%
figure('Position', [50, 100, 1400, 300], 'Color', 'w');

% Pattern 1: Stored pattern
subplot(1, 3, 1);
imagesc(stored_img); 
colormap(gca, [1 1 1; 1 0 0]); caxis([0 1]);
axis equal; axis tight; axis off;
title('Stored pattern', 'FontSize', 12, 'FontWeight', 'bold');

% Pattern 2: Noise-free probe
subplot(1, 3, 2);
imagesc(noise0_input_img);
colormap(gca, [1 1 1; 1 0 0]); caxis([0 1]);
axis equal; axis tight; axis off;
title('Noise-free Probe', 'FontSize', 12, 'FontWeight', 'bold');

% Pattern 3: Recall M (using M_indices_noise0)
m_idx = M_indices_noise0(1);  % Use first index from M_indices_noise0
recall_img = reshape((paper_data.noise0_recalls{m_idx} > 0), grid_size, grid_size);
subplot(1, 3, 3);
imagesc(recall_img);
colormap(gca, [1 1 1; 1 0 1]); caxis([0 1]);
axis equal; axis tight; axis off;
Q_mean = paper_data.qualities_noise0(m_idx) * 100;
Q_std = paper_data.std_noise0(m_idx) * 100;
title(sprintf('Recall Q: %.1f±%.1f%%\nM=%d', Q_mean, Q_std, M_list(m_idx)), ...
    'FontSize', 12, 'FontWeight', 'bold');

sgtitle('Pattern Retrieval: Noise = 0%', 'FontSize', 14, 'FontWeight', 'bold');

% Figure 2: Noise = 20%
figure('Position', [50, 450, 1800, 300], 'Color', 'w');

% Pattern 4: Noisy input (20% noise) - Red=original, Black=flipped
subplot(1, 4, 1);
imagesc(noise20_visualization);
colormap(gca, [1 1 1; 1 0 0; 0 0 0]); % White, Red, Black
caxis([0 2]);
axis equal; axis tight; axis off;
title('Noise: 20%', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Probe', 'FontSize', 10);

% Patterns 5-7: Recalls with M (using M_indices_noise20)
for i = 1:length(M_indices_noise20)
    m_idx = M_indices_noise20(i);
    subplot(1, 4, i+1);
    recall_img = reshape((paper_data.noise20_recalls{m_idx} > 0), grid_size, grid_size);
    imagesc(recall_img);
    colormap(gca, [1 1 1; 1 0 1]); caxis([0 1]);
    axis equal; axis tight; axis off;
    
    Q_mean = paper_data.qualities_noise20(m_idx) * 100;
    Q_std = paper_data.std_noise20(m_idx) * 100;
    title(sprintf('Q: %.1f±%.1f%%', Q_mean, Q_std), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(sprintf('M=%d', M_list(m_idx)), 'FontSize', 10, 'FontWeight', 'bold');
end

sgtitle('Pattern Retrieval: Noise = 20%', 'FontSize', 14, 'FontWeight', 'bold');

end

%% Generate quality vs noise plot
function generate_quality_plot(all_mean_quality, all_std_quality, noise_levels, M_list)

figure('Position', [100, 100, 800, 600], 'Color', 'w');
hold on;

colors = lines(length(M_list));

for m_idx = 1:length(M_list)
    errorbar(noise_levels*100, all_mean_quality(m_idx, :)*100, all_std_quality(m_idx, :)*100, ...
        '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', colors(m_idx, :), ...
        'Color', colors(m_idx, :), 'DisplayName', sprintf('M=%d', M_list(m_idx)));
end

grid on;
xlabel('Noise level, %', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Quality of recall Q, %', 'FontSize', 14, 'FontWeight', 'bold');
title('Quality of Recall vs Noise Level', 'FontSize', 16, 'FontWeight', 'bold');
xlim([-5, 105]);
ylim([0, 105]);
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
box on;

end