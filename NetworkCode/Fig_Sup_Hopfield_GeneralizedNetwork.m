function results = run_spiking_hopfield_v2(params)
% RUN_SPIKING_HOPFIELD_V2 - Hopfield retrieval with signed spillover (fixed)
% - Persistent state s ∈ {-1,+1}
% - Spillover is signed ring-convolution (K * s), kernel normalized
% - Safe default load P≈0.1N, classic 1/N Hebbian scaling
% - Primary metric: overlap m; report Q=(1+m)/2
% - Optional visualization of patterns, noisy inputs, and recalls
%


    if nargin < 1, params = default_params_v2(); else, params = merge_with_defaults_v2(params, default_params_v2()); end
    if isempty(params.P), params.P = max(1, floor(0.1*params.N)); end  % keep load below ~0.12
    rng(params.seed);

    fprintf('Generating patterns (bipolar)... ');
    patterns = generate_patterns_bipolar_v2(params.N, params.P);
    fprintf('done.\n');

    fprintf('Training Hebbian weights (1/N)... ');
    W = train_weights_v2(patterns, "hebb_1_over_N");
    fprintf('done.\n');

    fprintf('Precomputing normalized spillover kernels... ');
    kernel_bank = cell(numel(params.sigma_list), numel(params.M_list));
    for si = 1:numel(params.sigma_list)
        for mi = 1:numel(params.M_list)
            kernel_bank{si,mi} = make_spillover_kernel_weights_v2(params.N, params.sigma_list(si), params.M_list(mi));
        end
    end
    fprintf('done.\n\n');

    % Visualize stored patterns if requested
    if params.visualize_patterns
        visualize_stored_patterns(patterns, params);
    end
    
    % Save stored patterns data to CSV if requested
    if params.save_data_csv
        save_patterns_to_csv(patterns, params);
    end

    n_eta   = numel(params.etas);
    n_sigma = numel(params.sigma_list);
    n_M     = numel(params.M_list);

    m_all     = zeros(n_eta, n_sigma, n_M, params.trials);
    Q_all     = zeros(n_eta, n_sigma, n_M, params.trials);
    prec_all  = zeros(n_eta, n_sigma, n_M, params.trials);

    % Storage for visualization examples
    if params.visualize_patterns
        viz_data = struct();
        viz_data.examples = cell(min(3, n_eta), 1); % Store up to 3 examples
        viz_collected = 0;
    end

    total = n_eta*n_sigma*n_M*params.trials; prog = 0;
    fprintf('Running %d simulations...\n', total);

    for ei = 1:n_eta
        eta = params.etas(ei);
        for si = 1:n_sigma
            kw = kernel_bank{si,1}; %#ok<NASGU> % (just to keep linter quiet)
            for mi = 1:n_M
                kw = kernel_bank{si,mi};  % normalized weights for distances 1..M
                for tr = 1:params.trials
                    idx    = randi(params.P);
                    target = patterns(idx,:).';
                    s0     = apply_noise_bipolar_v2(target, eta);

                    s_final = simulate_hopfield_with_spillover_v2(W, s0, target, kw, params);
                    

                    m  = mean(s_final .* target);
                    Q  = (1 + m)/2;
                    pr = precision_against_active_bits_v2(s_final, target);

                    m_all(ei,si,mi,tr)    = m;
                    Q_all(ei,si,mi,tr)    = Q;
                    prec_all(ei,si,mi,tr) = pr;

                    % Collect visualization examples
                    if params.visualize_patterns && si==1 && mi==1 && tr==1
                        if viz_collected < min(3, n_eta) && (ei==1 || ei==round(n_eta/2) || ei==n_eta)
                            viz_collected = viz_collected + 1;
                            viz_data.examples{viz_collected} = struct(...
                                'target', target, 'noisy', s0, 'recalled', s_final, ...
                                'eta', eta, 'overlap', m);
                        end
                    end

                    prog = prog + 1;
                    if mod(prog, max(1,round(total/50)))==0, fprintf('.'); end
                end
                fprintf('  η=%.2f σ=%g M=%d | m≈%.3f Q≈%.3f\n', eta, params.sigma_list(si), params.M_list(mi), ...
                        mean(m_all(ei,si,mi,:),'all'), mean(Q_all(ei,si,mi,:),'all'));
            end
        end
    end
    fprintf('\nDone.\n');

    results = struct();
    results.m_mean    = mean(m_all, 4);
    results.m_std     = std(m_all, 0, 4);
    results.Q_mean    = mean(Q_all, 4);
    results.Q_std     = std(Q_all, 0, 4);
    results.prec_mean = mean(prec_all, 4);
    results.prec_std  = std(prec_all, 0, 4);
    results.params    = params;
    results.eta_values= params.etas;

    % Visualize recall examples if requested
    if params.visualize_patterns
        visualize_recall_examples(viz_data, params);
    end
    
    % Save recall examples data to CSV if requested
    if params.save_data_csv && params.visualize_patterns
        save_recall_examples_to_csv(viz_data, params);
    end

    if ~isempty(params.save_dir)
        if ~exist(params.save_dir,'dir'), mkdir(params.save_dir); end
        save(fullfile(params.save_dir,'results_v2.mat'), 'results');
        fprintf('Saved results to %s/results_v2.mat\n', params.save_dir);
    end

    plot_results_v2(results);
end

% ========================= Parameters & Utilities =========================

function p = default_params_v2()
    p = struct();
    p.N = 400;
    p.P = [3];                 % if empty, set to floor(0.1*N) at runtime

    p.etas   = 0:0.1:1;
    p.trials = 100;

    p.Tmax        = 200;
    p.halt_window = 5;
    p.beta        = 0.08;     % small cue clamp
    p.T_clamp     = 8;

    % Spillover field params (now SIGNED, normalized)
    p.sigma_list = [0, 2, 5, 10];
    p.M_list     = [200];
    p.alpha      = 4;      % small so it modulates, not dominates
   
    % Visualization options
    p.visualize_patterns = false;  % Set to true to enable pattern visualization
    p.n_patterns_show = 3;         % Number of patterns to show in visualization
    p.save_data_csv = false;       % Set to true to save plot data as CSV files

    p.seed     = 2111;
    p.save_dir = "results_v2";
end

function merged = merge_with_defaults_v2(user_params, defaults)
    merged = defaults;
    f = fieldnames(user_params);
    for i=1:numel(f)
        merged.(f{i}) = user_params.(f{i});
    end
end

function patterns = generate_patterns_bipolar_v2(N, P)
    assert(P <= min(0.8*N, 200), 'Too many patterns for this N');
    X = rand(P, N) > 0.5;
    patterns = 2*X - 1; % {-1,+1}
end

function W = train_weights_v2(patterns, ~)
    % classic Hopfield 1/N scaling
    [P, N] = size(patterns);
    W = (patterns' * patterns) / N;
    W = W - diag(diag(W));
    % no extra 1/sqrt(P) factor (keeps Hopfield signal strong)
end

function s_noisy = apply_noise_bipolar_v2(s, eta)
    N = numel(s); k = round(eta*N);
    idx = randperm(N, k);
    s_noisy = s; s_noisy(idx) = -s_noisy(idx);
end

function kw = make_spillover_kernel_weights_v2(N, sigma, M)
    %#ok<INUSD>
    
    if sigma <= 0 || M <= 0
        kw = zeros(1,0); return;
    end
    
    d = (1:M);
    kw = exp(-(d.^2)/(2*sigma^2));

    % Normalize so that total symmetric mass sum_d 2*kw(d) = 1
    s = 2*sum(kw);
    if s > 0, kw = kw / s; end
end




function s = simulate_hopfield_with_spillover_v2(W, s0, ~, kw, params)

    s = s0;
    N = numel(s);
    M = numel(kw);
    alpha = params.alpha;
    beta  = params.beta;
    
    % --- auto-calibrate alpha so spillover std is r * std(W*s) -------------
    r_target = alpha;                 % interpret params.alpha as ratio
    alpha = calibrate_alpha(W, kw, r_target, N);
    % -----------------------------------------------------------------------

    recent = zeros(params.halt_window, N);

    for t = 1:params.Tmax
        % Signed spillover via ring convolution over s (±1), kernel normalized
        if M > 0
            spill = zeros(N,1);
            for d = 1:M
                if kw(d) == 0, continue; end
                spill = spill + kw(d) * (circshift(s, d) + circshift(s, -d));
            end
            spill = alpha * spill;
        else
            spill = 0;
        end

        beta_t = (t <= params.T_clamp) * beta;

        h = W*s + spill + beta_t*s0;

        s_new = sign(h);
        z = (s_new==0); s_new(z) = s(z);

        recent(mod(t-1, params.halt_window)+1, :) = s_new.';
        if t >= params.halt_window
            if all(std(recent, 0, 1) == 0)
                s = s_new; break;
            end
        end
        s = s_new;
    end
end

function alpha_eff = calibrate_alpha(W, kw, r_target, N)
    % Estimate std of Ws and K*s on random states, set alpha so that
    % std(alpha * K*s) = r_target * std(W*s).
    if r_target <= 0
        alpha_eff = 0; return;
    end
    M = numel(kw);
    nSamples = 5;  % enough for a stable estimate, cheap
    stdWs = zeros(nSamples,1);
    stdKs = zeros(nSamples,1);
    for k = 1:nSamples
        s = sign(randn(N,1));
        hs = W*s;
        spill = zeros(N,1);
        for d = 1:M
            if kw(d)==0, continue; end
            spill = spill + kw(d) * (circshift(s,d) + circshift(s,-d));
        end
        stdWs(k) = std(hs);
        stdKs(k) = std(spill);
    end
    sW = mean(stdWs); sK = max(mean(stdKs), eps);
    alpha_eff = r_target * (sW / sK);
end

% ============================= Metrics & Plots ============================

function pr = precision_against_active_bits_v2(s_retrieved, s_target)
    R = (s_retrieved > 0);
    O = (s_target    > 0);
    TP = sum(R & O);
    P  = sum(O);
    if P == 0
        pr = double(sum(R)==0);
    else
        pr = TP / P;
    end
end

% ======================= Visualization Functions ==========================

function visualize_stored_patterns(patterns, params)
    [P, N] = size(patterns);
    n_show = min(params.n_patterns_show, P);
    
    % Determine grid size for visualization
    grid_size = ceil(sqrt(N));
    
    % Handle case where N is not a perfect square
    pad_size = grid_size^2 - N;
    
    figure('Position', [50 50 1200 300]);
    for i = 1:n_show
        subplot(1, n_show, i);
        % Pad pattern if necessary to make it square
        if pad_size > 0
            pattern_padded = [patterns(i,:), zeros(1, pad_size)];
        else
            pattern_padded = patterns(i,:);
        end
        % Reshape pattern into 2D grid for visualization
        pattern_2d = reshape(pattern_padded, grid_size, grid_size);
        imagesc(pattern_2d); colormap(gray); axis square; axis off;
        title(sprintf('Pattern %d', i));
        colorbar;
    end
    sgtitle(sprintf('Stored Patterns (N=%d, P=%d)', N, P));
    
    % Save figure
    if ~isempty(params.save_dir)
        if ~exist(params.save_dir,'dir'), mkdir(params.save_dir); end
        saveas(gcf, fullfile(params.save_dir, 'stored_patterns.png'));
        fprintf('Saved stored patterns visualization to %s/stored_patterns.png\n', params.save_dir);
    end
end

function visualize_recall_examples(viz_data, params)
    n_examples = length(viz_data.examples);
    if n_examples == 0, return; end
    
    N = length(viz_data.examples{1}.target);
    grid_size = ceil(sqrt(N));
    
    % Handle case where N is not a perfect square
    pad_size = grid_size^2 - N;
    
    figure('Position', [80 80 1400 400*n_examples]);
    
    for ex = 1:n_examples
        example = viz_data.examples{ex};
        
        % Pad patterns if necessary
        if pad_size > 0
            target_padded = [example.target(:); zeros(pad_size, 1)];
            noisy_padded = [example.noisy(:); zeros(pad_size, 1)];
            recalled_padded = [example.recalled(:); zeros(pad_size, 1)];
        else
            target_padded = example.target(:);
            noisy_padded = example.noisy(:);
            recalled_padded = example.recalled(:);
        end
        
        % Target pattern
        subplot(n_examples, 3, (ex-1)*3 + 1);
        target_2d = reshape(target_padded, grid_size, grid_size);
        imagesc(target_2d); colormap(gca, [1 0 0; 1 1 1]);  % 0 = w, 1 = b
    caxis([0 1]); axis square; axis off;
        title(sprintf('Original Pattern'));
        %colorbar;
        
        % Noisy pattern
        subplot(n_examples, 3, (ex-1)*3 + 2);
        noisy_2d = reshape(noisy_padded, grid_size, grid_size);
        NewNoiseyPatter=target_2d+noisy_2d;
        imagesc(NewNoiseyPatter);  axis square; axis off;
        colormap(gca, [1 0 0; 0 0 0; 1 1 1]);
caxis([-2 2]);       % [-2,0,2]
        
         title(sprintf('New noisy pattern (Noise=%.0f)', 100*example.eta));
        %colorbar;
        
        % Recalled pattern
        subplot(n_examples, 3, (ex-1)*3 + 3);
        recalled_2d = reshape(recalled_padded, grid_size, grid_size);
        imagesc(recalled_2d);  axis square; axis off;
        colormap(gca, [1 0 1; 1 1 1]);  % 0 = w, 1 = b
    caxis([0 1]);
        title(sprintf('Retrieved Pattern', example.overlap));
        %colorbar;
    end
    
    sgtitle('Pattern Recall Examples: Target → Noisy → Recalled');
    
    % Save figure
    if ~isempty(params.save_dir)
        if ~exist(params.save_dir,'dir'), mkdir(params.save_dir); end
        saveas(gcf, fullfile(params.save_dir, 'recall_examples.png'));
        fprintf('Saved recall examples visualization to %s/recall_examples.png\n', params.save_dir);
    end
end

function plot_results_v2(results)
    p = results.params; etas = results.eta_values;
    nM = numel(p.M_list); nS = numel(p.sigma_list); cols = lines(nS);

    % Q = (1+m)/2
    figure('Position',[80 80 1200 780]);
    for mi = 1:nM
        subplot(2, ceil(nM/2), mi); hold on;
        for si = 1:nS
            Qm = results.Q_mean(:,si,mi); Qs = results.Q_std(:,si,mi);
            errorbar(etas, Qm, Qs, 'LineWidth',2, 'Marker','o', 'Color', cols(si,:));
        end
        plot(etas, 0.5*ones(size(etas)), 'k--', 'LineWidth',1);
        xlabel('\eta (flip fraction)'); ylabel('Q = (1+m)/2');
        title(sprintf('M=%d, \\alpha=%.3f', p.M_list(mi), p.alpha));
        grid on; ylim([0 1]); xlim([0 1]);
        if mi==1
            leg = arrayfun(@(x) sprintf('\\sigma=%g', x), p.sigma_list, 'uni',0);
            leg{end+1}='Random baseline';
            legend(leg,'Location','best');
        end
    end
    sgtitle(sprintf('Hopfield retrieval with signed spillover: N=%d, P=%d', p.N, p.P));

    % Save figure
    if ~isempty(p.save_dir)
        if ~exist(p.save_dir,'dir'), mkdir(p.save_dir); end
        saveas(gcf, fullfile(p.save_dir,'Q_vs_noise_v2.png'));
    end

    % Overlap m
    figure('Position',[110 110 1200 400]);
    for mi = 1:nM
        subplot(1, nM, mi); hold on;
        for si = 1:nS
            mm = results.m_mean(:,si,mi); ms = results.m_std(:,si,mi);
            errorbar(etas, mm, ms, 'LineWidth',2, 'Marker','s', 'Color', cols(si,:));
        end
        plot(etas, zeros(size(etas)), 'k--', 'LineWidth',1);
        xlabel('\eta'); ylabel('Overlap m');
        title(sprintf('M=%d', p.M_list(mi)));
        grid on; ylim([-0.1 1]); xlim([0 1]);
        if mi==1
            leg = arrayfun(@(x) sprintf('\\sigma=%g', x), p.sigma_list, 'uni',0);
            leg{end+1}='m=0 baseline';
            legend(leg,'Location','best');
        end
    end

    % Save figure
    if ~isempty(p.save_dir)
        saveas(gcf, fullfile(p.save_dir,'m_vs_noise_v2.png'));
    end

    % Precision (secondary)
    figure('Position',[140 140 1200 400]);
    for mi = 1:nM
        subplot(1, nM, mi); hold on;
        for si = 1:nS
            pm = results.prec_mean(:,si,mi); ps = results.prec_std(:,si,mi);
            errorbar(etas, pm, ps, 'LineWidth',2, 'Marker','^', 'Color', cols(si,:));
        end
        plot(etas, 0.5*ones(size(etas)), 'k--', 'LineWidth',1);
        xlabel('\eta'); ylabel('Precision (active bits)');
        title(sprintf('M=%d', p.M_list(mi)));
        grid on; ylim([0 1]); xlim([0 1]);
        if mi==1
            leg = arrayfun(@(x) sprintf('\\sigma=%g', x), p.sigma_list, 'uni',0);
            leg{end+1}='Random baseline';
            legend(leg,'Location','best');
        end
    end

    % Save figure
    if ~isempty(p.save_dir)
        saveas(gcf, fullfile(p.save_dir,'precision_vs_noise_v2.png'));
        fprintf('Saved plots to %s/\n', p.save_dir);
    end
    
    % Save plot data to CSV if requested
    if p.save_data_csv
        save_plot_data_to_csv(results, p);
    end
end

% ======================= CSV Data Export Functions ========================

function save_patterns_to_csv(patterns, params)
    if ~exist(params.save_dir,'dir'), mkdir(params.save_dir); end
    
    [P, N] = size(patterns);
    
    % Save patterns matrix
    csvwrite(fullfile(params.save_dir, 'stored_patterns.csv'), patterns);
    
    % Save metadata
    fid = fopen(fullfile(params.save_dir, 'stored_patterns_metadata.csv'), 'w');
    fprintf(fid, 'Parameter,Value\n');
    fprintf(fid, 'Number_of_Patterns,%d\n', P);
    fprintf(fid, 'Pattern_Size,%d\n', N);
    fprintf(fid, 'Grid_Size,%d\n', ceil(sqrt(N)));
    fclose(fid);
    
    fprintf('Saved stored patterns data to %s/stored_patterns.csv\n', params.save_dir);
end

function save_recall_examples_to_csv(viz_data, params)
    if ~exist(params.save_dir,'dir'), mkdir(params.save_dir); end
    
    n_examples = length(viz_data.examples);
    if n_examples == 0, return; end
    
    N = length(viz_data.examples{1}.target);
    grid_size = ceil(sqrt(N));
    pad_size = grid_size^2 - N;
    
    for ex = 1:n_examples
        example = viz_data.examples{ex};
        
        % Pad patterns if necessary and reshape all at once
        if pad_size > 0
            patterns_combined = [example.target(:); zeros(pad_size,1); ...
                                example.noisy(:); zeros(pad_size,1); ...
                                example.recalled(:); zeros(pad_size,1)];
        else
            patterns_combined = [example.target(:); example.noisy(:); example.recalled(:)];
        end
        
        % Reshape into 3 columns of grid_size x grid_size matrices
        data = reshape(patterns_combined, grid_size, grid_size*3);
        
        % Save to CSV
        filename = sprintf('recall_example_%d_eta_%.2f.csv', ex, example.eta);
        csvwrite(fullfile(params.save_dir, filename), data);
        
        % Save metadata for this example
        meta_filename = sprintf('recall_example_%d_metadata.csv', ex);
        fid = fopen(fullfile(params.save_dir, meta_filename), 'w');
        fprintf(fid, 'Parameter,Value\n');
        fprintf(fid, 'Example_Number,%d\n', ex);
        fprintf(fid, 'Eta,%.4f\n', example.eta);
        fprintf(fid, 'Overlap,%.4f\n', example.overlap);
        fprintf(fid, 'Pattern_Size,%d\n', N);
        fprintf(fid, 'Grid_Size,%d\n', grid_size);
        fprintf(fid, 'Padded_Size,%d\n', grid_size^2);
        fprintf(fid, 'Data_Format,Target|Noisy|Recalled (3 blocks of %dx%d)\n', grid_size, grid_size);
        fclose(fid);
    end
    
    fprintf('Saved recall examples data to %s/recall_example_*.csv\n', params.save_dir);
end

function save_plot_data_to_csv(results, params)
    if ~exist(params.save_dir,'dir'), mkdir(params.save_dir); end
    
    etas = results.eta_values;
    n_eta = length(etas);
    n_sigma = numel(params.sigma_list);
    n_M = numel(params.M_list);
    
    % Save Q vs noise data
    for mi = 1:n_M
        % Prepare data matrix
        data = zeros(n_eta, 1 + 2*n_sigma);
        data(:,1) = etas(:);
        
        for si = 1:n_sigma
            data(:, 1 + 2*(si-1) + 1) = results.Q_mean(:,si,mi);
            data(:, 1 + 2*(si-1) + 2) = results.Q_std(:,si,mi);
        end
        
        % Create header
        header = 'eta';
        for si = 1:n_sigma
            header = [header, sprintf(',Q_mean_sigma_%.4f,Q_std_sigma_%.4f', ...
                      params.sigma_list(si), params.sigma_list(si))];
        end
        
        % Save to file
        filename = sprintf('Q_vs_noise_M_%d.csv', params.M_list(mi));
        fid = fopen(fullfile(params.save_dir, filename), 'w');
        fprintf(fid, '%s\n', header);
        fclose(fid);
        dlmwrite(fullfile(params.save_dir, filename), data, '-append', 'precision', 6);
    end
    
    % Save overlap m vs noise data
    for mi = 1:n_M
        data = zeros(n_eta, 1 + 2*n_sigma);
        data(:,1) = etas(:);
        
        for si = 1:n_sigma
            data(:, 1 + 2*(si-1) + 1) = results.m_mean(:,si,mi);
            data(:, 1 + 2*(si-1) + 2) = results.m_std(:,si,mi);
        end
        
        header = 'eta';
        for si = 1:n_sigma
            header = [header, sprintf(',m_mean_sigma_%.4f,m_std_sigma_%.4f', ...
                      params.sigma_list(si), params.sigma_list(si))];
        end
        
        filename = sprintf('m_vs_noise_M_%d.csv', params.M_list(mi));
        fid = fopen(fullfile(params.save_dir, filename), 'w');
        fprintf(fid, '%s\n', header);
        fclose(fid);
        dlmwrite(fullfile(params.save_dir, filename), data, '-append', 'precision', 6);
    end
    
    % Save precision vs noise data
    for mi = 1:n_M
        data = zeros(n_eta, 1 + 2*n_sigma);
        data(:,1) = etas(:);
        
        for si = 1:n_sigma
            data(:, 1 + 2*(si-1) + 1) = results.prec_mean(:,si,mi);
            data(:, 1 + 2*(si-1) + 2) = results.prec_std(:,si,mi);
        end
        
        header = 'eta';
        for si = 1:n_sigma
            header = [header, sprintf(',precision_mean_sigma_%.4f,precision_std_sigma_%.4f', ...
                      params.sigma_list(si), params.sigma_list(si))];
        end
        
        filename = sprintf('precision_vs_noise_M_%d.csv', params.M_list(mi));
        fid = fopen(fullfile(params.save_dir, filename), 'w');
        fprintf(fid, '%s\n', header);
        fclose(fid);
        dlmwrite(fullfile(params.save_dir, filename), data, '-append', 'precision', 6);
    end
    
    fprintf('Saved all plot data to CSV files in %s/\n', params.save_dir);
end