% New LIF model of Hopfiled network for the spillover paper. The same
% result as for non spikingh network

function results = run_spiking_hopfield_LIF_v3(params)
% RUN_SPIKING_HOPFIELD_LIF_V3
% Hopfield memory retrieval using LIF neurons with signed spillover
%
% NOISE-DEPENDENT QUALITY ANALYSIS:
% - Initial cue is corrupted by noise level η ∈ [0,1] (fraction of flipped bits)
% - Retrieval quality Q(η) = (1+m(η))/2, where m(η) is noise-dependent overlap
% - Expected behavior: Q(0) ≈ 1 (perfect), Q(1) ≈ 0.5 (random)
%
% LIF dynamics:
% - Leaky Integrate-and-Fire neurons with recurrent excitation/inhibition
% - Persistent attractor recall through Hebbian weights
% - Optional local "spillover" field (ring convolution) for lateral interactions
%
% Example:
%   p = default_params_LIF();
%   results = run_spiking_hopfield_LIF_v3(p);

    if nargin < 1, params = default_params_LIF(); 
    else, params = merge_with_defaults(params, default_params_LIF()); end
    if isempty(params.P), params.P = max(1, floor(0.1*params.N)); end
    rng(params.seed);

    fprintf('=== NOISE-DEPENDENT HOPFIELD RETRIEVAL ===\n');
    fprintf('Network size N=%d, Patterns P=%d\n', params.N, params.P);
    fprintf('Noise levels η: [%.2f, %.2f] in %d steps\n\n', ...
        min(params.etas), max(params.etas), numel(params.etas));

    fprintf('Generating patterns... ');
    patterns = generate_patterns_bipolar(params.N, params.P);
    fprintf('done.\n');

    fprintf('Training weights (Hebb)... ');
    W = train_weights(patterns, params);
    fprintf('done. Weight scale: %.2f\n', mean(abs(W(:))));

    fprintf('Preparing spillover kernels... ');
    kernel_bank = cell(numel(params.sigma_list), numel(params.M_list));
    for si = 1:numel(params.sigma_list)
        for mi = 1:numel(params.M_list)
            kernel_bank{si,mi} = make_spillover_kernel(params.N, ...
                params.sigma_list(si), params.M_list(mi));
        end
    end
    fprintf('done.\n\n');

    n_eta   = numel(params.etas);
    n_sigma = numel(params.sigma_list);
    n_M     = numel(params.M_list);

    % Storage for noise-dependent metrics
    m_all     = zeros(n_eta, n_sigma, n_M, params.trials);
    Q_all     = zeros(n_eta, n_sigma, n_M, params.trials);
    prec_all  = zeros(n_eta, n_sigma, n_M, params.trials);
    
    % Track noise corruption and retrieval success
    initial_overlap_all = zeros(n_eta, n_sigma, n_M, params.trials);

    total = n_eta*n_sigma*n_M*params.trials; prog = 0;
    fprintf('Running %d simulations across noise levels...\n', total);

    for ei = 1:n_eta
        eta = params.etas(ei);  % Current noise level
        
        for si = 1:n_sigma
            for mi = 1:n_M
                kw = kernel_bank{si,mi};
                
                for tr = 1:params.trials
                    idx    = randi(params.P);
                    target = patterns(idx,:).';
                    
                    % Apply noise: η fraction of bits are flipped
                    s0 = apply_noise_bipolar(target, eta);
                    
                    % Compute initial overlap (before retrieval)
                    m_initial = mean(s0 .* target);
                    initial_overlap_all(ei,si,mi,tr) = m_initial;

                    % ---- LIF simulation: retrieve from noisy cue ----
                    s_final = simulate_hopfield_LIF_with_spillover(W, s0, kw, params);

                    % Compute final overlap and quality as function of noise
                    m_final = mean(s_final .* target);
                    Q_final = (1 + m_final)/2;  % Q(η) ∈ [0,1]
                    pr = precision_against_active_bits(s_final, target);

                    m_all(ei,si,mi,tr)    = m_final;
                    Q_all(ei,si,mi,tr)    = Q_final;
                    prec_all(ei,si,mi,tr) = pr;

                    prog = prog + 1;
                    if mod(prog, max(1,round(total/50)))==0, fprintf('.'); end
                end
                
                % Report noise-dependent performance
                fprintf('  η=%.2f (σ=%.1f,M=%d): m_init=%.3f → m_final=%.3f, Q(η)=%.3f\n', ...
                    eta, params.sigma_list(si), params.M_list(mi), ...
                    mean(initial_overlap_all(ei,si,mi,:),'all'), ...
                    mean(m_all(ei,si,mi,:),'all'), ...
                    mean(Q_all(ei,si,mi,:),'all'));
            end
        end
    end
    fprintf('\nDone.\n\n');

    % Analyze noise-quality relationship
    fprintf('=== NOISE-QUALITY RELATIONSHIP ===\n');
    for mi = 1:n_M
        for si = 1:n_sigma
            Q_vs_noise = mean(Q_all(:,si,mi,:), 4);
            
            % Find critical noise level (where Q drops below 0.75)
            critical_idx = find(Q_vs_noise < 0.75, 1);
            if ~isempty(critical_idx)
                eta_critical = params.etas(critical_idx);
                fprintf('σ=%.2f, M=%d: Critical noise η_c=%.2f (Q<0.75)\n', ...
                    params.sigma_list(si), params.M_list(mi), eta_critical);
            else
                fprintf('σ=%.2f, M=%d: Robust to all noise levels\n', ...
                    params.sigma_list(si), params.M_list(mi));
            end
        end
    end

    % Package results with noise-dependent analysis
    results = struct();
    results.m_mean    = mean(m_all, 4);
    results.m_std     = std(m_all, 0, 4);
    results.Q_mean    = mean(Q_all, 4);  % Q as function of η
    results.Q_std     = std(Q_all, 0, 4);
    results.prec_mean = mean(prec_all, 4);
    results.prec_std  = std(prec_all, 0, 4);
    results.initial_overlap_mean = mean(initial_overlap_all, 4);
    results.params    = params;
    results.eta_values= params.etas;
    
    % Compute noise robustness metrics
    results.noise_robustness = compute_noise_robustness(results);

    plot_results_with_noise_analysis(results);
end

% =======================================================================
function p = default_params_LIF()
    p.N = 400;  % Smaller network for better retrieval
    p.P = 3;    % Few patterns for strong attractors

    % Noise levels: η=0 (no noise) to η=0.5 (high noise)
    p.etas   = 0:0.1:1;
    p.trials = 100;

    % LIF neuron parameters - tuned for bistability
    p.dt        = 0.5;      % ms
    p.T         = 500;      % ms - longer for convergence
    p.tau_m     = 20;       % membrane time constant
    p.v_rest    = 0;        % Simplified
    p.v_th      = 1.0;      % Threshold
    p.v_reset   = 0;        % Reset
    p.ref_period= 3;        % ms
    p.noise_std = 0.02;     % Small intrinsic noise

    % Cue current - strong initial bias
    p.beta        = 0.5;     % Strong cue
    p.T_clamp     = 100;     % Longer cue period
    p.I_cue_gain  = 3.0;     % Strong gain

    % Recurrent weight scaling
    p.w_scale     = 8.0;     % Strong recurrence for attractors

    % Spillover
    p.sigma_list = [0.0,  1, 2, 3];
    p.M_list     = [100];
    p.alpha      = 50;

    p.seed     = 123;
    p.save_dir = "results_LIF";
end

function merged = merge_with_defaults(user_params, defaults)
    merged = defaults;
    f = fieldnames(user_params);
    for i=1:numel(f)
        merged.(f{i}) = user_params.(f{i});
    end
end

function patterns = generate_patterns_bipolar(N, P)
    X = rand(P, N) > 0.5;
    patterns = 2*X - 1;
end

function W = train_weights(patterns, params)
    [P, N] = size(patterns);
    % Standard Hebbian with strong scaling
    W = params.w_scale * (patterns' * patterns) / N;
    W = W - diag(diag(W));  % No self-connections
end

function s_noisy = apply_noise_bipolar(s, eta)
    % Apply noise by flipping η fraction of bits
    % η=0: no noise (perfect cue)
    % η=0.5: half bits flipped
    N = numel(s); 
    k = round(eta*N);  % Number of bits to flip
    if k > 0
        idx = randperm(N, k);
        s_noisy = s; 
        s_noisy(idx) = -s_noisy(idx);
    else
        s_noisy = s;
    end
end

function kw = make_spillover_kernel(N, sigma, M)
    if sigma <= 0 || M <= 0
        kw = zeros(1,0); return;
    end
    d = (1:M);
    kw = exp(-(d.^2)/(2*sigma^2));
    s = 2*sum(kw);
    if s > 0, kw = kw / s; end
end

% =======================================================================
% ---- LIF dynamics with spillover --------------------------------------
function s_final = simulate_hopfield_LIF_with_spillover(W, s0, kw, params)

    N = numel(s0);
    dt = params.dt;
    steps = round(params.T / dt);

    % Initialize membrane potentials
    v = params.v_rest * ones(N,1);
    tref = zeros(N,1);
    
    % Use a continuous rate variable instead of discrete spikes
    rate = zeros(N,1);  % Firing rate estimate
    
    % Time constant for rate dynamics
    tau_rate = 50;  % ms - smooth filtering
    decay_rate = exp(-dt/tau_rate);

    % Calibrate spillover scale
    alpha_eff = calibrate_alpha(W, kw, params.alpha, N);

    % For decoding: track time-averaged activity
    avg_rate = zeros(N,1);
    avg_count = 0;
    
    for t = 1:steps
        time_ms = t * dt;
        
        % Cue current: strong initially, then fades
        if time_ms <= params.T_clamp
            cue_strength = 1.0 - (time_ms / params.T_clamp) * 0.5;  % Fade gradually
            I_cue = params.I_cue_gain * params.beta * s0 * cue_strength;
        else
            I_cue = zeros(N,1);
        end
        
        % Recurrent input based on smoothed firing rates
        I_syn   = W * rate;
        I_spill = alpha_eff * ring_convolution(rate, kw);

        % Small intrinsic noise
        I_noise = params.noise_std * randn(N,1);

        % Total input
        I_total = I_syn + I_spill + I_cue + I_noise;
        
        % LIF dynamics
        active = (tref <= 0);
        dv = ((-v + I_total) / params.tau_m) * dt;
        v(active) = v(active) + dv(active);

        % Check for spikes
        fired = active & (v >= params.v_th);
        v(fired) = params.v_reset;
        tref(fired) = params.ref_period;
        
        % Decrement refractory counters
        tref(~active) = tref(~active) - dt;

        % Update firing rate with low-pass filter
        spike_signal = double(fired);
        rate = rate * decay_rate + spike_signal * (1 - decay_rate);
        
        % Accumulate for time average (after initial transient)
        if time_ms > params.T_clamp + 50
            avg_rate = avg_rate + rate;
            avg_count = avg_count + 1;
        end
    end

    % Decode pattern from time-averaged activity
    if avg_count > 0
        avg_rate = avg_rate / avg_count;
    end
    
    % Use recurrent field for decoding
    h = W * avg_rate;
    s_final = sign(h);
    s_final(s_final==0) = 1;
end

function spill = ring_convolution(x, kw)
    N = numel(x);
    M = numel(kw);
    spill = zeros(N,1);
    for d = 1:M
        if kw(d)==0, continue; end
        spill = spill + kw(d)*(circshift(x,d) + circshift(x,-d));
    end
end

function alpha_eff = calibrate_alpha(W, kw, r_target, N)
    if r_target <= 0, alpha_eff = 0; return; end
    M = numel(kw);
    if M == 0, alpha_eff = 0; return; end
    
    nSamples = 5;
    stdWs = zeros(nSamples,1);
    stdKs = zeros(nSamples,1);
    for k = 1:nSamples
        s = sign(randn(N,1));
        hs = W*s;
        spill = zeros(N,1);
        for d = 1:M
            if kw(d)==0, continue; end
            spill = spill + kw(d)*(circshift(s,d)+circshift(s,-d));
        end
        stdWs(k) = std(hs);
        stdKs(k) = std(spill);
    end
    sW = mean(stdWs); 
    sK = max(mean(stdKs), eps);
    alpha_eff = r_target * (sW / sK);
end

% =======================================================================
function pr = precision_against_active_bits(s_retrieved, s_target)
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

% =======================================================================
function robustness = compute_noise_robustness(results)
    % Compute metrics that quantify noise robustness
    % - Area under Q(η) curve
    % - Critical noise level η_c where Q drops below threshold
    % - Slope of Q vs η at low noise
    
    etas = results.eta_values;
    [n_eta, n_sigma, n_M] = size(results.Q_mean);
    
    robustness = struct();
    robustness.AUC = zeros(n_sigma, n_M);  % Area under curve
    robustness.eta_critical = zeros(n_sigma, n_M);  % Critical noise
    robustness.slope_low_noise = zeros(n_sigma, n_M);  % Initial slope
    
    for si = 1:n_sigma
        for mi = 1:n_M
            Q_curve = results.Q_mean(:, si, mi);
            
            % Area under Q(η) curve (integral approximation)
            robustness.AUC(si, mi) = trapz(etas, Q_curve);
            
            % Critical noise (where Q < 0.75)
            idx = find(Q_curve < 0.75, 1);
            if ~isempty(idx)
                robustness.eta_critical(si, mi) = etas(idx);
            else
                robustness.eta_critical(si, mi) = max(etas);
            end
            
            % Slope at low noise (first 3 points)
            if n_eta >= 3
                p = polyfit(etas(1:min(3,n_eta)), Q_curve(1:min(3,n_eta)), 1);
                robustness.slope_low_noise(si, mi) = p(1);
            end
        end
    end
end

% =======================================================================
function plot_results_with_noise_analysis(results)
    p = results.params; 
    etas = results.eta_values;
    nM = numel(p.M_list); 
    nS = numel(p.sigma_list); 
    cols = lines(nS);

    figure('Position',[80 80 1400 900]);
    
    % Main plot: Q(η) for different spillover parameters
    for mi = 1:nM
        subplot(2, ceil(nM/2), mi); hold on;
        
        for si = 1:nS
            Qm = results.Q_mean(:,si,mi); 
            Qs = results.Q_std(:,si,mi);
            h = errorbar(etas, Qm, Qs, 'LineWidth',2,'Marker','o',...
                'Color',cols(si,:), 'MarkerSize',6, 'MarkerFaceColor',cols(si,:));
            
            % Add legend label with robustness metric
            AUC = results.noise_robustness.AUC(si, mi);
            h.DisplayName = sprintf('σ=%.2f (AUC=%.2f)', p.sigma_list(si), AUC);
        end
        
        % Reference lines
        plot(etas, 0.5*ones(size(etas)), 'k--', 'LineWidth',1.5, 'DisplayName', 'Chance');
        plot(etas, 0.75*ones(size(etas)), 'r:', 'LineWidth',1, 'DisplayName', 'Threshold');
        
        % Theoretical initial overlap
        m_theory = 1 - 2*etas;
        Q_theory = (1 + m_theory)/2;
        plot(etas, Q_theory, 'k:', 'LineWidth',1, 'DisplayName', 'Initial Q');
        
        xlabel('Noise Level η (flip fraction)', 'FontSize',11); 
        ylabel('Retrieval Quality Q(η)', 'FontSize',11);
        title(sprintf('M=%d, α=%.2f, w-scale=%.1f', p.M_list(mi), p.alpha, p.w_scale), 'FontSize',12);
        legend('Location','southwest', 'FontSize',8);
        grid on; 
        ylim([0.0 1.05]);
        xlim([min(etas) max(etas)]);
    end
    
    sgtitle(sprintf('Noise-Dependent Retrieval Quality: N=%d, P=%d, LIF Hopfield', ...
        p.N, p.P), 'FontSize',14, 'FontWeight','bold');
    
    % Additional analysis figure
    if nS > 1 && nM > 1
        figure('Position',[100 100 1000 400]);
        
        % Plot 1: AUC comparison
        subplot(1,2,1);
        AUC_matrix = results.noise_robustness.AUC;
        imagesc(AUC_matrix');
        colorbar;
        xlabel('Spillover Width σ');
        ylabel('Spillover Range M');
        title('Area Under Q(η) Curve');
        set(gca, 'XTick', 1:nS, 'XTickLabel', p.sigma_list);
        set(gca, 'YTick', 1:nM, 'YTickLabel', p.M_list);
        colormap(hot);
        
        % Plot 2: Critical noise level
        subplot(1,2,2);
        eta_c_matrix = results.noise_robustness.eta_critical;
        imagesc(eta_c_matrix');
        colorbar;
        xlabel('Spillover Width σ');
        ylabel('Spillover Range M');
        title('Critical Noise Level η_c (Q<0.75)');
        set(gca, 'XTick', 1:nS, 'XTickLabel', p.sigma_list);
        set(gca, 'YTick', 1:nM, 'YTickLabel', p.M_list);
        colormap(parula);
        
        sgtitle('Noise Robustness Analysis', 'FontSize',14, 'FontWeight','bold');
    end
end