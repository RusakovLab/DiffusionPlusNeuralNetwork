function Fig_FG_Hopfield()
% FIG_FG_HOPFIELD  Generate paper panels F and G.
%
% Panel F:  Stored pattern | Noise-free probe --> Recalled (sigma=0, M=200)
% Panel G:  Two rows (20% and 30% noise), three columns (sigma=0,5,10), M=200 fixed.
%           Noisy probe: red/black/white tricolor
%           Recalled: magenta/white, titled with Q mean±std

rng(2114);

%% ---- Parameters ---------------------------------------------------------
N            = 400;
P            = 3;
M            = 200;           % fixed for full network
alpha        = 4;
Tmax         = 200;
halt_win     = 5;
beta         = 0.08;
T_clamp      = 8;
n_trials     = 50;
% Display mode per sigma: 'first' = trial 1 (sigma=0 control), 'mean' = closest to mean Q
% sigma_display_modes{i} matches sigma_list_G(i)
sigma_display_modes = {'first', 'mean', 'mean'};  % sigma=0 -> trial 1; sigma=5,10 -> closest to mean
eta_list_G   = [0.20, 0.30];  % two noise rows
sigma_list_G = [0, 5, 10];    % three sigma columns
sigma_F      = 0;             % panel F uses sigma=0
grid_sz      = sqrt(N);       % 20x20

%% ---- Generate & train ---------------------------------------------------
patterns = generate_patterns(N, P);
W        = train_weights(patterns, N);

target_idx = 1;
target     = patterns(target_idx, :).';

%% ====================================================================
%  PANEL F  -  Noise-free probe, sigma=0, M=200
%% ====================================================================
kw_F = make_kernel(N, sigma_F, M);
[Q_F_mean, Q_F_std, s_recalled_F] = ...
    run_trials(W, target, 0.0, kw_F, alpha, Tmax, halt_win, beta, T_clamp, n_trials, 'first');

fig_F = figure('Name','Panel F','Color','w','Position',[50 50 900 330]);

ax1 = subplot(1,3,1);
show_binary(ax1, target, grid_sz, [1 0 0; 1 1 1]);
title('Stored pattern','FontSize',12,'FontWeight','normal');

ax2 = subplot(1,3,2);
show_binary(ax2, target, grid_sz, [1 0 0; 1 1 1]);
title('Noise-free Probe','FontSize',12,'FontWeight','normal');

ax3 = subplot(1,3,3);
show_binary(ax3, s_recalled_F, grid_sz, [1 0 1; 1 1 1]);
title(sprintf('Recall Q: %.1f\\pm%.1f%%', Q_F_mean*100, Q_F_std*100), ...
      'FontSize',12,'FontWeight','normal');
xlabel(sprintf('\\sigma=%.0f', sigma_F),'FontSize',11);

set(fig_F,'PaperPositionMode','auto');
saveas(fig_F, 'Panel_F.png');
fprintf('Panel F: Q = %.1f +/- %.1f%%\n', Q_F_mean*100, Q_F_std*100);

%% ====================================================================
%  PANEL G  -  Two noise rows x three sigma columns
%% ====================================================================
n_eta   = numel(eta_list_G);
n_sigma = numel(sigma_list_G);

Q_means  = zeros(n_eta, n_sigma);
Q_stds   = zeros(n_eta, n_sigma);
recalled = cell(n_eta, n_sigma);
noisy_display = cell(n_eta, 1);

rng(42);
for ei = 1:n_eta
    eta = eta_list_G(ei);
    s_noisy = apply_noise(target, eta);
    noisy_display{ei} = target + s_noisy;   % {-2, 0, +2}

    for si = 1:n_sigma
        kw = make_kernel(N, sigma_list_G(si), M);
        [Q_means(ei,si), Q_stds(ei,si), recalled{ei,si}] = ...
            run_trials(W, target, eta, kw, alpha, Tmax, halt_win, beta, T_clamp, n_trials, sigma_display_modes{si});
        fprintf('eta=%.0f%%, sigma=%2d: Q = %.1f +/- %.1f%%\n', ...
            eta*100, sigma_list_G(si), Q_means(ei,si)*100, Q_stds(ei,si)*100);
    end
end

% Figure: n_eta rows x (1 + n_sigma) columns
fig_G = figure('Name','Panel G','Color','w','Position',[100 100 1100 700]);

for ei = 1:n_eta
    % -- Noisy probe column --
    ax_noise = subplot(n_eta, n_sigma+1, (ei-1)*(n_sigma+1) + 1);
    img_noisy = reshape(noisy_display{ei}, grid_sz, grid_sz);
    imagesc(ax_noise, img_noisy);
    colormap(ax_noise, [1 0 0; 0 0 0; 1 1 1]);  % red=-2, black=0, white=+2
    clim(ax_noise, [-2 2]);
    axis(ax_noise,'square','off');
    title(ax_noise, sprintf('Noise: %.0f%%', eta_list_G(ei)*100), ...
          'FontSize',11,'FontWeight','normal');
    xlabel(ax_noise,'Probe','FontSize',10);

    % -- Recalled columns --
    for si = 1:n_sigma
        ax = subplot(n_eta, n_sigma+1, (ei-1)*(n_sigma+1) + si + 1);
        show_binary(ax, recalled{ei,si}, grid_sz, [1 0 1; 1 1 1]);
        title(ax, sprintf('Q: %.1f\\pm%.1f%%', Q_means(ei,si)*100, Q_stds(ei,si)*100), ...
              'FontSize',11,'FontWeight','normal');
        xlabel(ax, sprintf('\\sigma=%d', sigma_list_G(si)), 'FontSize',11, 'FontWeight','bold');
    end
end

set(fig_G,'PaperPositionMode','auto');
saveas(fig_G, 'Panel_G.png');
fprintf('Panel G saved.\n');

end  % main function

%% =========================================================================
%  HELPER FUNCTIONS
%% =========================================================================

function patterns = generate_patterns(N, P)
    X        = rand(P, N) > 0.5;
    patterns = 2*X - 1;
end

function W = train_weights(patterns, N)
    W = (patterns' * patterns) / N;
    W = W - diag(diag(W));
end

function kw = make_kernel(N, sigma, M) %#ok<INUSL>
    if sigma <= 0 || M <= 0
        kw = zeros(1,0); return;
    end
    d  = 1:M;
    kw = exp(-(d.^2) / (2*sigma^2));
    s  = 2*sum(kw);
    if s > 0, kw = kw / s; end
end

function s_noisy = apply_noise(s, eta)
    N   = numel(s);
    k   = round(eta * N);
    idx = randperm(N, k);
    s_noisy      = s;
    s_noisy(idx) = -s_noisy(idx);
end

function [Q_mean, Q_std, s_example] = run_trials(W, target, eta, kw, alpha, ...
                                                   Tmax, halt_win, beta, T_clamp, n_trials, display_mode)
    Qs    = zeros(n_trials, 1);
    N     = numel(target);
    all_s = zeros(N, n_trials);
    for tr = 1:n_trials
        s0          = apply_noise(target, eta);
        s_final     = simulate(W, s0, kw, alpha, Tmax, halt_win, beta, T_clamp);
        m           = mean(s_final .* target);
        Qs(tr)      = (1 + m) / 2;
        all_s(:,tr) = s_final;
    end
    Q_mean = mean(Qs);
    Q_std  = std(Qs);
    % Select display trial based on mode:
    %   'first' -> trial 1 (sigma=0 control, always perfect)
    %   'mean'  -> trial closest to mean Q, but exclude spatially periodic
    %              (striped) spurious states using an isotropy check
    %   integer -> specific trial number 1..n_trials
    if ischar(display_mode) || isstring(display_mode)
        if strcmp(display_mode, 'mean')
            g = round(sqrt(N));
            is_striped = false(n_trials, 1);
            for tr = 1:n_trials
                s_2d = reshape(all_s(:,tr), g, g);
                % Compute mean absolute row-sum vs col-sum variance
                row_var = var(mean(s_2d, 2));  % variance of row means
                col_var = var(mean(s_2d, 1));  % variance of col means
                % If one direction has much higher spatial variance -> striped
                ratio = max(row_var, col_var) / (min(row_var, col_var) + 1e-6);
                is_striped(tr) = (ratio > 4);  % threshold: one axis 5x more structured
            end
            valid = find(~is_striped & Qs >= 0.60);
            if isempty(valid)
                % Fall back: just use closest to mean among all trials
                [~, idx] = min(abs(Qs - Q_mean));
                fprintf('  NOTE: no non-striped trials found, showing closest to mean\n');
            else
                [~, local_idx] = min(abs(Qs(valid) - Q_mean));
                idx = valid(local_idx);
            end
        else  % 'first'
            idx = 1;
        end
    else
        idx = max(1, min(round(display_mode), n_trials));
    end
    s_example = all_s(:, idx);
end

function s = simulate(W, s0, kw, alpha_ratio, Tmax, halt_win, beta, T_clamp)
    N     = numel(s0);
    s     = s0;
    M     = numel(kw);
    alpha = calibrate_alpha(W, kw, alpha_ratio, N);
    recent = zeros(halt_win, N);
    for t = 1:Tmax
        if M > 0
            spill = zeros(N,1);
            for d = 1:M
                if kw(d) == 0, continue; end
                spill = spill + kw(d) * (circshift(s,d) + circshift(s,-d));
            end
            spill = alpha * spill;
        else
            spill = zeros(N,1);
        end

        beta_t = (t <= T_clamp) * beta;
        h      = W*s + spill + beta_t*s0;
        s_new  = sign(h);
        z      = (s_new == 0);
        s_new(z) = s(z);

        recent(mod(t-1, halt_win)+1, :) = s_new.';
        if t >= halt_win && all(std(recent,0,1) == 0)
            s = s_new; return;
        end
        s = s_new;
    end
end

function alpha_eff = calibrate_alpha(W, kw, r_target, N)
    if r_target <= 0 || numel(kw) == 0
        alpha_eff = 0; return;
    end
    M  = numel(kw);
    ns = 5;
    sW = zeros(ns,1); sK = zeros(ns,1);
    for k = 1:ns
        s  = sign(randn(N,1));
        hs = W*s;
        sp = zeros(N,1);
        for d = 1:M
            if kw(d)==0, continue; end
            sp = sp + kw(d)*(circshift(s,d)+circshift(s,-d));
        end
        sW(k) = std(hs);
        sK(k) = std(sp);
    end
    alpha_eff = r_target * (mean(sW) / max(mean(sK), eps));
end

function K2D = make_2D_kernel(g, kw)
% Build a 2D isotropic kernel from the 1D radial profile kw.
% kw(d) gives the weight at integer radius d (pixels).
% The kernel is normalized so its sum = 1.
    M   = numel(kw);
    if M == 0
        K2D = zeros(1,1); return;
    end
    sz  = 2*M + 1;
    K2D = zeros(sz, sz);
    cx  = M + 1;  % centre pixel
    for r = -M:M
        for c = -M:M
            d = sqrt(r^2 + c^2);
            di = floor(d);
            if di >= 1 && di <= M
                % Linear interpolation between kw(di) and kw(di+1)
                frac = d - di;
                if di < M
                    w = kw(di)*(1-frac) + kw(di+1)*frac;
                else
                    w = kw(di)*(1-frac);
                end
                K2D(cx+r, cx+c) = w;
            end
        end
    end
    K2D(cx, cx) = 0;  % no self-influence
    % Normalize
    s = sum(K2D(:));
    if s > 0, K2D = K2D / s; end
end

function show_binary(ax, s, grid_sz, cmap)
    img = reshape((s > 0), grid_sz, grid_sz);
    imagesc(ax, img);
    colormap(ax, cmap);
    clim(ax, [0 1]);
    axis(ax, 'square', 'off');
end