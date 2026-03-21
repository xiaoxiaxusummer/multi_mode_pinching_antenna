% ------------------------------------------------------------
% Multi-mode PASS simulation:
%   Sum-rate vs number of PAs/antennas N at fixed transmit power.
%   Methods compared (all averaged over Monte-Carlo user drops):
%   - Proposed Case-1: multi-mode PASS (mode selection) 
%   - Proposed Case-2: multi-mode PASS (mode combining) 
%   - Proposed Uniform Mode Combining: fixed PA propagation constant betaPA = (beta1+beta2)/2
%   - Baseline: single-mode PASS + TDMA
%   - Baseline: conventional MISO with I antennas + Hybrid BF (RF+BB)
% ------------------------------------------------------------

clc; clear; close all;
addpath(genpath(pwd));

%% ====================== Global parameters ======================
c = 3e8;
fc = 28e9;                         % frequency
lambda = c/fc;
k0 = 2*pi/lambda;

% Waveguide / geometry (required)
Lwg  = 20;                         % length of waveguide
xmin = 0;  xmax = Lwg;

% PA
Lpa  = 6e-3;                       % length of PA
hPA  = 2.5;                        % height of PA
dmin = lambda/2;                   % min spacing

% Noise
B = 100e6;
N0_dBmHz = -174;
sigma2_dBm = N0_dBmHz + 10*log10(B);
sigma2 = 1e-3 * 10^(sigma2_dBm/10);


% -------- Fixed transmit power --------
P_tx_dBm = 25;                    
Pmax = 1e-3 * 10^(P_tx_dBm/10);

%% ====================== Two guided modes ======================
epsr_core = 4.0;
n_core = sqrt(epsr_core);
n_clad = 1.0;
h = 8e-3;      % rectangular strip height = 7 mm
w = 4e-3;      % width  = 4 mm

% --- Effective Index Method (EIM) ---
% Step 1: vertical slab (thickness h), TE order py=0 -> n_eff_y0
py0 = 0;
n_eff_y0 = slab_neff_TE(n_core, n_clad, h, fc, py0);

% Step 2: horizontal slab (width w), TE orders px=0 and px=1
px0 = 0;   % quasi-TE0
px1 = 1;   % quasi-TE1
n_eff_qTE0 = slab_neff_TE(n_eff_y0, n_clad, w, fc, px0);
n_eff_qTE1 = slab_neff_TE(n_eff_y0, n_clad, w, fc, px1);

beta1 = k0 * n_eff_qTE0;   % quasi-TE0
beta2 = k0 * n_eff_qTE1;   % quasi-TE1
beta  = [beta1; beta2];

fprintf('--- Mode propagation constants at fc=%.1f GHz (dielectric rectangular strip) ---\n', fc/1e9);
fprintf('epsr_core=%.2f (n_core=%.3f), n_clad=%.3f\n', epsr_core, n_core, n_clad);
fprintf('h = %.1f mm, w = %.1f mm\n', h*1e3, w*1e3);

if isnan(n_eff_y0) || isnan(n_eff_qTE0) || isnan(n_eff_qTE1)
    error('EIM solver returned NaN (mode not guided for this geometry).');
else
    fprintf('n_eff(quasi-TE0) = %.6f\n', n_eff_qTE0);
    fprintf('n_eff(quasi-TE1) = %.6f\n', n_eff_qTE1);
    fprintf('beta(quasi-TE0)  = %.6f rad/m\n', beta1);
    fprintf('beta(quasi-TE1)  = %.6f rad/m\n', beta2);
end

%% ====================== Scenario config ======================
K = 2; M = 2;


%% ====================== PSO params ======================
PSO = struct();
% PSO.P   = 50;              % population size
% PSO.T   = 100;              % maximum iterations
% PSO.w   = 0.72;
% PSO.c1  = 1.6;
% PSO.c2  = 1.5;
% PSO.vmax_x_factor    = 1;
% PSO.vmax_beta_factor = 1;
% PSO.vmax_loglambda   = 1;
% PSO.vmax_logp        = 1;
% PSO.radius           = inf;     % no local clamp
% PSO.warmstart_case2    = true;
% PSO.case2_elite_frac = 0.20;
% % KPBF numerical bounds
% cfg.lambdaBF_min = 1e-4;
% cfg.lambdaBF_max = 1e4;

PSO.P   = 50;              % population size
PSO.T   = 200;              % maximum iterations
PSO.w   = 0.72;
PSO.c1  = 1.8;
PSO.c2  = 1.2;
PSO.vmax_x_factor    = 1;
PSO.vmax_beta_factor = 1;
PSO.vmax_loglambda   = 1;
PSO.vmax_logp        = 1;
PSO.radius           = inf;     % no local clamp
PSO.warmstart_case2    = true;
PSO.case2_elite_frac = 0.20;
% KPBF numerical bounds
cfg.lambdaBF_min = 1e-3;
cfg.lambdaBF_max = 1e3;

%% ====================== Sweep over I ======================
N_list = 4:2:16;
baseSeed = 6;
test_num = 100; % Monte Carlo trials
rng(baseSeed,'twister');
x_user_min = 5;
a_all = x_user_min*ones(2,test_num) + rand(2,test_num)*(Lwg-x_user_min);
y_all = rand(2,test_num)*5;
% Uniform mode combining
betaPA_mid = (beta1 + beta2)/2;
% Result arrays (I x trial)
R_case1 = zeros(numel(N_list), test_num);
R_case2 = zeros(numel(N_list), test_num);
R_midbeta = zeros(numel(N_list), test_num);
R_tdma1 = zeros(numel(N_list), test_num);
R_mimo_hybrid = zeros(numel(N_list), test_num);

for it = 1:test_num
    a = a_all(:,it);                          % user x 
    y = y_all(:,it);                          % user y
    
    parfor iN = 1:numel(N_list)
    % for iN = 1:numel(N_list)
        N = N_list(iN);

        x_arr = dmin*(0:1:N-1);
        tx_pos = [x_arr(:), zeros(N,1), hPA*ones(N,1)];

        % D = 0.2;   % fixed aperture for fixed MIMO
        % x_arr = linspace(-D/2, D/2, N);
        % tx_pos = [x_arr(:), zeros(N,1), hPA*ones(N,1)];

        cfg = struct();
        cfg.lambda=lambda; cfg.k0=k0; cfg.beta=beta;
        cfg.xmin=xmin; cfg.xmax=xmax; cfg.dmin=dmin;
        cfg.hPA=hPA; cfg.a=a; cfg.y=y;
        cfg.Pmax=Pmax; cfg.sigma2=sigma2; 
        cfg.N=N; cfg.K=K; cfg.M=M;
        cfg.Lpa=Lpa;
        cfg.kappa_fixed = (pi/6)/cfg.Lpa;  % sin(kappa*Lpa)=sin(pi/6)=0.5 is the maximum radiation ratio

        cfg.kappa_match = cfg.kappa_fixed;

        % bounds for betaPA (used by case1/case2)
        cfg.betaPA_fixed = [];
        [cfg.betaPA_min, cfg.betaPA_max] = default_betaPA_bounds(cfg);

        % -------------------- Proposed: Case 1 Mode selection--------------------
        x_init = init_random_allx(cfg);
        beta_init = cfg.betaPA_min + rand(cfg.N,1).*(cfg.betaPA_max - cfg.betaPA_min);
        lambda_init = ones(cfg.K,1);
        p_init = ones(cfg.K,1);
        PSO_local = PSO; 
        out1 = pso_KPBF_multiPA_case1_optx_beta_lambda_p(cfg, PSO_local, x_init, beta_init, lambda_init, p_init);
        R_case1(iN,it) = out1.best_sr;

        % -------------------- Proposed: Case 2 Mode combining --------------------
        if PSO.warmstart_case2
            x_init2 = out1.best_x;
            beta_init2 = out1.best_betaPA;
            lambda_init2 = out1.best_lambda;
            p_init2 = out1.best_p_rel;
        else
            x_init2 = x_init;
            beta_init2 = cfg.betaPA_min + rand(cfg.N,1).*(cfg.betaPA_max - cfg.betaPA_min);
            lambda_init2 = ones(cfg.K,1);
            p_init2 = ones(cfg.K,1);
        end
        out2 = pso_KPBF_multiPA_case2_optx_beta_lambda_p(cfg, PSO_local, x_init2, beta_init2, lambda_init2, p_init2);
        R_case2(iN,it) = out2.best_sr;

        % -------------------- Proposed: (fixed) uniform mode combining --------------------
        cfg.betaPA_fixed = betaPA_mid * ones(cfg.N,1);
        cfg.betaPA_min = cfg.betaPA_fixed;  % lock
        cfg.betaPA_max = cfg.betaPA_fixed;

        cfg.case_id = 1;
        lambda_init_b = ones(cfg.K,1);
        p_init_b = ones(cfg.K,1);
        out_mid = pso_KPBF_multiPA_optx_lambda_p_fixedbeta(cfg, PSO_local, x_init, lambda_init, p_init);
        R_midbeta(iN,it) = out_mid.best_sr;

        % -------------------- TDMA baseline --------------------
        x_TDMA = cfg.a;
        R_tdma1(iN,it) = tdma_singlemode_PASS(cfg, x_TDMA);

        % -------------------- Conventional MISO (Hybrid BF) baseline --------------------
        % Hmimo = build_H_free(cfg,tx_pos(:,1));
        Hmimo = build_H_miso(cfg,tx_pos(:,1));
        [~, F_RF_best, F_BB_best] = hybridBF_MIMO(Hmimo, Pmax, sigma2);
        R_mimo_hybrid(iN,it) = sumrate_from_HW(Hmimo, F_RF_best*F_BB_best, sigma2);

        fprintf('[it=%02d/%02d] N=%2d, Ptx=%4.1f: C1=%.3f, C2=%.3f, C1-mid=%.3f, TDMA=%.3f, MISO=%.3f\n', ...
            it, test_num, N, Pmax, R_case1(iN,it), R_case2(iN,it), R_midbeta(iN,it), R_tdma1(iN,it), R_mimo_hybrid(iN,it));
    end
end

%% ====================== Plot ======================
figure('Color',[1,1,1]); hold on; box on;
plot(N_list, mean(R_case1,2), '--ob', 'LineWidth', 2.5);
plot(N_list, mean(R_case2,2), '-sb', 'LineWidth', 2.5);
plot(N_list, mean(R_midbeta,2), '--hb', 'LineWidth', 2.5);
plot(N_list, mean(R_tdma1,2), '-.^r', 'LineWidth', 2.5);
plot(N_list, mean(R_mimo_hybrid,2), '-dk', 'LineWidth', 2.5);

xlabel('Number of PAs / antennas, I','Interpreter','latex');
ylabel('Sum rate (bps/Hz)','Interpreter','latex');
% title(sprintf('Dual-mode PASS @ %.0f GHz, $P_{tx}$ = %.1f dBm', fc/1e9, P_tx_dBm), 'Interpreter','latex');

legend({ ...
    'Multi-mode PASS (mode selection)', ...
    'Multi-mode PASS (mode combining)', ...
    'Multi-mode PASS (uniform mode combining)', ...
    'Single-mode PASS (TDMA)', ...
    'Conventional MISO (Hybrid BF)' ...
    }, 'Location','best');

set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [22, 16], 'PaperPositionMode','auto','Renderer','painters');
set(gca,'FontSize',14,'FontName','Times New Roman');

if ~exist('figs','dir'), mkdir('figs'); end
savefig(sprintf('figs/rate_N_Ptx%.0fdBm_multiPA.fig', P_tx_dBm));
print(gcf, sprintf('figs/rate_N_Ptx%.0fdBm_multiPA.pdf', P_tx_dBm), '-dpdf','-painters');

% Save raw data
save(sprintf('figs/rate_N_Ptx%.0fdBm.mat', P_tx_dBm), ...
    'N_list','P_tx_dBm','Pmax','R_case1','R_case2','R_midbeta','R_tdma1','R_mimo_hybrid');

function out = pso_KPBF_multiPA_case1_optx_beta_lambda_p(cfg, PSO, x0, beta0, lambda0, p0)
N = cfg.N; K = cfg.K;
% --- per-variable velocity caps (dimension-aware) ---
vmax_x   = PSO.vmax_x_factor * abs(cfg.xmax - cfg.xmin);
vmax_beta= PSO.vmax_beta_factor * abs(cfg.beta(2)-cfg.beta(1)) + 1e-9;
vmax_ll  = PSO.vmax_loglambda;
vmax_lp  = PSO.vmax_logp;
% lambda bounds (for numerical stability)
if ~isfield(cfg,'lambdaBF_min') || isempty(cfg.lambdaBF_min), cfg.lambdaBF_min = 1e-4; end
if ~isfield(cfg,'lambdaBF_max') || isempty(cfg.lambdaBF_max), cfg.lambdaBF_max = 1e4;  end

Xx = zeros(PSO.P, N); Vx = zeros(PSO.P, N);
Xb = zeros(PSO.P, N); Vb = zeros(PSO.P, N);
Xl = zeros(PSO.P, K); Vl = zeros(PSO.P, K);
Xp = zeros(PSO.P, K); Vp = zeros(PSO.P, K);
Xx0 = repmat(x0(:).', PSO.P, 1); % common clamp center (original behavior)

Xx0 = zeros(PSO.P, N); % per-particle init center for optional radius clamp

for pp=1:PSO.P
    % x init
    jitter_x = (2*rand(1,N)-1) * (0.5*cfg.lambda);
    Xx(pp,:) = proj_and_repair(x0 + jitter_x, cfg.xmin, cfg.xmax, cfg.dmin);
    Vx(pp,:) = (2*rand(1,N)-1) * (0.1*cfg.lambda);

    % betaPA init
    jitter_b = (2*rand(1,N)-1) .* (0.05*abs(cfg.beta(2)-cfg.beta(1)) + 1e-9);
    Xb(pp,:) = clamp_vec(beta0(:).' + jitter_b, cfg.betaPA_min, cfg.betaPA_max);
    Vb(pp,:) = (2*rand(1,N)-1) .* (0.02*abs(cfg.beta(2)-cfg.beta(1)) + 1e-9);

    % lambda init (optimize in log-domain)
    l0 = max(cfg.lambdaBF_min, min(cfg.lambdaBF_max, lambda0(:)));
    z0 = log(l0(:).');
    Xl(pp,:) = z0 + (2*rand(1,K)-1)*0.2;      % log-lambda
    Vl(pp,:) = (2*rand(1,K)-1)*0.1;

    % power init (log-domain; mapped to simplex by softmax)
    p0v = max(1e-8, p0(:));
    zP0 = log(p0v(:).');
    Xp(pp,:) = zP0 + (2*rand(1,K)-1)*0.2;
    Vp(pp,:) = (2*rand(1,K)-1)*0.1;
end

pbest_x = Xx; pbest_b = Xb; pbest_l = Xl; pbest_p = Xp;
pbest_val = -inf(PSO.P,1);
gbest_x = Xx(1,:); gbest_b = Xb(1,:); gbest_l = Xl(1,:); gbest_p = Xp(1,:);
gbest_val = -inf;

for t=1:PSO.T
    for pp=1:PSO.P
        lambda_vec = exp(clamp_vec(Xl(pp,:), log(cfg.lambdaBF_min), log(cfg.lambdaBF_max))).';
        p_rel = softmax_vec(Xp(pp,:));
        f = fitness_sumrate_KPBF(cfg, Xx(pp,:), Xb(pp,:), lambda_vec, p_rel);

        if f > pbest_val(pp)
            pbest_val(pp) = f;
            pbest_x(pp,:) = Xx(pp,:);
            pbest_b(pp,:) = Xb(pp,:);
            pbest_l(pp,:) = Xl(pp,:);
            pbest_p(pp,:) = Xp(pp,:);
        end
        if f > gbest_val
            gbest_val = f;
            gbest_x = Xx(pp,:);
            gbest_b = Xb(pp,:);
            gbest_l = Xl(pp,:);
            gbest_p = Xp(pp,:);
        end
    end

    for pp=1:PSO.P
        r1x = rand(1,N); r2x = rand(1,N);
        r1k = rand(1,K); r2k = rand(1,K);

        % ---- x update ----
        Vx(pp,:) = PSO.w*Vx(pp,:) ...
            + PSO.c1*(r1x.*(pbest_x(pp,:)-Xx(pp,:))) ...
            + PSO.c2*(r2x.*(gbest_x -Xx(pp,:)));
        Vx(pp,:) = clip_vec(Vx(pp,:), vmax_x);
        Xx(pp,:) = proj_and_repair(Xx(pp,:) + Vx(pp,:), cfg.xmin, cfg.xmax, cfg.dmin);
        Xx(pp,:) = clamp_around_init(Xx(pp,:), Xx0(pp,:), PSO.radius);
        Xx(pp,:) = proj_and_repair(Xx(pp,:), cfg.xmin, cfg.xmax, cfg.dmin);

        % ---- betaPA update ----
        Vb(pp,:) = PSO.w*Vb(pp,:) ...
            + PSO.c1*(r1x.*(pbest_b(pp,:)-Xb(pp,:))) ...
            + PSO.c2*(r2x.*(gbest_b -Xb(pp,:)));
        Vb(pp,:) = clip_vec_beta(Vb(pp,:), vmax_beta);
        Xb(pp,:) = clamp_vec(Xb(pp,:) + Vb(pp,:), cfg.betaPA_min, cfg.betaPA_max);

        % ---- log-lambda update ----
        Vl(pp,:) = PSO.w*Vl(pp,:) ...
            + PSO.c1*(r1k.*(pbest_l(pp,:)-Xl(pp,:))) ...
            + PSO.c2*(r2k.*(gbest_l -Xl(pp,:)));
        Vl(pp,:) = clamp_vec(Vl(pp,:), -vmax_ll, vmax_ll);
        Xl(pp,:) = Xl(pp,:) + Vl(pp,:);
        Xl(pp,:) = clamp_vec(Xl(pp,:), log(cfg.lambdaBF_min), log(cfg.lambdaBF_max));

        % ---- log-power update ----
        Vp(pp,:) = PSO.w*Vp(pp,:) ...
            + PSO.c1*(r1k.*(pbest_p(pp,:)-Xp(pp,:))) ...
            + PSO.c2*(r2k.*(gbest_p -Xp(pp,:)));
        Vp(pp,:) = clamp_vec(Vp(pp,:), -vmax_lp, vmax_lp);
        Xp(pp,:) = Xp(pp,:) + Vp(pp,:);
        % no hard bounds; softmax will normalize
    end
end

best_lambda = exp(clamp_vec(gbest_l, log(cfg.lambdaBF_min), log(cfg.lambdaBF_max))).';
best_p_rel  = softmax_vec(gbest_p);
[best_sr, best_W, best_sinr, Heff] = eval_KPBF_new(cfg, gbest_x, gbest_b, best_lambda, best_p_rel);

out.best_x = gbest_x;
out.best_betaPA = gbest_b(:);
out.best_lambda = best_lambda;
out.best_p_rel  = best_p_rel(:);
out.best_sr = best_sr;
out.best_W = best_W;
out.best_sinr = best_sinr;
out.best_Heff = Heff;
end



function out = pso_KPBF_multiPA_case2_optx_beta_lambda_p(cfg, PSO, x0, beta0, lambda0, p0)
% Case-2: PSO jointly optimizes x, betaPA, lambda_k, p_k, with an "elite + random" initialization:
%   - First Ne elite particles: jittered warm-start around the provided (x0,beta0,lambda0,p0)
%                             (typically Case-1 best), to refine the known good basin.
%   - Remaining particles: random initialization for global exploration.
%
% This is usually better than "all warm-start" (swarm collapse) or "all random" (no baseline).

I = cfg.N; K = cfg.K;


% --- per-variable velocity caps (dimension-aware) ---
vmax_x   = PSO.vmax_x_factor * cfg.lambda;
vmax_beta= PSO.vmax_beta_factor * abs(cfg.beta(2)-cfg.beta(1)) + 1e-9;
vmax_ll  = PSO.vmax_loglambda;
vmax_lp  = PSO.vmax_logp;
% lambda bounds (for numerical stability)
if ~isfield(cfg,'lambdaBF_min') || isempty(cfg.lambdaBF_min), cfg.lambdaBF_min = 1e-4; end
if ~isfield(cfg,'lambdaBF_max') || isempty(cfg.lambdaBF_max), cfg.lambdaBF_max = 1e4;  end

Xx = zeros(PSO.P, I); Vx = zeros(PSO.P, I);
Xb = zeros(PSO.P, I); Vb = zeros(PSO.P, I);
Xl = zeros(PSO.P, K); Vl = zeros(PSO.P, K);
Xp = zeros(PSO.P, K); Vp = zeros(PSO.P, K);
Xx0 = zeros(PSO.P, I); % per-particle init center for optional radius clamp

% ---- Elite warm-start particles (jittered around x0/beta0/lambda0/p0) ----
if ~isfield(PSO,'case2_elite_frac') || isempty(PSO.case2_elite_frac)
    PSO.case2_elite_frac = 0.08;
end
Ne = max(1, min(PSO.P, round(PSO.P * PSO.case2_elite_frac)));

l0 = max(cfg.lambdaBF_min, min(cfg.lambdaBF_max, lambda0(:)));
z0 = log(l0(:).');
p0v = max(1e-8, p0(:));
zP0 = log(p0v(:).');

for pp = 1:Ne
    % x
    jitter_x = (2*rand(1,I)-1) * (0.15*cfg.lambda);
    Xx(pp,:) = proj_and_repair(x0(:).' + jitter_x, cfg.xmin, cfg.xmax, cfg.dmin);
    Xx0(pp,:) = Xx(pp,:);
    Vx(pp,:) = (2*rand(1,I)-1) * (0.05*cfg.lambda);

    % betaPA
    jitter_b = (2*rand(1,I)-1) .* (0.03*abs(cfg.beta(2)-cfg.beta(1)) + 1e-9);
    Xb(pp,:) = clamp_vec(beta0(:).' + jitter_b, cfg.betaPA_min, cfg.betaPA_max);
    Vb(pp,:) = (2*rand(1,I)-1) .* (0.015*abs(cfg.beta(2)-cfg.beta(1)) + 1e-9);

    % lambda (log-domain)
    Xl(pp,:) = z0 + (2*rand(1,K)-1)*0.10;
    Vl(pp,:) = (2*rand(1,K)-1)*0.05;

    % power (log-domain)
    Xp(pp,:) = zP0 + (2*rand(1,K)-1)*0.10;
    Vp(pp,:) = (2*rand(1,K)-1)*0.05;
end

% ---- Remaining particles: random init (broad exploration) ----
for pp=(Ne+1):PSO.P
    % x random
    x_rand = init_random_allx(cfg);
    jitter_x = (2*rand(1,I)-1) * (0.5*cfg.lambda);
    Xx(pp,:) = proj_and_repair(x_rand(:).' + jitter_x, cfg.xmin, cfg.xmax, cfg.dmin);
    Xx0(pp,:) = Xx(pp,:);
    Vx(pp,:) = (2*rand(1,I)-1) * (0.1*cfg.lambda);

    % betaPA random
    b_rand = cfg.betaPA_min + rand(1,I).*(cfg.betaPA_max - cfg.betaPA_min);
    jitter_b = (2*rand(1,I)-1) .* (0.05*abs(cfg.beta(2)-cfg.beta(1)) + 1e-9);
    Xb(pp,:) = clamp_vec(b_rand + jitter_b, cfg.betaPA_min, cfg.betaPA_max);
    Vb(pp,:) = (2*rand(1,I)-1) .* (0.02*abs(cfg.beta(2)-cfg.beta(1)) + 1e-9);

    % lambda random (log-domain)
    l_rand = cfg.lambdaBF_min * (cfg.lambdaBF_max/cfg.lambdaBF_min) .^ rand(1,K); % log-uniform
    Xl(pp,:) = log(l_rand) + (2*rand(1,K)-1)*0.2;
    Vl(pp,:) = (2*rand(1,K)-1)*0.1;

    % power random (log-domain; softmax -> simplex)
    Xp(pp,:) = randn(1,K);
    Vp(pp,:) = (2*rand(1,K)-1)*0.1;
end

pbest_x = Xx; pbest_b = Xb; pbest_l = Xl; pbest_p = Xp;
pbest_val = -inf(PSO.P,1);
gbest_x = Xx(1,:); gbest_b = Xb(1,:); gbest_l = Xl(1,:); gbest_p = Xp(1,:);
gbest_val = -inf;

for t=1:PSO.T
    for pp=1:PSO.P
        lambda_vec = exp(clamp_vec(Xl(pp,:), log(cfg.lambdaBF_min), log(cfg.lambdaBF_max))).';
        p_rel = softmax_vec(Xp(pp,:));
        f = fitness_sumrate_KPBF(cfg, Xx(pp,:), Xb(pp,:), lambda_vec, p_rel);

        if f > pbest_val(pp)
            pbest_val(pp) = f;
            pbest_x(pp,:) = Xx(pp,:);
            pbest_b(pp,:) = Xb(pp,:);
            pbest_l(pp,:) = Xl(pp,:);
            pbest_p(pp,:) = Xp(pp,:);
        end
        if f > gbest_val
            gbest_val = f;
            gbest_x = Xx(pp,:);
            gbest_b = Xb(pp,:);
            gbest_l = Xl(pp,:);
            gbest_p = Xp(pp,:);
        end
    end

    for pp=1:PSO.P
        r1x = rand(1,I); r2x = rand(1,I);
        r1k = rand(1,K); r2k = rand(1,K);

        % ---- x update ----
        Vx(pp,:) = PSO.w*Vx(pp,:) ...
            + PSO.c1*(r1x.*(pbest_x(pp,:)-Xx(pp,:))) ...
            + PSO.c2*(r2x.*(gbest_x -Xx(pp,:)));
        Vx(pp,:) = clip_vec(Vx(pp,:), vmax_x);
        Xx(pp,:) = proj_and_repair(Xx(pp,:) + Vx(pp,:), cfg.xmin, cfg.xmax, cfg.dmin);
        Xx(pp,:) = clamp_around_init(Xx(pp,:), Xx0(pp,:), PSO.radius);
        Xx(pp,:) = proj_and_repair(Xx(pp,:), cfg.xmin, cfg.xmax, cfg.dmin);

        % ---- betaPA update ----
        Vb(pp,:) = PSO.w*Vb(pp,:) ...
            + PSO.c1*(r1x.*(pbest_b(pp,:)-Xb(pp,:))) ...
            + PSO.c2*(r2x.*(gbest_b -Xb(pp,:)));
        Vb(pp,:) = clip_vec_beta(Vb(pp,:), vmax_beta);
        Xb(pp,:) = clamp_vec(Xb(pp,:) + Vb(pp,:), cfg.betaPA_min, cfg.betaPA_max);

        % ---- log-lambda update ----
        Vl(pp,:) = PSO.w*Vl(pp,:) ...
            + PSO.c1*(r1k.*(pbest_l(pp,:)-Xl(pp,:))) ...
            + PSO.c2*(r2k.*(gbest_l -Xl(pp,:)));
        Vl(pp,:) = clamp_vec(Vl(pp,:), -vmax_ll, vmax_ll);
        Xl(pp,:) = Xl(pp,:) + Vl(pp,:);
        Xl(pp,:) = clamp_vec(Xl(pp,:), log(cfg.lambdaBF_min), log(cfg.lambdaBF_max));

        % ---- log-power update ----
        Vp(pp,:) = PSO.w*Vp(pp,:) ...
            + PSO.c1*(r1k.*(pbest_p(pp,:)-Xp(pp,:))) ...
            + PSO.c2*(r2k.*(gbest_p -Xp(pp,:)));
        Vp(pp,:) = clamp_vec(Vp(pp,:), -vmax_lp, vmax_lp);
        Xp(pp,:) = Xp(pp,:) + Vp(pp,:);
        % no hard bounds; softmax will normalize
    end
end

best_lambda = exp(clamp_vec(gbest_l, log(cfg.lambdaBF_min), log(cfg.lambdaBF_max))).';
best_p_rel  = softmax_vec(gbest_p);
[best_sr, best_W, best_sinr, Heff] = eval_KPBF_new(cfg, gbest_x, gbest_b, best_lambda, best_p_rel);

out.best_x = gbest_x;
out.best_betaPA = gbest_b(:);
out.best_lambda = best_lambda;
out.best_p_rel  = best_p_rel(:);
out.best_sr = best_sr;
out.best_W = best_W;
out.best_sinr = best_sinr;
out.best_Heff = Heff;
end



function f = fitness_sumrate_KPBF(cfg, x, betaPA, lambda_vec, p_rel)
[sumrate,~,~,~] = eval_KPBF_new(cfg, x, betaPA, lambda_vec, p_rel);
f = sumrate;
end


function [sumrate, W, sinr, Heff] = eval_KPBF_new(cfg, x, betaPA, lambda_vec, p_rel)
% H: N x K, G: N x M -> Heff = H^H G: K x M
H = build_H_free(cfg, x);
G = build_G_multiPA_sequential(cfg, x, betaPA);
Heff = (H') * G;  % K x M

W = kkt_parameterized_precoder_withP(Heff, cfg.Pmax, cfg.sigma2, lambda_vec, p_rel);

K = cfg.K;
sinr = zeros(K,1);
for k=1:K
    hk = Heff(k,:); % 1xM
    sig = abs(hk * W(:,k))^2;
    interf = 0;
    for j=1:K
        if j==k, continue; end
        interf = interf + abs(hk * W(:,j))^2;
    end
    sinr(k) = sig / (interf + cfg.sigma2);
end
sumrate = sum(log2(1 + sinr));
end


function W = kkt_parameterized_precoder_withP(Heff, Pmax, sigma2, lambda_vec, p_rel)
K = size(Heff,1);
M = size(Heff,2);
lambda_vec = max(lambda_vec(:), 0);
if numel(lambda_vec) ~= K
    error('lambda_vec must be Kx1 (K = number of users).');
end

p_rel = max(p_rel(:), 0);
if numel(p_rel) ~= K
    error('p_rel must be Kx1 (K = number of users).');
end
if sum(p_rel) <= 0
    p_rel = ones(K,1)/K;
else
    p_rel = p_rel / sum(p_rel);
end

Lambda = diag(lambda_vec);
A = eye(M) + (1/sigma2) * (Heff') * Lambda * Heff; % MxM
Wdir = A \ (Heff');                                % MxK
Wtilde = Wdir * diag(sqrt(p_rel));                 % MxK

pow = real(trace(Wtilde*Wtilde'));
if pow <= 0
    W = zeros(M,K);
    return;
end
W = Wtilde * sqrt(Pmax / pow);
end


function p = softmax_vec(z)
% Map unconstrained real vector z to positive simplex (sum=1)
z = z(:);
z = z - max(z);
q = exp(z);
if sum(q) <= 0
    p = ones(size(q))/numel(q);
else
    p = q / sum(q);
end
end


%% ======================================================================
%% ===================== Channel H and CMT-based G =======================
%% ======================================================================

function H = build_H_free(cfg, x)
% h_{i,k} = (lambda/(4*pi)) * exp(j k0 R_{i,k}) / R_{i,k}
N = cfg.N; K = cfg.K;
H = zeros(N,K);
for i=1:N
    for k=1:K
        R = sqrt((x(i)-cfg.a(k))^2 + cfg.y(k)^2 + cfg.hPA^2);
        H(i,k) = (cfg.lambda/(4*pi)) * exp(-1j*cfg.k0*R) / R;
    end
end
end

function H = build_H_miso(cfg, tx)
% Fair compact-array MISO channel:
% common path loss from array centroid + exact spherical-wave phase

I = numel(tx);
K = cfg.K;
H = zeros(I,K);

xc = mean(tx);   % array centroid

for k = 1:K
    Rc = sqrt((xc - cfg.a(k))^2 + cfg.y(k)^2 + cfg.hPA^2);
    beta_k = (cfg.lambda/(4*pi*Rc))^2;   % common path loss for user k

    for i = 1:I
        Ri = sqrt((tx(i) - cfg.a(k))^2 + cfg.y(k)^2 + cfg.hPA^2);
        H(i,k) = sqrt(beta_k) * exp(-1j * cfg.k0 * Ri);
    end
end
end

function H = build_H_miso_fixed_aperture(cfg, x, D)
% Fixed aperture + aperture-preserving discretization
% x : antenna x-coordinates over a fixed aperture D
% D : total aperture length

I = numel(x);
K = cfg.K;
H = zeros(I,K);

if I == 1
    d = D;
else
    d = D / (I-1);
end

d_ref = cfg.lambda / 2;
scale = sqrt(d / d_ref);

for i = 1:I
    for k = 1:K
        R = sqrt((x(i)-cfg.a(k))^2 + cfg.y(k)^2 + cfg.hPA^2);
        H(i,k) = scale * (cfg.lambda/(4*pi)) * exp(-1j*cfg.k0*R) / R;
    end
end

end

function G = build_G_multiPA_sequential(cfg, x, betaPA)
N = cfg.N;
M = cfg.M;
beta = cfg.beta(:).';      % 1 x M
x = x(:);                  % I x 1
betaPA = betaPA(:);        % I x 1
Lpa = cfg.Lpa;

% --- kappa settings ---
kappa_fixed = cfg.kappa_fixed;
eta = zeros(N,M);
for n=1:N
    for m=1:M
        Delta = betaPA(n) - beta(m);
        phi = sqrt(abs(kappa_fixed)^2 + (abs(Delta)^2)/4);
        eta(n,m) = (kappa_fixed/phi) * sin(phi*Lpa) * exp(-1j*(Lpa/2)*Delta);
    end
end

% clamp |eta|
eta_abs = abs(eta);
mask = eta_abs > 0.999999;
eta(mask) = eta(mask) .* (0.999999 ./ (eta_abs(mask) + 1e-15));

% --- sequential product in ascending x order ---
[~, ord] = sort(x, 'ascend');

G = zeros(N, M);
prod_remain = ones(1, M);

for kk = 1:N
    n = ord(kk);
    for m=1:M
        G(n,m) = eta(n,m) * prod_remain(m) * exp(-1j*beta(m)*x(n));
        prod_remain(m) = prod_remain(m) * sqrt(max(0, 1 - abs(eta(n,m))^2));
    end
end
end

%% ======================================================================
%% ======================== Baselines ===================================
%% ======================================================================

function Rtdma = tdma_singlemode_PASS(cfg, x)
% single-mode PASS + TDMA:
% slot-1 (1/2 time): use mode-1 only, user-1 only
% slot-2 (1/2 time): use mode-1 only, user-2 only
% Build full H and G
% Heff = build_H_free_same_law(cfg, x);
% G = build_G_multiPA_CMT(cfg, x);
% Heff = (H') * G; % K x M

R1 = sqrt((x(1)-cfg.a(1))^2 + cfg.y(1)^2 + cfg.hPA^2);
R2 = sqrt((x(2)-cfg.a(2))^2 + cfg.y(2)^2 + cfg.hPA^2);
% Slot 1: user1 with mode1 only
h1 = (cfg.lambda/(4*pi)) * exp(1j*cfg.k0*R1) / R1;
snr1 = cfg.Pmax * abs(h1)^2 / cfg.sigma2;

% Slot 2: user2 with mode2 only
h2 = (cfg.lambda/(4*pi)) * exp(1j*cfg.k0*R2) / R2;
snr2 = cfg.Pmax * abs(h2)^2 / cfg.sigma2;

Rtdma = 0.5*log2(1+snr1) + 0.5*log2(1+snr2);
end


function Rsum = sumrate_from_HW(H, W, sigma2)
K = size(H,2);
Rsum = 0;
for k=1:K
    hk = H(:,k);
    sig = abs(hk' * W(:,k))^2;
    int = 0;
    for j=1:K
        if j~=k
            int = int + abs(hk' * W(:,j))^2;
        end
    end
    sinr = sig / (int + sigma2);
    Rsum = Rsum + log2(1+sinr);
end
end

%% ======================================================================
%% ============================== Projection ============================
%% ======================================================================
function x = proj_and_repair(x, xmin, xmax, dmin)
x = sort(x(:).');
x = min(max(x, xmin), xmax);

% forward spacing
for i=2:numel(x)
    if x(i) < x(i-1) + dmin
        x(i) = x(i-1) + dmin;
    end
end

% shift back if overflow
if x(end) > xmax
    shift = x(end) - xmax;
    x = x - shift;
    x = max(x, xmin);
    x = sort(x);
    for i=2:numel(x)
        if x(i) < x(i-1) + dmin
            x(i) = x(i-1) + dmin;
        end
    end
end

x = min(max(x, xmin), xmax);
end

function v = clip_vec(v, vmax)
v = max(min(v, vmax), -vmax);
end

function x = clamp_around_init(x, x0, radius)
x = min(max(x, x0 - radius), x0 + radius);
end


%% ======================================================================
%% ===================== Helpers for betaPA optimization =================
%% ======================================================================
function [bmin, bmax] = default_betaPA_bounds(cfg)
bmin = min(cfg.beta);
bmax = max(cfg.beta);
end

function x = clamp_vec(x, xmin, xmax)
% Element-wise clamp
x = min(max(x, xmin), xmax);
end

function v = clip_vec_beta(v, vmax)
% Clip beta velocity element-wise to [-vmax, vmax]
v = min(max(v, -vmax), vmax);
end




%% ======================================================================
%% =============== Baseline: PSO-Parameterized BF (fixed betaPA) =========
%% ======================================================================
function out = pso_KPBF_multiPA_optx_lambda_p_fixedbeta(cfg, PSO, x0, lambda0, p0)
% Baseline: Case-1 objective, but betaPA is fixed to cfg.betaPA_fixed (I x 1).
% PSO jointly optimizes:
%   - PA positions x (continuous)
%   - Parameterized-BF variables lambda_k (positive, k=1..K)
%   - Power-loading variables p_k (relative, positive; normalized inside)

I = cfg.N; K = cfg.K;

assert(isfield(cfg,'betaPA_fixed') && ~isempty(cfg.betaPA_fixed), ...
    'cfg.betaPA_fixed must be provided for fixed-betaPA baseline.');

% velocity caps
vmax_x   = PSO.vmax_x_factor * abs(cfg.xmax - cfg.xmin);
vmax_ll  = PSO.vmax_loglambda;
vmax_lp  = PSO.vmax_logp;

% lambda bounds
if ~isfield(cfg,'lambdaBF_min') || isempty(cfg.lambdaBF_min), cfg.lambdaBF_min = 1e-4; end
if ~isfield(cfg,'lambdaBF_max') || isempty(cfg.lambdaBF_max), cfg.lambdaBF_max = 1e4;  end

Xx = zeros(PSO.P, I); Vx = zeros(PSO.P, I);
Xl = zeros(PSO.P, K); Vl = zeros(PSO.P, K);
Xp = zeros(PSO.P, K); Vp = zeros(PSO.P, K);
Xx0 = zeros(PSO.P, I);

beta_fixed_row = cfg.betaPA_fixed(:).';

for pp=1:PSO.P
    jitter_x = (2*rand(1,I)-1) * (0.5*cfg.lambda);
    Xx(pp,:) = proj_and_repair(x0 + jitter_x, cfg.xmin, cfg.xmax, cfg.dmin);
    Vx(pp,:) = (2*rand(1,I)-1) * (0.1*cfg.lambda);

    l0 = max(cfg.lambdaBF_min, min(cfg.lambdaBF_max, lambda0(:)));
    z0 = log(l0(:).');
    Xl(pp,:) = z0 + (2*rand(1,K)-1)*0.2;
    Vl(pp,:) = (2*rand(1,K)-1)*0.1;

    p0v = max(1e-8, p0(:));
    zP0 = log(p0v(:).');
    Xp(pp,:) = zP0 + (2*rand(1,K)-1)*0.2;
    Vp(pp,:) = (2*rand(1,K)-1)*0.1;
end

pbest_x = Xx; pbest_l = Xl; pbest_p = Xp;
pbest_val = -inf(PSO.P,1);
gbest_x = Xx(1,:); gbest_l = Xl(1,:); gbest_p = Xp(1,:);
gbest_val = -inf;

for t=1:PSO.T
    for pp=1:PSO.P
        lambda_vec = exp(clamp_vec(Xl(pp,:), log(cfg.lambdaBF_min), log(cfg.lambdaBF_max))).';
        p_rel = softmax_vec(Xp(pp,:));
        f = fitness_sumrate_KPBF(cfg, Xx(pp,:), beta_fixed_row, lambda_vec, p_rel);

        if f > pbest_val(pp)
            pbest_val(pp) = f;
            pbest_x(pp,:) = Xx(pp,:);
            pbest_l(pp,:) = Xl(pp,:);
            pbest_p(pp,:) = Xp(pp,:);
        end
        if f > gbest_val
            gbest_val = f;
            gbest_x = Xx(pp,:);
            gbest_l = Xl(pp,:);
            gbest_p = Xp(pp,:);
        end
    end

    for pp=1:PSO.P
        r1x = rand(1,I); r2x = rand(1,I);
        r1k = rand(1,K); r2k = rand(1,K);

        Vx(pp,:) = PSO.w*Vx(pp,:) ...
            + PSO.c1*(r1x.*(pbest_x(pp,:)-Xx(pp,:))) ...
            + PSO.c2*(r2x.*(gbest_x -Xx(pp,:)));
        Vx(pp,:) = clip_vec(Vx(pp,:), vmax_x);
        Xx(pp,:) = proj_and_repair(Xx(pp,:) + Vx(pp,:), cfg.xmin, cfg.xmax, cfg.dmin);
        Xx(pp,:) = clamp_around_init(Xx(pp,:), Xx0(pp,:), PSO.radius);
        Xx(pp,:) = proj_and_repair(Xx(pp,:), cfg.xmin, cfg.xmax, cfg.dmin);

        Vl(pp,:) = PSO.w*Vl(pp,:) ...
            + PSO.c1*(r1k.*(pbest_l(pp,:)-Xl(pp,:))) ...
            + PSO.c2*(r2k.*(gbest_l -Xl(pp,:)));
        Vl(pp,:) = clamp_vec(Vl(pp,:), -vmax_ll, vmax_ll);
        Xl(pp,:) = Xl(pp,:) + Vl(pp,:);
        Xl(pp,:) = clamp_vec(Xl(pp,:), log(cfg.lambdaBF_min), log(cfg.lambdaBF_max));

        Vp(pp,:) = PSO.w*Vp(pp,:) ...
            + PSO.c1*(r1k.*(pbest_p(pp,:)-Xp(pp,:))) ...
            + PSO.c2*(r2k.*(gbest_p -Xp(pp,:)));
        Vp(pp,:) = clamp_vec(Vp(pp,:), -vmax_lp, vmax_lp);
        Xp(pp,:) = Xp(pp,:) + Vp(pp,:);
    end
end

best_lambda = exp(clamp_vec(gbest_l, log(cfg.lambdaBF_min), log(cfg.lambdaBF_max))).';
best_p_rel  = softmax_vec(gbest_p);

out = struct();
out.best_sr    = gbest_val;
out.best_x     = gbest_x(:);
out.best_betaPA= cfg.betaPA_fixed(:);
out.best_lambda= best_lambda(:);
out.best_p_rel = best_p_rel(:);
end

