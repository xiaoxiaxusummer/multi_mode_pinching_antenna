function [sumrate, F_RF, F_BB] = hybridBF_MIMO(H, P, sigma2)
    if nargin < 4 || isempty(alpha)
        alpha = ones(size(H,2),1);
    end
    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'maxOuter'),  opts.maxOuter  = 50; end
    if ~isfield(opts,'maxInner'),  opts.maxInner  = 20; end
    if ~isfield(opts,'tolOuter'),  opts.tolOuter  = 1e-5; end
    if ~isfield(opts,'tolInner'),  opts.tolInner  = 1e-6; end
    if ~isfield(opts,'lambdaTol'), opts.lambdaTol = 1e-8; end
    if ~isfield(opts,'lambdaMax'), opts.lambdaMax = 1e6; end
    if ~isfield(opts,'padEps'),    opts.padEps    = 1e-10; end
    if ~isfield(opts,'init'),      opts.init      = 'svd'; end

    [Nt, K] = size(H);
    NRF = K; Ns = K;
    assert(length(alpha)==K, 'alpha must be Kx1.');

    % ----- init F_RF -----
    F_RF = (1/sqrt(Nt)) * exp(1j*2*pi*rand(Nt,NRF));

    % init target full-digital directions: matched filter (H)
    F0 = H; % Nt x K
    F_BB = pinv(F_RF) * F0; % K x K
    F_BB = normalize_power(F_RF, F_BB, P);

    wsr_hist = zeros(opts.maxOuter,1);
    prev_wsr = -inf;

    for it = 1:opts.maxOuter
        % Effective precoder
        F = F_RF * F_BB; % Nt x K

        % ----- WMMSE update u_k, w_k (single-antenna receivers) -----
        u = zeros(K,1);
        w = zeros(K,1);

        % Precompute effective gains G(k,j) = h_k^H f_j
        Heff = H' * F; % (K x K)

        for k = 1:K
            Ik = sum(abs(Heff(k,:)).^2) + sigma2;      % total received power + noise
            u(k) = Heff(k,k) / Ik;                     % MMSE equalizer (scalar)
            ek = 1 - 2*real(u(k)*Heff(k,k)) + (abs(u(k))^2)*Ik;
            w(k) = alpha(k) / max(ek, opts.padEps);
        end

        % ----- full-digital WMMSE precoder F_FD -----
        % A = sum_k w_k |u_k|^2 h_k h_k^H
        wu2 = w .* (abs(u).^2);          % Kx1
        A = H * diag(wu2) * H';          % Nt x Nt

        % B = [w1 u1* h1, ..., wK uK* hK]
        B = H * diag(w .* conj(u));      % Nt x K

        F_FD = solve_fd_with_power(A, B, P, opts.lambdaTol, opts.lambdaMax);

        % ----- inner: fit F_FD ≈ F_RF F_BB with phase-only F_RF -----
        for in = 1:opts.maxInner
            F_prev = F_RF * F_BB;

            % (i) update F_BB (LS)
            Gm = (F_RF' * F_RF) + opts.padEps*eye(NRF);
            F_BB = Gm \ (F_RF' * F_FD);            % (KxK)

            % (ii) update F_RF via phase projection
            C = (F_BB * F_BB') + opts.padEps*eye(Ns);
            Z = F_FD * (F_BB') / C;               % Nt x K
            F_RF = (1/sqrt(Nt)) * exp(1j*angle(Z));

            % (iii) normalize power
            F_BB = normalize_power(F_RF, F_BB, P);

            F_now = F_RF * F_BB;
            rel = norm(F_now - F_prev,'fro')^2 / max(norm(F_prev,'fro')^2, opts.padEps);
            if rel < opts.tolInner
                break;
            end
        end

        % ----- evaluate weighted sum-rate -----
        F = F_RF * F_BB;
        Heff = H' * F_RF; % KxK
        rates = zeros(K,1); SINR = zeros(K,1);
        for k=1:K
            hk = Heff(k,:); % 1xM
            sig = abs(hk * F_BB(:,k))^2;
            interf = 0;
            for j=1:K
                if j==k, continue; end
                interf = interf + abs(hk * F_BB(:,j))^2;
            end
            SINR(k) = sig / (interf + sigma2);
            rates(k) = log2(1 + SINR(k));
        end
        wsr = sum(alpha .* rates);
        sumrate = sum(rates);
        wsr_hist(it) = wsr;

        if abs(wsr - prev_wsr) < opts.tolOuter * max(1, abs(prev_wsr))
            wsr_hist = wsr_hist(1:it);
            break;
        end
        prev_wsr = wsr;
    end
    wsr = wsr_hist(end);
end

% ---------- helpers ----------
function F_BB = normalize_power(F_RF, F_BB, P)
    F = F_RF * F_BB;
    pow = norm(F,'fro')^2;
    if pow > 0
        F_BB = sqrt(P) * F_BB / sqrt(pow);
    end
end

function F = solve_fd_with_power(A, B, P, tol, lambdaMax)
    Nt = size(A,1);
    I = eye(Nt);

    % F0 = (A + 0*I) \ B;
    % F0 = lsqminnorm(A, B);   % 最小范数解，数值稳定
    F0 = pinv(A + 0*I)*B;

    p0 = norm(F0,'fro')^2;

    if p0 <= P
        if p0 > 0
            F = sqrt(P) * F0 / sqrt(p0);
        else
            F = F0;
        end
        return;
    end

    lo = 0;
    hi = lambdaMax;

    % grow hi if needed
    for t = 1:50
        Fhi = (A + hi*I) \ B;
        if norm(Fhi,'fro')^2 < P
            break;
        end
        hi = hi * 10;
        if hi > 1e12
            break;
        end
    end

    for t = 1:80
        mid = 0.5*(lo+hi);
        Fm = (A + mid*I) \ B;
        pm = norm(Fm,'fro')^2;

        if abs(pm - P) / P < tol
            F = Fm;
            return;
        end

        if pm > P
            lo = mid;
        else
            hi = mid;
        end
    end

    F = (A + hi*I) \ B;
end
