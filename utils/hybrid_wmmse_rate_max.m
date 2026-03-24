function [sumRate, F_RF, F_BB] = hybrid_wmmse_rate_max(H, Pmax, sigma2, opts)
% Two-stage method:
% 1) full-digital WMMSE
% 2) hybrid decomposition by AltMin
    
    if nargin < 5
        opts = struct();
    end

    [Nt, K] = size(H);
    NRF = K;

    fdOpts = struct();
    if isfield(opts, 'fdOpts'); fdOpts = opts.fdOpts; end

    hybOpts = struct();
    if isfield(opts, 'hybOpts'); hybOpts = opts.hybOpts; end

    [Vfd, rateFD, infoFD] = full_digital_wmmse(H, Pmax, sigma2, fdOpts);
    [F_RF, F_BB, info] = decompose_hybrid_altmin(Vfd, NRF, Pmax, hybOpts);
    
    Vhyb = F_RF * F_BB;
    [sumRate, rate_user] = compute_sum_rate_fd(H, Vhyb, sigma2);
end

function [sr,rate] = compute_sum_rate_fd(H, V, sigma2)
    [~, K] = size(H);
    for k = 1:K
        hk = H(:,k);
        hkV = hk' * V;
        sig = abs(hkV(k))^2;
        interf = sum(abs(hkV).^2) - sig;
        rate(k) = log2(1 + sig / max(interf + sigma2, 1e-12));
    end
    sr = sum(rate);
end

function [V, sumRate, info] = full_digital_wmmse(H, Pmax, sigma2, opts)
% H: Nt x K, column k is hk
% V: Nt x K full-digital precoder

    if nargin < 4
        opts = struct();
    end
    opts = set_default(opts, 'maxIter', 100);
    opts = set_default(opts, 'tol', 1e-5);
    opts = set_default(opts, 'verbose', false);

    [Nt, K] = size(H);

    % RZF initialization
    V = H / (H' * H + (K * sigma2 / Pmax) * eye(K));
    V = scale_fd_power(V, Pmax);

    rateHist = zeros(opts.maxIter,1);
    prevRate = -inf;

    for it = 1:opts.maxIter
        u = zeros(K,1);
        w = zeros(K,1);

        % update u, w
        for k = 1:K
            hk = H(:,k);
            hkV = hk' * V;
            T = sum(abs(hkV).^2) + sigma2;
            u(k) = hkV(k) / T;
            ek = 1 - 2*real(conj(u(k))*hkV(k)) + abs(u(k))^2 * T;
            ek = max(real(ek), 1e-12);
            w(k) = 1 / ek;
        end

        % update V
        A = zeros(Nt, Nt);
        B = zeros(Nt, K);
        for k = 1:K
            hk = H(:,k);
            A = A + w(k) * abs(u(k))^2 * (hk * hk');
            B(:,k) = w(k) * conj(u(k)) * hk;
        end
        A = (A + A') / 2;

        % bisection for lambda
        V0 = safe_solve(A + 1e-8*eye(Nt), B);
        if real(trace(V0' * V0)) <= Pmax
            V = V0;
        else
            lamL = 0;
            lamU = 1;
            for t = 1:60
                Vt = safe_solve(A + lamU*eye(Nt) + 1e-8*eye(Nt), B);
                if real(trace(Vt' * Vt)) <= Pmax
                    break;
                end
                lamU = 2 * lamU;
            end
            for t = 1:60
                lam = 0.5 * (lamL + lamU);
                Vt = safe_solve(A + lam*eye(Nt) + 1e-8*eye(Nt), B);
                if real(trace(Vt' * Vt)) > Pmax
                    lamL = lam;
                else
                    lamU = lam;
                end
            end
            V = safe_solve(A + lamU*eye(Nt) + 1e-8*eye(Nt), B);
        end

        rateHist(it) = compute_sum_rate_fd(H, V, sigma2);

        if opts.verbose
            fprintf('FD-WMMSE iter %3d, rate = %.6f\n', it, rateHist(it));
        end

        if it > 1
            if abs(rateHist(it)-prevRate)/max(1,abs(prevRate)) < opts.tol
                rateHist = rateHist(1:it);
                break;
            end
        end
        prevRate = rateHist(it);
    end

    sumRate = rateHist(end);
    info.rateHistory = rateHist;
end

function V = scale_fd_power(V, Pmax)
    p = real(trace(V' * V));
    p = max(p, 1e-12);
    V = sqrt(Pmax / p) * V;
end


function [F_RF, F_BB, info] = decompose_hybrid_altmin(Vfd, NRF, Pmax, opts)
    if nargin < 4
        opts = struct();
    end

    validateattributes(Vfd, {'double'}, {'2d', 'nonempty'});
    validateattributes(NRF, {'double'}, {'real','scalar','integer','>=',1});
    validateattributes(Pmax, {'double'}, {'real','scalar','positive'});

    [Nt, Ns] = size(Vfd);
    if NRF > Nt
        error('NRF must satisfy NRF <= Nt.');
    end

    opts = set_default(opts, 'maxIter', 200);
    opts = set_default(opts, 'tol', 1e-6);
    opts = set_default(opts, 'verbose', false);
    opts = set_default(opts, 'numRestart', 5);
    opts = set_default(opts, 'initMethod', 'mixed');
    opts = set_default(opts, 'muOrth', 0.02);
    opts = set_default(opts, 'regMin', 1e-8);
    opts = set_default(opts, 'regScale', 1e-6);
    opts = set_default(opts, 'reSpreadThresh', 0.90);
    opts = set_default(opts, 'acceptWorseStep', false);

    bestErr = inf;
    bestF_RF = [];
    bestF_BB = [];
    bestErrHist = [];
    bestPowHist = [];
    restartErrors = inf(opts.numRestart, 1);
    bestRestart = 1;

    for rr = 1:opts.numRestart
        % ===== Initialization =====
        F_RF = initialize_frf(Vfd, NRF, Nt, rr, opts);

        % Initial digital LS
        G = F_RF' * F_RF;
        G = hermitize(G);
        regG = adaptive_reg(G, opts);
        F_BB = safe_solve(G + regG * eye(NRF), F_RF' * Vfd);

        % Normalize power
        F_BB = scale_hybrid_power(F_RF, F_BB, Pmax);

        errHist = zeros(opts.maxIter, 1);
        powHist = zeros(opts.maxIter, 1);
        prevErr = inf;

        for it = 1:opts.maxIter
            % ===== Step 1: Update F_BB =====
            G1 = F_RF' * F_RF;
            G1 = hermitize(G1);
            reg1 = adaptive_reg(G1, opts);
            F_BB_new = safe_solve(G1 + reg1 * eye(NRF), F_RF' * Vfd);

            % ===== Step 2: Update F_RF =====
            G2 = F_BB_new * F_BB_new';
            G2 = hermitize(G2);
            reg2 = adaptive_reg(G2, opts);

            % LS analog update
            F_RF_ls = (Vfd * F_BB_new') / (G2 + reg2 * eye(NRF));

            % Constant-modulus projection
            F_RF_new = project_constant_modulus(F_RF_ls, Nt);

            % ===== Step 3: Anti-collapse nudging =====
            if opts.muOrth > 0
                F_RF_new = anti_collapse_step(F_RF_new, Nt, opts.muOrth);
            end

            % ===== Step 4: Re-spread highly correlated columns =====
            F_RF_new = respread_if_needed(F_RF_new, Nt, opts.reSpreadThresh);

            % ===== Step 5: Recompute F_BB after analog update =====
            G3 = F_RF_new' * F_RF_new;
            G3 = hermitize(G3);
            reg3 = adaptive_reg(G3, opts);
            F_BB_new = safe_solve(G3 + reg3 * eye(NRF), F_RF_new' * Vfd);

            % ===== Step 6: Power scaling =====
            F_BB_new = scale_hybrid_power(F_RF_new, F_BB_new, Pmax);

            % ===== Step 7: Evaluate =====
            curErr = norm(Vfd - F_RF_new * F_BB_new, 'fro')^2;
            curPow = hybrid_power(F_RF_new, F_BB_new);

            % Optional safeguard: reject harmful step
            if ~opts.acceptWorseStep
                oldErr = norm(Vfd - F_RF * F_BB, 'fro')^2;
                if curErr > oldErr * (1 + 1e-10)
                    % Increase regularization and re-solve once
                    reg1b = 10 * reg1;
                    reg2b = 10 * reg2;

                    F_BB_try = safe_solve(G1 + reg1b * eye(NRF), F_RF' * Vfd);
                    G2b = hermitize(F_BB_try * F_BB_try');
                    F_RF_ls_b = (Vfd * F_BB_try') / (G2b + reg2b * eye(NRF));
                    F_RF_try = project_constant_modulus(F_RF_ls_b, Nt);
                    F_RF_try = anti_collapse_step(F_RF_try, Nt, opts.muOrth);
                    F_RF_try = respread_if_needed(F_RF_try, Nt, opts.reSpreadThresh);

                    G3b = hermitize(F_RF_try' * F_RF_try);
                    reg3b = adaptive_reg(G3b, opts);
                    F_BB_try = safe_solve(G3b + reg3b * eye(NRF), F_RF_try' * Vfd);
                    F_BB_try = scale_hybrid_power(F_RF_try, F_BB_try, Pmax);

                    tryErr = norm(Vfd - F_RF_try * F_BB_try, 'fro')^2;
                    if tryErr <= curErr
                        F_RF_new = F_RF_try;
                        F_BB_new = F_BB_try;
                        curErr = tryErr;
                        curPow = hybrid_power(F_RF_new, F_BB_new);
                    end
                end
            end

            F_RF = F_RF_new;
            F_BB = F_BB_new;

            errHist(it) = curErr;
            powHist(it) = curPow;

            if opts.verbose
                fprintf(['Restart %2d | Iter %3d | err = %.6e | power = %.6f | ' ...
                         'cond(F_RF''F_RF)=%.3e\n'], ...
                        rr, it, curErr, curPow, cond_safe(F_RF' * F_RF));
            end

            if it > 1
                relChange = abs(prevErr - curErr) / max(1, abs(prevErr));
                if relChange < opts.tol
                    errHist = errHist(1:it);
                    powHist = powHist(1:it);
                    break;
                end
            end
            prevErr = curErr;
        end

        finalErr = errHist(end);
        restartErrors(rr) = finalErr;

        if finalErr < bestErr
            bestErr = finalErr;
            bestF_RF = F_RF;
            bestF_BB = F_BB;
            bestErrHist = errHist;
            bestPowHist = powHist;
            bestRestart = rr;
        end
    end

    F_RF = bestF_RF;
    F_BB = bestF_BB;

    info = struct();
    info.bestError = bestErr;
    info.errorHistory = bestErrHist;
    info.powerHistory = bestPowHist;
    info.restartErrors = restartErrors;
    info.bestRestart = bestRestart;
    info.options = opts;
end

% ======================================================================
function F_RF = initialize_frf(Vfd, NRF, Nt, rr, opts)
    [U, ~, ~] = svd(Vfd, 'econ');
    ncol = min(size(U,2), NRF);

    switch lower(opts.initMethod)
        case 'svd'
            U0 = zeros(Nt, NRF);
            U0(:,1:ncol) = U(:,1:ncol);
            if NRF > ncol
                tmp = randn(Nt, NRF - ncol) + 1j * randn(Nt, NRF - ncol);
                [Qextra, ~] = qr(tmp, 0);
                U0(:, ncol+1:NRF) = Qextra(:,1:(NRF-ncol));
            end
            F_RF = project_constant_modulus(U0, Nt);

        case 'random'
            rng(rr, 'twister');
            tmp = randn(Nt, NRF) + 1j * randn(Nt, NRF);
            [Q, ~] = qr(tmp, 0);
            F_RF = project_constant_modulus(Q(:,1:NRF), Nt);

        otherwise % 'mixed'
            if rr == 1
                U0 = zeros(Nt, NRF);
                U0(:,1:ncol) = U(:,1:ncol);
                if NRF > ncol
                    tmp = randn(Nt, NRF - ncol) + 1j * randn(Nt, NRF - ncol);
                    [Qextra, ~] = qr(tmp, 0);
                    U0(:, ncol+1:NRF) = Qextra(:,1:(NRF-ncol));
                end
                F_RF = project_constant_modulus(U0, Nt);
            else
                rng(rr, 'twister');
                tmp = randn(Nt, NRF) + 1j * randn(Nt, NRF);
                [Q, ~] = qr(tmp, 0);

                % Mix SVD subspace and random orthogonal basis
                beta = 0.7;
                U0 = zeros(Nt, NRF);
                U0(:,1:ncol) = U(:,1:ncol);

                X = beta * U0 + (1 - beta) * Q(:,1:NRF);
                F_RF = project_constant_modulus(X, Nt);
            end
    end

    F_RF = respread_if_needed(F_RF, Nt, 0.90);
end

% ======================================================================
function reg = adaptive_reg(M, opts)
    M = hermitize(M);
    d = size(M,1);
    tr = real(trace(M));
    base = opts.regScale * max(tr / max(d,1), 1);
    reg = max(opts.regMin, base);

    c = cond_safe(M);
    if ~isfinite(c)
        c = 1e16;
    end

    if c > 1e12
        reg = max(reg, 1e-3);
    elseif c > 1e10
        reg = max(reg, 1e-4);
    elseif c > 1e8
        reg = max(reg, 1e-5);
    end
end

% ======================================================================
function X = safe_solve(M, B)
    M = hermitize(M);
    r = rcond(M);

    if ~isfinite(r) || r < 1e-12
        delta = max(1e-8, 1e-6 * norm(M, 'fro'));
        M = M + delta * eye(size(M));
    end

    [L, p] = chol(M, 'lower');
    if p == 0
        X = L' \ (L \ B);
    else
        X = pinv(M) * B;
    end
end

% ======================================================================
function F_RF = anti_collapse_step(F_RF, Nt, muOrth)
    if muOrth <= 0
        return;
    end

    Rf = hermitize(F_RF' * F_RF);
    NRF = size(Rf,1);
    F_RF = F_RF - muOrth * F_RF * (Rf - eye(NRF));
    F_RF = project_constant_modulus(F_RF, Nt);
end

% ======================================================================
function F_RF = respread_if_needed(F_RF, Nt, thresh)
    NRF = size(F_RF, 2);
    Gram = abs(hermitize(F_RF' * F_RF));

    for i = 1:NRF
        for j = i+1:NRF
            denom = sqrt(real(Gram(i,i) * Gram(j,j)) + 1e-12);
            corr_ij = Gram(i,j) / denom;
            if corr_ij > thresh
                % Reinitialize column j with a phase-orthogonalized random vector
                tmp = randn(size(F_RF,1),1) + 1j*randn(size(F_RF,1),1);
                for t = 1:(j-1)
                    tmp = tmp - F_RF(:,t) * (F_RF(:,t)' * tmp);
                end
                if norm(tmp) < 1e-10
                    tmp = randn(size(F_RF,1),1) + 1j*randn(size(F_RF,1),1);
                end
                F_RF(:,j) = project_constant_modulus(tmp, Nt);
            end
        end
    end
end

% ======================================================================
function X = project_constant_modulus(X, Nt)
    X = exp(1j * angle(X)) / sqrt(Nt);
end

% ======================================================================
function p = hybrid_power(F_RF, F_BB)
    V = F_RF * F_BB;
    p = real(trace(V' * V));
end

% ======================================================================
function F_BB = scale_hybrid_power(F_RF, F_BB, Pmax)
    p = hybrid_power(F_RF, F_BB);
    p = max(p, 1e-12);
    F_BB = sqrt(Pmax / p) * F_BB;
end

% ======================================================================
function X = hermitize(X)
    X = (X + X') / 2;
end

% ======================================================================
function c = cond_safe(M)
    M = hermitize(M);
    s = svd(M);
    if isempty(s) || s(end) < 1e-16
        c = inf;
    else
        c = s(1) / s(end);
    end
end

% ======================================================================
function s = set_default(s, name, value)
    if ~isfield(s, name) || isempty(s.(name))
        s.(name) = value;
    end
end
