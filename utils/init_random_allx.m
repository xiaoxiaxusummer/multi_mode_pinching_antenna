function x0 = init_random_allx(cfg)
% Random feasible initialization for all x with min spacing dmin in [xmin, xmax]
% Generates sorted x with guaranteed spacing >= dmin.
%
% Method: random gaps (Dirichlet-like) over the available slack length.

N = cfg.N;
xmin = cfg.xmin; xmax = cfg.xmax; dmin = cfg.dmin;

L = xmax - xmin;
Lmin = (N-1)*dmin;

if L < Lmin
    error('Infeasible: xmax-xmin < (I-1)*dmin');
end

slack = L - Lmin;          % free length to distribute randomly
g = rand(N+1,1);
g = g / sum(g) * slack;    % gaps sum to slack

x0 = zeros(1,N);
pos = xmin + g(1);
x0(1) = pos;
for i = 2:N
    pos = pos + dmin + g(i);
    x0(i) = pos;
end

% Small numerical safety + clamp
x0 = min(max(x0, xmin), xmax);
x0 = proj_and_repair(x0, xmin, xmax, dmin);
end
