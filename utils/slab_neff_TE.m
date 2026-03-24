%% ====================== Helper: symmetric slab TE neff solver ======================
function neff = slab_neff_TE(ncore, nclad, t, fc, m)
c = 3e8;
k0 = 2*pi*fc/c;

if t <= 0 || ncore <= nclad
    neff = nan; return;
end

V = k0*(t/2)*sqrt(ncore^2 - nclad^2);

if V <= (m*pi/2)
    neff = nan; return;
end

uL = m*pi/2 + 1e-6;
uU = min((m+1)*pi/2 - 1e-6, V - 1e-6);
if uU <= uL
    neff = nan; return;
end

if mod(m,2)==0
    F = @(u) u.*tan(u) - sqrt(max(V^2 - u.^2, 0));
else
    F = @(u) -u./tan(u) - sqrt(max(V^2 - u.^2, 0));
end

uu = linspace(uL, uU, 400);
ff = F(uu);
idx = find(sign(ff(1:end-1)).*sign(ff(2:end)) <= 0, 1, 'first');
if isempty(idx)
    neff = nan; return;
end

a = uu(idx); b = uu(idx+1);

try
    u0 = fzero(F, [a b]);
catch
    neff = nan; return;
end

bparam = 1 - (u0/V)^2;
neff2 = nclad^2 + bparam*(ncore^2 - nclad^2);
neff = sqrt(max(neff2, nan));
end