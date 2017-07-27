function [g, Jrow] = eqn_matBalSolid_Poly(j,k,P,x,sol)

% extraction of required states
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.nx:P.nx*P.nj, :);

c = sol(P.idx_c:P.ns:P.ns*P.nj, :);
q = sol(P.idx_q:P.ns:P.ns*P.nj, :);

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);

if j == 1;
    a = [];
end
if j == P.nj
    d = [];
end

if j <= P.bnd_sep_neg
    Rp = P.Rp_neg;
    Ds = P.Ds_neg;
elseif j >= P.bnd_pos_sep
    Rp = P.Rp_pos;
    Ds = P.Ds_pos;
else % separator not defined
    g = -css(j,k);
    b(P.idx_css) = 1;
    Jrow = [a b d];
    return;
end

if P.dt == 0
    g = -css(j,k) + css(j,k-1);
    
    b(P.idx_css) = 1;
    Jrow = [a b d];
    return;
end

g = (-css(j,k) + c(j,k-1)  + Rp/35*8*q(j,k-1) - 48/7*P.dt*Ds/Rp*q(j,k-1)) ...
    /(36/7*P.dt/Rp + Rp/35/Ds + P.dt*3/Rp) - jn(j,k);

b(P.idx_css) = 1/(36/7*P.dt/Rp + Rp/35/Ds + P.dt*3/Rp);
b(P.idx_jn)  = 1;


% concatenation of the the row vectors for the band jacobian matrix
Jrow = [a b d];