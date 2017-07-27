function [g, Jrow] = eqn_matBalSolid_Eigen(j,k,P,x,sol,ts)

% extraction of required states
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.nx:P.nx*P.nj, :);

c = sol(P.idx_c:P.ns:P.ns*P.nj, :);
Q = sol;
Q(P.idx_c:P.ns:P.ns*P.nj, :) = [];

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);

if j == 1
    a = [];
end
if j == P.nj
    d = [];
end

if j <= P.bnd_sep_neg
    Rp = P.Rp_neg;
    Ds = P.Ds_neg;
    csmax = P.csmax_neg;
elseif j >= P.bnd_pos_sep
    Rp = P.Rp_pos;
    Ds = P.Ds_pos;
    csmax = P.csmax_pos;
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

lambda = P.lambda;
N = length(lambda)-1;

tau = Ds*ts(k)/Rp^2;
dtau = Ds*P.dt/Rp^2;

delta = jn(j,k)*Rp/csmax/Ds;

sumCol1 = sum((Q((j-1)*(P.ns-1)+(P.idx_q:P.idx_qN)-1,k-1)-2*dtau*delta)./(1+dtau*lambda(1:N).^2));


sumCol2 = -2*delta*((1/10-sum(1./lambda(1:N).^2))*(1-exp(-lambda(N+1)^2*tau))+sqrt(tau/pi)*erfc(lambda(N+1)*sqrt(tau)));

c(j,k) = c(j,k-1)-3*P.dt/Rp*jn(j,k);
sumMod = csmax*sum((Q((j-1)*(P.ns-1)+(P.idx_q:P.idx_qN)-1,k-1))./(1+dtau*lambda(1:N).^2));

g = jn(j,k) + (-css(j,k) + c(j,k) + sumMod)/(-3*P.dt/Rp ...
    + Rp/Ds*sum(-2*dtau./(1+dtau*lambda(1:N).^2))...
    - 2*Rp/Ds*((1/10-sum(1./lambda(1:N).^2))*(1-exp(-lambda(N+1)^2*tau))+sqrt(tau/pi)*erfc(lambda(N+1)*sqrt(tau))));

b(P.idx_css) = 1/(-3*P.dt/Rp ...
    + Rp/Ds*sum(-2*dtau./(1+dtau*lambda(1:N).^2))...
    - 2*Rp/Ds*((1/10-sum(1./lambda(1:N).^2))*(1-exp(-lambda(N+1)^2*tau))+sqrt(tau/pi)*erfc(lambda(N+1)*sqrt(tau))));
b(P.idx_jn)  = -1;


% concatenation of the the row vectors for the band jacobian matrix
Jrow = [a b d];