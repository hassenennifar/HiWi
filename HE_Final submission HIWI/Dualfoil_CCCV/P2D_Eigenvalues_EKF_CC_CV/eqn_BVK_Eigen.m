function [g, Jrow] = eqn_BVK_Eigen(j,k,P,x,sol,ts)

% extraction of required states
cl = x(P.idx_cl:P.nx:P.nx*P.nj, :);
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.nx:P.nx*P.nj, :);
phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);
phil = x(P.idx_phil:P.nx:P.nx*P.nj, :);

c = sol(P.idx_c:P.ns:P.ns*P.nj, :);
Q = sol;
Q(P.idx_c:P.ns:P.ns*P.nj, :) = [];

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);


if j <= P.bnd_sep_neg % negative electrode
    Rf = P.Rf_neg;
    csmax = P.csmax_neg;
    rka = P.rka_neg;
    Rp = P.Rp_neg;
    Ds = P.Ds_neg;
    [Eeq, dEeq] = Eeq_neg(css(j,k)/csmax, P);
elseif j >= P.bnd_pos_sep % positive electrode
    Rf = P.Rf_pos;
    csmax = P.csmax_pos;
    rka = P.rka_pos;
    Rp = P.Rp_pos;
    Ds = P.Ds_pos;
    [Eeq, dEeq] = Eeq_pos(css(j,k)/csmax, P);
else
    g = -jn(j,k);
    b(P.idx_jn) = 1;
    Jrow = [a b d];
    return;
end

eta = phis(j,k) - phil(j,k) - Eeq - P.F*Rf*jn(j,k);

% exchange current density and its derivatives wrt css and cl
i0 = P.F*rka*sqrt(csmax-css(j,k))*sqrt(css(j,k))*sqrt(cl(j,k));
di0_dcss = i0/2*(1/css(j,k)-1/(csmax-css(j,k)));
di0_dcl  = i0/2/cl(j,k);

RTF2 = 2*P.R*P.T/P.F;

coef = 80*0;

if(P.T > P.Tam)
    coef = 1;
end

denominator = 1+coef/cl(j,k)*exp(-eta/RTF2);

% derivative of jn wrt css
if P.dt == 0
    if j == 1
        a = [];
    end
    if j == P.nj
        d = [];
    end
    
    g = jn(j,k) - i0/P.F*(exp(eta/RTF2) - exp(-eta/RTF2))/denominator; g=-g;
    
    b(P.idx_cl) = - 1/P.F*di0_dcl*(exp(eta/RTF2) - exp(-eta/RTF2))/denominator ...
        - i0/P.F*(exp(eta/RTF2) - exp(-eta/RTF2))/denominator^2 ...
        *(coef/cl(j,k)^2*exp(-eta/RTF2));

    b(P.idx_phil) = i0/P.F/RTF2*(exp(eta/RTF2) + exp(-eta/RTF2))/denominator ...
        - i0/P.F*(exp(eta/RTF2) - exp(-eta/RTF2))/denominator^2 ...
        *(-coef/cl(j,k)/RTF2*exp(-eta/RTF2));

    b(P.idx_css) = -(-1/P.F*di0_dcss*(exp(eta/RTF2) - exp(-eta/RTF2)) ...
        + i0/P.F/RTF2*dEeq*(exp(eta/RTF2) + exp(-eta/RTF2)))/denominator;

    b(P.idx_jn) = 1 + i0/P.F/RTF2*P.F*Rf*(exp(eta/RTF2) + exp(-eta/RTF2))/denominator^2;

    b(P.idx_phis) = -i0/P.F/RTF2*(exp(eta/RTF2) + exp(-eta/RTF2))/denominator ...
        + i0/P.F*(exp(eta/RTF2) - exp(-eta/RTF2))/denominator^2 ...
        *(-coef/cl(j,k)/RTF2*exp(-eta/RTF2));
    
    % concatenation of the the row vectors for the band jacobian matrix
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

% zero equation
g = -jn(j,k) + (i0/P.F*(exp(eta/RTF2) - exp(-eta/RTF2)))/denominator ...
    + jn(j,k) + (-css(j,k) + c(j,k) + sumMod)/(-3*P.dt/Rp ...
    + Rp/Ds*sum(-2*dtau./(1+dtau*lambda(1:N).^2))...
    - 2*Rp/Ds*((1/10-sum(1./lambda(1:N).^2))*(1-exp(-lambda(N+1)^2*tau))+sqrt(tau/pi)*erfc(lambda(N+1)*sqrt(tau))));
    
djn = 1/(-3*P.dt/Rp ...
    + Rp/Ds*sum(-2*dtau./(1+dtau*lambda(1:N).^2))...
    - 2*Rp/Ds*((1/10-sum(1./lambda(1:N).^2))*(1-exp(-lambda(N+1)^2*tau))+sqrt(tau/pi)*erfc(lambda(N+1)*sqrt(tau))));

if j == 1
    a = [];
end
if j == P.nj
    d = [];
end

% jacobian row for BVK equation
b(P.idx_cl) = -1/P.F*di0_dcl*(exp(eta/RTF2)-exp(-eta/RTF2))/denominator ...
    - i0/P.F*(exp(eta/RTF2)-exp(-eta/RTF2))/denominator^2 ...
        *(coef/cl(j,k)^2*exp(-eta/RTF2));

b(P.idx_phil) = i0/P.F/RTF2*(exp(eta/RTF2)+exp(-eta/RTF2))/denominator ...
        - i0/P.F*(exp(eta/RTF2) - exp(-eta/RTF2))/denominator^2 ...
        *(-coef/cl(j,k)/RTF2*exp(-eta/RTF2));

b(P.idx_css) = djn +(- 1/P.F*di0_dcss*(exp(eta/RTF2)-exp(-eta/RTF2)) ...
    + i0/P.F/RTF2*(dEeq+P.F*Rf*djn)*(exp(eta/RTF2)+exp(-eta/RTF2)))/denominator;

b(P.idx_jn) = -1 - 0*i0/P.F/RTF2*P.F*Rf*(exp(eta/RTF2)+exp(-eta/RTF2)) +1;

b(P.idx_phis) = -i0/P.F/RTF2*(exp(eta/RTF2)+exp(-eta/RTF2))/denominator ...
        + i0/P.F*(exp(eta/RTF2) - exp(-eta/RTF2))/denominator^2 ...
        *(-coef/cl(j,k)/RTF2*exp(-eta/RTF2));

% concatenation of the the row vectors for the band jacobian matrix
Jrow = [a b d];