function eta = overpot(x, P, T, kk)

phil = x(P.idx_phil:P.nx:P.nx*P.nj, :);
phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.nx:P.nx*P.nj, :); 
Eeq = zeros(size(phis));
Eeq(1:P.bnd_sep_neg,:) = Eeq_neg(css(1:P.bnd_sep_neg,:)/P.csmax_neg, P);
Eeq(P.bnd_pos_sep:P.nj,:) = Eeq_pos(css(P.bnd_pos_sep:P.nj,:)/P.csmax_pos, P);


T = T+273.15;
phiRf = zeros(size(jn));
for k = 1:kk
    Rf_neg = P.Rf_neg0*exp((P.Ebarr_neg)*(298.15-T(k))/(T(k)*298.15));
    Rf_pos = P.Rf_pos0*exp((P.Ebarr_pos)*(298.15-T(k))/(T(k)*298.15));
    phiRf(1:P.bnd_sep_neg,k) = P.F*Rf_neg*jn(1:P.bnd_sep_neg,k);
    phiRf(P.bnd_pos_sep:P.nj,k) = P.F*Rf_pos*jn(P.bnd_pos_sep:P.nj,k);
end
eta = phis-phil-Eeq-phiRf;
eta = eta*1e3;