function eta = overpot(x, k, P)

phil = x(P.idx_phil:P.nx:P.nx*P.nj, k);
phis = x(P.idx_phis:P.nx:P.nx*P.nj, k);
jn = x(P.idx_jn:P.nx:P.nx*P.nj, k);
css = x(P.idx_css:P.nx:P.nx*P.nj, k);
Eeq = zeros(size(phis));
Eeq(1:P.bnd_sep_neg) = Eeq_neg(css(1:P.bnd_sep_neg)/P.csmax_neg, P);
Eeq(P.bnd_pos_sep:P.nj) = Eeq_pos(css(P.bnd_pos_sep:P.nj)/P.csmax_pos, P);



phiRf = zeros(size(jn));
phiRf(1:P.bnd_sep_neg) = P.F*P.Rf_neg*jn(1:P.bnd_sep_neg);
phiRf(P.bnd_pos_sep:P.nj) = P.F*P.Rf_pos*jn(P.bnd_pos_sep:P.nj);
eta = phis-phil-Eeq-phiRf;
% eta = eta*1e3;