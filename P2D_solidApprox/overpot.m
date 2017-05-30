function eta = overpot(x, P, T)

T = T+273.15;
phil = x(P.idx_phil:P.nx:P.nx*P.nj, :);
phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.nx:P.nx*P.nj, :); 
Eeq = zeros(size(phis));
Eeq(1:P.bnd_sep_neg,:) = Eeq_neg(css(1:P.bnd_sep_neg,:)/P.csmax_neg, P);
Eeq(P.bnd_pos_sep:P.nj,:) = Eeq_pos(css(P.bnd_pos_sep:P.nj,:)/P.csmax_pos, P);



phiRf = zeros(size(jn));
phiRf(1:P.bnd_sep_neg,:) = P.F*P.Rf_neg0*repmat(exp((P.Ebarr_neg)*(298.15-T)./(T*298.15)),P.bnd_sep_neg,1).*jn(1:P.bnd_sep_neg,:);
phiRf(P.bnd_pos_sep:P.nj,:) = P.F*P.Rf_pos0*repmat(exp((P.Ebarr_pos)*(298.15-T)./(T*298.15)),P.nj-P.bnd_pos_sep+1,1).*jn(P.bnd_pos_sep:P.nj,:);
eta = phis-phil-Eeq-phiRf;
eta = eta*1e3;