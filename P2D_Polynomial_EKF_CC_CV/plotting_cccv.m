xm = mat.x;
curm = mat.cur;
vm = mat.v;
Pm = mat.P;
Tm = mat.T; Tm = [Tm, Tm(end)];

figure(10)
ax(1) = subplot(221);
cla
hold on
plot(ts(1:k), vm(1:k), 'b')
plot(ts(1:k), v_est(1:k), '--r')
% plot(ts(1:k), v(1:k), 'g')

ax(2) = subplot(223);
cla
plot(ts(1:k), curm(1:k), 'b')

ax(3) = subplot(222);
cla
css = x(P.idx_css:P.nx:P.nx*P.nj, 1:k);
csse= x_est(P.idx_css:P.nx:P.nx*P.nj, 1:k);
cssm = xm(Pm.idx_css:Pm.nx:P.nx*Pm.nj,1:k);

hold on
plot(ts(1:k), cssm(end,1:k)/Pm.csmax_pos, 'b')
plot(ts(1:k), csse(end,1:k)/P.csmax_pos,'--r')
% plot(ts(1:k), css(end,1:k)/P.csmax_pos, 'g')

plot(ts(1:k), cssm(1,1:k)/Pm.csmax_neg, 'b')
plot(ts(1:k), csse(1,1:k)/P.csmax_neg, '--r')
% plot(ts(1:k), css(P.bnd_sep_neg,1:k)/P.csmax_neg, 'g')

ax(4) = subplot(224);
cla
phis = x(P.idx_phis:P.nx:P.nx*P.nj, 1:k);
phism = xm(Pm.idx_phis:Pm.nx:Pm.nx*Pm.nj, 1:k);
phise = x_est(P.idx_phis:P.nx:P.nx*P.nj, 1:k);

hold on
plot(ts(1:k), phism(Pm.bnd_sep_neg,1:k), 'b')
plot(ts(1:k), phise(P.bnd_sep_neg,1:k), '--r')

plot(ts(1:k), phism(1,1:k), 'b')
plot(ts(1:k), phise(1,1:k), '--r')
% plot(ts(1:k), phis(P.bnd_sep_neg,1:k), 'g')

linkaxes(ax, 'x')
pause(0.001)