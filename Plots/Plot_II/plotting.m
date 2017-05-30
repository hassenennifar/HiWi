methods = {'Duhamel', 'Polynomial', 'Eigenvalues'};
colors = {'b', 'r--', 'k-.'};

close all
figure
subplot(211)
hold on
xlabel('Time /s')
ylabel('c_{ss}/c_{s,max} at x = L^{neg}')
subplot(212)
hold on
xlabel('Time /s')
ylabel('cl /mol m^{-3} at x = L^{neg}')
for i = 1:length(methods)
    filename = ['pulse_', methods{i}, '.mat'];
    load(filename)
    css = x(P.idx_css:P.nx:P.nx*P.nj, :);
    cl = x(P.idx_cl:P.nx:P.nx*P.nj, :);
    
    subplot(211)
    plot(ts, css(P.bnd_sep_neg,:)/P.csmax_neg, colors{i})
    
    subplot(212)
    plot(ts, cl(P.bnd_sep_neg,:)/P.csmax_neg, colors{i}) 
end