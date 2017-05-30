clear

mat{1} = load('Duhamel_CC_CV_2.5C.mat');
mat{2} = load('Polynomial_CC_CV_2.5C.mat');
mat{3} = load('Eigenvalues_CC_CV_2.5C.mat');

colors = {'b', 'r--', 'k-.'};
close all


figure(1)
hold on
xlabel('Time /s')
ylabel('Voltage /V')

figure(2)
subplot(211)
hold on
xlabel('Time /s')
ylabel('c_{ss}/c_{s,max}')

subplot(212)
hold on
xlabel('Time /s')
ylabel('\eta /mV')


figure(3)

figure(4)


for i = 1:3
    switch i
        case 1
            ts = mat{i}.ts;
            v = mat{i}.v;
            x = mat{i}.x;
            P = mat{i}.P;
            cur = mat{i}.cur;
        case 2
            ts = mat{i}.ts;
            v = mat{i}.v;
            x = mat{i}.x; 
            P = mat{i}.P;
            cur = mat{i}.cur;
        case 3
            ts = mat{i}.ts;
            v = mat{i}.v;
            x = mat{i}.x;
            P = mat{i}.P;
            cur = mat{i}.cur;
    end
    
    node = P.bnd_sep_neg;
    figure(1)
    plot(ts, v, colors{i})
    
    figure(2)
    subplot(211)
    css = x(P.idx_css:P.nx:P.nx*P.nj,:);
    plot(ts, css(node,:), colors{i})
    
    subplot(212)
    eta = overpot(x, P);
    plot(ts, eta(node,:), colors{i})
    
    figure(3)
    hold on
    cl = x(P.idx_cl:P.nx:P.nx*P.nj,:);
    plot(ts, cl(node,:), colors{i})
    ylabel('c_l /mol m^{-3}')
    
    
    figure(4)
    hold on
    jn = x(P.idx_jn:P.nx:P.nx*P.nj,:);
    plot(ts, jn(node,:), colors{i})
    ylabel('j_n /mol m^{-2} s^{-1}')
    
    figure(5)
    hold on
    phil = x(P.idx_phil:P.nx:P.nx*P.nj,:);
    plot(ts, phil(node,:), colors{i})
    ylabel('\Phi_l /V')
    
    figure(6)
    hold on
    plot(ts, cur, colors{i})
    ylabel('Applied current density /Ah m^{-2}')
    
    figure(7)
    hold on
    phis = x(P.idx_phis:P.nx:P.nx*P.nj,:);
    plot(ts, phis(node,:), colors{i})
    ylabel('\Phi_s /V')
end
           