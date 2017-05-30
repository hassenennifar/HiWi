clear

mat{1} = load('Duhamel_2.5C.mat');
mat{2} = load('Duhamel_2.5C_EKF_Polynomial.mat');
mat{3} = load('Duhamel_2.5C_EKF_Eigenvalues.mat');

colors = {'b', 'r--', 'k-.'};
close all



for i = 1:3
    switch i
        case 1
            ts = mat{i}.ts;
            v = mat{i}.v;
            x = mat{i}.x;
            P = mat{i}.P;
            cur = mat{i}.cur;
            T = mat{i}.T;
            k = mat{i}.k;
        case 2
            ts = mat{i}.ts;
            v = mat{i}.v_est;
            x = mat{i}.x_est; 
            P = mat{i}.P;
            cur = mat{i}.I;
            T = mat{i}.T_est;
            k = mat{i}.k;
        case 3
            ts = mat{i}.ts;
            v = mat{i}.v_est;
            x = mat{i}.x_est;
            P = mat{i}.P;
            cur = mat{i}.I;
            T = mat{i}.T_est;
            k = mat{i}.k;
    end
    
    node = P.bnd_sep_neg;
    figure(1)
    hold on
    plot(ts, v, colors{i})
    xlabel('Time /s')
    ylabel('Voltage /V')
    
    figure(2)
    hold on
    css = x(P.idx_css:P.nx:P.nx*P.nj,:);
    plot(ts, css(node,:), colors{i})
    xlabel('Time /s')
    ylabel('c_{ss}/c_{s,max}')
    
    figure(10)
    hold on
    eta = overpot(x, P, T, k);
    plot(ts(2:end), eta(node,2:end), colors{i})
    xlabel('Time /s')
    ylabel('\eta /mV')
    
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
           