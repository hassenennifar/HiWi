clear

vcut = 3; % stop condition for voltage
cur = -21.5*8; % applied current (A/m2)

%% Initialization

P = init_param;
P.dt = 0; % time step

x = init_states(P,cur);

k = 1;
ts(k) = 0; % time 

v(k) = terminalVoltage(k,x,P,cur); % terminal voltage
qq(k) = 0;
T(k) = P.T-273.15;
utz = zeros(P.nj,k);
utz(1:P.bnd_sep_neg,k) = P.cs0_neg/P.csmax_neg;
utz(P.bnd_pos_sep:P.nj,k) = P.cs0_pos/P.csmax_pos;
eta = zeros(P.nj,k);
% total lithium in salt
P.totLiold = 0.07;
P.mode = 'CC';
tsave = 1e10;
boost = false;
%% run simulation

while 1 %vcut
    k = k+1;
    ts(k) = ts(k-1)+P.dt;
    v(k) = v(k-1);
    x(:,k) = x(:,k-1);
    cur(k) = cur(k-1);
%     if ts(k)>0 && mod(ts(k),10) == 0 && isequal(P.mode, 'CC')
%         cur(k) = -abs(cur(k-1)+21.5*10);
%         x(end-1,k) = cur(k);
%     end
    
    
    tic
    [x, P, iter(k-1)] = runModel(k,P,x,cur(k),v(k),ts,T);
    titer(k-1) = toc;
    
    cur(k) = x(end-1,k);
    v(k) = x(end,k);
    css = x(P.idx_css:P.nx:P.nx*P.nj, :);
    
%     jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
%     phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);
%     if isequal(P.mode, 'CC')
%         cur(k) = -2*21.5;
%         x(end-1,k) = cur(k);
%         v(k) = x(end,k);
%     elseif isequal(P.mode, 'CV')
%         cur(k) = x(end-1,k);
%         v(k) = x(end,k);
%     end
    
    utz(:,k) = zeros(P.nj,1);
    [utz(:,k-1:k), qq(k), P] = temperature(x(:,k-1:k), utz(:,k-1:k), 2, cur(k), v(k), P);
    T(k) = P.T-273.15;
    P.Tk = T;
    eta = overpot(x, P, T, k);
    
    if P.dt <= 0
        P.dt = 0.025;
    elseif ts(k) <1
        P.dt = 0.025;
    else
        P.dt =1;
    end
    
    if abs(cur(k))<21.5*0.01 && isequal(P.mode, 'CV')
%         cur(k) = 0;
%         P.mode = 'RX';
%         tsave = ts(k);
%         P.dt = 10;
        break;
%     elseif css(P.bnd_sep_neg,k)/P.csmax_neg<0.8
%         P.mode = 'CC';k
%     elseif eta(P.bnd_sep_neg, k) > 0
%         il = x(P.idx_il:P.nx:P.nx*P.nj, k); 
%         phis = x(P.idx_phis:P.nx:P.nx*P.nj, k);
%         jn = x(P.idx_jn:P.nx:P.nx*P.nj, k);
%         phil = x(P.idx_phil:P.nx:P.nx*P.nj, k);
%         css = x(P.idx_css:P.nx:P.nx*P.nj, k);
%         Eeq = Eeq_neg(css(P.bnd_sep_neg)/P.csmax_neg, P);
%         
%         hx = P.L_pos/P.n_pos;
%         il_pos = hx*trapz(il(P.bnd_pos_sep:P.nj))/P.sigma_pos;
%         
%         hx = P.L_neg/P.n_neg;
%         il_neg = hx*trapz(il(1:P.bnd_sep_neg))/P.sigma_neg;
%         
%         phis_bnd_sep_neg = 0.001+phil(P.bnd_sep_neg)+Eeq+P.Rf_neg*P.F*jn(P.bnd_sep_neg);
%         
%         cur(k) = (phis(P.bnd_pos_sep)+il_pos+il_neg-phis_bnd_sep_neg-v(k))/(P.L_pos/P.sigma_pos+P.L_neg/P.sigma_neg+P.RG);
%         cur(k) = 0;
       
    elseif css(P.bnd_sep_neg,k)/P.csmax_neg>0.8 % v(k)>= 4.1 && isequal(P.mode,'CC')%css(P.bnd_sep_neg,k)/P.csmax_neg>0.8 %v(k)>= 4.1 && cur(k)<0 %&& isequal(P.mode, 'CC')
%         if boost == false
%             cur(k) = -0.5*21.5;
%             boost = true;
%         else
        P.mode = 'CV';
%         end
%     elseif isequal(P.mode, 'RX') && v(k) > 3 && (ts(k)-tsave)>60*30 && css(1,k)/P.csmax_neg > 0.6
%         cur(k) = 21.5;
%         P.mode = 'DC';
%     elseif isequal(P.mode, 'DC') && v(k) < 3
%         cur(k) = 0;
%         P.mode = 'RX';
%         tsave = ts(k);
%         P.dt = 10;
%     elseif isequal(P.mode, 'RX') && (ts(k)-tsave)>60*30
%         P.mode = 'CC';
%         cur(k) = -21.5*2.5;
    end
    x(end-1,k) = cur(k);
    
%     figure(10)
% %     il = x(P.idx_il:P.nx:P.nx*P.nj, k);
% %     dxn = P.L_neg/P.n_neg;
% %     h1 = P.L_neg;
% %     dxs = P.L_sep/P.n_sep;
% %     h2 = P.L_neg+P.L_sep;
% %     dxp = P.L_pos/P.n_pos;
% %     h3 = P.L_neg+P.L_sep+P.L_pos;
% %     h = [0:dxn:h1, h1+dxs:dxs:h2, h2+dxp:dxp:h3];
% %     h =h*1e6;
%     subplot(311)
%     cla
%     plot(ts(1:k), cur(1:k))
%     subplot(312)
%     cla
%     plot(ts(1:k), v(1:k))
%     subplot(313)
%     cla
%     eta(:,k) = eta(:,k-1);
%     eta(:,k) = overpot(x,k,P);
%     plot(ts(1:k), eta(1:P.bnd_sep_neg,1:k))
% %     css = x(P.idx_css:P.nx:P.nx*P.nj, :);
% %     plot(ts, css(1:P.bnd_sep_neg,:))
%     
%     pause(.001)
end