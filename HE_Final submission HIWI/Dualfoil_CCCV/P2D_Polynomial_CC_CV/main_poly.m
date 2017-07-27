clear

vcut = 3; % stop condition for voltage
cur = -21.5*2.5; % applied current (A/m2)

%% Initialization

P = init_param;
P.dt = 0; % time step

[x, sol] = init_states(P,cur);

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
%% run simulation
tic
while v(k) > vcut
    k = k+1;
    ts(k) = ts(k-1)+P.dt;
    cur(k) = cur(k-1);
    v(k) = v(k-1);
    
%     cur = Icur(k); % uncomment for dynamic current
    
    x(:,k) = x(:,k-1);
    sol(:,k) = sol(:,k-1);

    tic
    [x, P, iter(k-1)] = runModel(k,P,x,cur(k),ts,sol);
    titer(k-1) = toc;
    
    sol = approxSolid(k,P,x,sol);

%     v(k) = terminalVoltage(k,x,P,cur);

    cur(k) = x(end-1);
    v(k) = x(end);
    
    utz(:,k) = zeros(P.nj,1);
    [utz(:,k-1:k), qq(k), P] = temperature(x(:,k-1:k), utz(:,k-1:k), 2, cur(k), v(k), P);
    T(k) = P.T-273.15;
    
    if P.dt <= 0
        P.dt = 0.025;
    elseif ts(k) < 1
        P.dt = 0.025;
    else
        P.dt = 1;
    end
    
    if abs(cur(k))<21.5*0.01
        break;
%     elseif css(P.bnd_sep_neg,k)/P.csmax_neg<0.8
%         P.mode = 'CC';
    elseif  v(k)>= 4.1 %css(P.bnd_sep_neg,k)/P.csmax_neg>0.8 %v(k)>= 4.1 && cur(k)<0 %&& isequal(P.mode, 'CC')
        P.mode = 'CV';
    end
    
    figure(10)
%     il = x(P.idx_il:P.nx:P.nx*P.nj, k);
%     dxn = P.L_neg/P.n_neg;
%     h1 = P.L_neg;
%     dxs = P.L_sep/P.n_sep;
%     h2 = P.L_neg+P.L_sep;
%     dxp = P.L_pos/P.n_pos;
%     h3 = P.L_neg+P.L_sep+P.L_pos;
%     h = [0:dxn:h1, h1+dxs:dxs:h2, h2+dxp:dxp:h3];
%     h =h*1e6;
    subplot(311)
    cla
    plot(ts(1:k), cur(1:k))
    subplot(312)
    cla
    plot(ts(1:k), v(1:k))
    subplot(313)
    cla
    eta(:,k) = eta(:,k-1);
    eta(:,k) = overpot(x,k,P);
    plot(ts(1:k), eta(1:P.bnd_sep_neg,1:k))
%     css = x(P.idx_css:P.nx:P.nx*P.nj, :);
%     plot(ts, css(1:P.bnd_sep_neg,:))
    
    pause(.001)
end