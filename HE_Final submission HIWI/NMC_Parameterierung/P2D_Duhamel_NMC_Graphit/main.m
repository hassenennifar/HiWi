clear

vcut = 2.5; % stop condition for voltage
cur = 12.5; % applied current (A/m2)

%% Initialization

P = init_param;
P.dt = 0; % time step

P.Eeq_neg = load('Eeq_neg.mat');
P.Eeq_pos = load('Eeq_pos.mat');

P.dEeqdT_neg = load('dEeqdT_neg.mat');
P.dEeqdT_pos = load('dEeqdT_pos.mat');

x = init_states(P,cur);

k = 1;
ts(k) = 0; % time 

v(k) = terminalVoltage(k,x,P,cur); % terminal voltage
qq(k) = 0;
T(k) = P.T-273.15;
utz = zeros(P.nj,k);
utz(1:P.bnd_sep_neg,k) = P.cs0_neg/P.csmax_neg;
utz(P.bnd_pos_sep:P.nj,k) = P.cs0_pos/P.csmax_pos;

% mat = load('css_const.mat');
% I = mat.cur;
% ts = mat.ts;
%% run simulation
while v(k) > vcut
    k = k+1;
    ts(k) = ts(k-1)+P.dt;
%     P.dt = ts(k)-ts(k-1);
%     cur = I(k);
    
    x(:,k) = x(:,k-1);
    
    tic
    [x, P, iter(k-1)] = runModel(k,P,x,cur,ts);
    titer(k-1) = toc;
    
    v(k) = terminalVoltage(k,x,P,cur);
    
    utz(:,k) = zeros(P.nj,1);
    [utz(:,k-1:k), qq(k), P] = temperature(x(:,k-1:k), utz(:,k-1:k), 2, cur, v(k), P);
    T(k) = P.T-273.15;
    
    if P.dt <= 0
        P.dt = 0.025;
    elseif ts(k) < 1
        P.dt = 0.025;
    else
        P.dt = 5;
    end

end