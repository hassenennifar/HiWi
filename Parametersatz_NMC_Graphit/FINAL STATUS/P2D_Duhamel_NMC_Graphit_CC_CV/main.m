clear

vcut = 2.5; % stop condition for voltage
cur = -12.5*8; % applied current (A/m2)

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
eta = zeros(P.nj,k);
% total lithium in salt
P.totLiold = 0.07;
P.mode = 'CC';
tsave = 1e10;
boost = false;
P.vmax = 3.7;
%% run simulation

while 1 %vcut
    k = k+1;
    ts(k) = ts(k-1)+P.dt;
    v(k) = v(k-1);
    x(:,k) = x(:,k-1);
    cur(k) = cur(k-1);
    
    
    tic
    [x, P, iter(k-1)] = runModel(k,P,x,cur(k),v(k),ts,T);
    titer(k-1) = toc;
    
    cur(k) = x(end-1,k);
    v(k) = x(end,k);
    css = x(P.idx_css:P.nx:P.nx*P.nj, :);
    
    utz(:,k) = zeros(P.nj,1);
    [utz(:,k-1:k), qq(k), P] = temperature(x(:,k-1:k), utz(:,k-1:k), 2, cur(k), v(k), P);
    T(k) = P.T-273.15;
    P.Tk = T;
    
    if P.dt <= 0
        P.dt = 0.025;
    elseif ts(k) <1
        P.dt = 0.025;
    else
        P.dt =1;
    end
    
    if abs(cur(k))<21.5*0.01 && isequal(P.mode, 'CV')
        break;
    elseif v(k)>= P.vmax && isequal(P.mode,'CC') % css(P.bnd_sep_neg,k)/P.csmax_neg>0.8
        P.mode = 'CV';
    elseif v(k) < vcut
        break;
    end
    x(end-1,k) = cur(k);

end