clear

vcut = 2.5; % stop condition for voltage
cur = 12.5; % applied current (A/m2)

%% Initialization

P = init_param;
P.dt = 0.01; % time step

P.Eeq_neg = load('Eeq_neg.mat');
P.Eeq_pos = load('Eeq_pos.mat');

P.dEeqdT_neg = load('dEeqdT_neg.mat');
P.dEeqdT_pos = load('dEeqdT_pos.mat');


P_est=P;

[x, sol,P] = init_states(P,cur);

err = 0.2;
[x_est, sol_est,P_est] = init_states_EKF(P_est,cur,err);

x1 = x_est;
sol1 = sol_est;

k = 1;
ts(k) = 0; % time 


v(k) = terminalVoltage(k,x,P,cur); % terminal voltage
v_est(k) = terminalVoltage(k,x_est,P_est,cur);
ym(k) = v(k)+randn(1)/100;

qq(k) = 0;
T(k) = P.T-273.15;
utz = zeros(P.nj,k);
utz(1:P.bnd_sep_neg,k) = P.cs0_neg/P.csmax_neg;
utz(P.bnd_pos_sep:P.nj,k) = P.cs0_pos/P.csmax_pos;

qq_est(k) = 0;
T_est(k) = P_est.T-273.15;
utz_est = zeros(P.nj,k);
utz_est(1:P.bnd_sep_neg,k) = (1+err)*P.cs0_neg/P.csmax_neg;
utz_est(P.bnd_pos_sep:P.nj,k) = (1+err)*P.cs0_pos/P.csmax_pos;



% Kalman parameters
P.Pk = diag(repmat([1e1; 1e2; 1e7], P.nj,1));
P.Qk = P.Pk;
P.Rk = 1e-4;

P_est.Pk = P.Pk;
P_est.Qk = P.Qk;
P_est.Rk = P.Rk;

% load('outcomsol8.mat')
% com = outcomsol8;
% ts = com.d1(:,1);
% ym = com.d2(:,1);

mat = load('measurement.mat');
ts = mat.ts(1,2:end);
ym = mat.v(1,2:end);
% I = mat.I(1,2:end);

dxn = P.L_neg/P.n_neg;
h1 = P.L_neg;
dxs = P.L_sep/P.n_sep;
h2 = P.L_neg+P.L_sep;
dxp = P.L_pos/P.n_pos;
h3 = P.L_neg+P.L_sep+P.L_pos;
h = [0:dxn:h1, h1+dxs:dxs:h2, h2+dxp:dxp:h3];
h =h*1e6;

%% run simulation
while k<length(ts)%v(k) > vcut
    k = k+1;
%     ts(k) = ts(k-1)+P.dt;
%     ts(k) = tvec(k);
    P.dt = ts(k)-ts(k-1);
    P_est.dt = P.dt;
%     cur = I(k); % uncomment for dynamic current
    
    % KALMAN PARAMETERS
    P.Qk = diag(repmat([(0.6*P.dt*cur/P.F/P.L_neg/P.epsl_neg)^2; ...
        (3*P.dt*cur/P.F/P.Rp_pos)^2; ...
        0*(45/2*P.dt*cur/P.F/P.Rp_pos^2)^2], P.nj,1)); % process noise cov
    P.Rk = 5e-3; % measurement noise cov matrix (voltage noise)

    P_est.Pk = P.Pk;
    P_est.Qk = P.Qk;
    P_est.Rk = P.Rk;
    
    x(:,k) = x(:,k-1);
    sol(:,k) = sol(:,k-1);
    
    [x, P] = runModel(k,P,x,cur,ts,sol);
    
    sol = approxSolid(k,P,x,sol);

    v(k) = terminalVoltage(k,x,P,cur);
    
    utz(:,k) = zeros(P.nj,1);
%     [utz(:,k-1:k), qq(k), P] = temperature(x(:,k-1:k), utz(:,k-1:k), 2, cur, v(k), P);
    T(k) = P.T-273.15;
    
    P_est.T = P.T;
    
    % add noise to measurement
    ym1(k) = ym(k)+randn(1)/1000;
    
    x_est(:,k) = x_est(:,k-1);
    sol_est(:,k) = sol_est(:,k-1);

    [x_est, sol_est,P_est] = runEKF(k,P_est,x_est,cur,ts,sol_est,ym1(k));

    v_est(k) = terminalVoltage(k,x_est,P_est,cur);
    
    utz_est(:,k) = zeros(P_est.nj,1);
%     [utz_est(:,k-1:k), qq_est(k), P_est] = temperature(x_est(:,k-1:k), utz_est(:,k-1:k), 2, cur, v_est(k), P_est);
    T_est(k) = P_est.T-273.15;
    
    
    % check mass conservation in solid and liquid phase
    check_mass
    
    
    % plot gradients of states
%     plot_states
    
%     % Adapt time stamp (comment for constant dt)
%     if P.dt <= 0
%         P.dt = 0.025;
%     elseif ts(k) < 1
%         P.dt = 0.025;
%     else
%         P.dt = .1;
%     end
%     P_est.dt = P.dt;
end