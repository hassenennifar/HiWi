function [x, sol, P] = init_states_EKF(P,cur,err, V0)

% initial states
cl = P.cl0*ones(P.nj,1)*(1+err);
phil = zeros(P.nj,1);
il = zeros(P.nj,1);
jn = zeros(P.nj,1);

% calculate consistent intial stochiometric coeff. from V0 and nLi
iter = 0;
c = [P.cs0_neg; P.cs0_pos];
while 1
    iter = iter+1;
    cs0_neg = c(1);
    cs0_pos = c(2);
    [Ua, dUa] = Eeq_neg(cs0_neg/P.csmax_neg, P);
    [Uc, dUc] = Eeq_pos(cs0_pos/P.csmax_pos, P);
    nLis = P.nLi_solid_matlab;
    g = [-V0 + Uc - Ua;  -P.nLi_solid_matlab + P.epss_pos*P.L_pos*cs0_pos + P.epss_neg*P.L_neg*cs0_neg];
    J = [-dUa, dUc; P.epss_neg*P.L_neg, P.epss_pos*P.L_pos];
    dc = -J\g;
    cnew = c + dc;
    
    % check consistency
    if cnew(1) < c(1)/1000
        cnew(1) = c(1)/1000;
    end
    if cnew(1) > c(1)*1000
           cnew(1) = 100*c(1);
    end
    if cnew(2) < c(2)/1000
        cnew(2) = c(2)/1000;
    end
    if cnew(2) > c(2)*1000
        cnew(2) = c(2)*1000;
    end
    
    c = cnew;
    
    if max(abs(dc)) < P.tolabs || iter > 500
        break;
    end
end

css = [c(1)*ones(P.n_neg+1,1); zeros(P.n_sep-1,1); c(2)*ones(P.n_pos+1,1)];
Ua = Eeq_neg(c(1)/P.csmax_neg, P);
Uc = Eeq_pos(c(2)/P.csmax_pos, P);
phis = [Ua*ones(P.n_neg+1,1); zeros(P.n_sep-1,1); Uc*ones(P.n_pos+1,1)];

% css = [P.cs0_neg*ones(P.n_neg+1,1)-err*P.csmax_neg; zeros(P.n_sep-1,1); P.cs0_pos*ones(P.n_pos+1,1)+err*P.csmax_pos];
% jn = zeros(P.nj,1);

% Ua = Eeq_neg(P.cs0_neg*(1-err)/P.csmax_neg, P);
% Uc = Eeq_pos(P.cs0_pos*(1+err)/P.csmax_pos, P);
% phis = [Ua*ones(P.n_neg+1,1); zeros(P.n_sep-1,1); Uc*ones(P.n_pos+1,1)];

% total Lithium in salt
m = zeros(1,P.nj);
m(1) = 0.5*P.epsl_neg*P.L_neg/P.n_neg;
m(2:P.bnd_sep_neg-1) = P.epsl_neg*P.L_neg/P.n_neg;
m(P.bnd_sep_neg) = 0.5*(P.epsl_neg*P.L_neg/P.n_neg+P.epsl_sep*P.L_sep/P.n_sep);
m(P.bnd_sep_neg+1:P.bnd_pos_sep-1) = P.epsl_sep*P.L_sep/P.n_sep;
m(P.bnd_pos_sep) = 0.5*(P.epsl_pos*P.L_pos/P.n_pos+P.epsl_sep*P.L_sep/P.n_sep);
m(P.bnd_pos_sep+1:P.nj-1) = P.epsl_pos*P.L_pos/P.n_pos;
m(P.nj) = 0.5*P.epsl_pos*P.L_pos/P.n_pos;
P.totLiold_est = m*cl;
% % linear spatial distribution of current density
% kappa = fun_kappa(cl, P);
% phil(P.nj) = 0;
% area = 3*P.epss_pos/P.Rp_pos;
% for j = P.nj:-1:P.bnd_pos_sep % positive
%     il(j) = cur*(P.nj-j)/P.n_pos;
%     jn(j) = -cur/P.F/P.L_pos/area;
%     if(j ~= P.nj)
%         phil(j) = phil(j+1)+P.L_pos/P.n_pos*il(j)/kappa(j); % ohmic drop is not accurate
%     end
%     phis(j) = Uc+phil(j)+jn(j)*P.F*P.Rf_pos;
% end
% 
% for j = P.bnd_pos_sep-1:-1:P.bnd_sep_neg+1 % separator
%     il(j) = cur;
%     jn(j) = 0;
%     phil(j) = phil(j+1)+P.L_sep/P.n_sep*il(j)/kappa(j); % ohmic drop is not accurate
%     phis(j) = 0;
% end
% 
% area = 3*P.epss_neg/P.Rp_neg;
% for j = P.bnd_sep_neg:-1:1 % negative electrode
%     il(j) = cur*(j-1)/P.n_neg;
%     jn(j) = cur/P.F/P.L_neg/area;
%     phil(j) = phil(j+1)+P.L_neg/P.n_neg*il(j)/kappa(j); % ohmic drop is not accurate
%     phis(j) = Ua*0+phil(j)+jn(j)*P.F*P.Rf_neg;
% end

% intial state vector
x = zeros(P.nx*P.nj, 1);
k = 1;
for i = 1:P.nx:P.nx*P.nj
    x(i+P.idx_cl-1, k) = cl((i-1)/P.nx+1);
    x(i+P.idx_phil-1, k) = phil((i-1)/P.nx+1);
    x(i+P.idx_css-1, k) = css((i-1)/P.nx+1);
    x(i+P.idx_il-1, k) = il((i-1)/P.nx+1);
    x(i+P.idx_jn-1, k) = jn((i-1)/P.nx+1);
    x(i+P.idx_phis-1, k) = phis((i-1)/P.nx+1);
end

sol = zeros(P.ns*P.nj, 1);
for i = 1:P.ns:P.ns*P.nj
    sol(i+P.idx_c-1, k) = css((i-1)/P.ns+1);
    sol(i+P.idx_q-1, k) = 0;
end