function [x_est, sol_est, P, iter] = runEKF(k,P,x_est,cur,ts,sol_est,ym)


%% Prediction step
[x_est, P] = runModel_EKF(k,P,x_est,cur,ts,sol_est);

sol_est = approxSolid_Eigen(k,P,x_est,sol_est);

y = terminalVoltage(k,x_est,P,cur);


Pk = P.Pk;

%% Correction step
% Compute state tansition matrix Fk and output matix Hk
[Fk, Hk] = state_space_P2D_Eigen(x_est, k, P);

% error covariance matrix update
Pk = Fk*Pk*Fk' + P.Qk;

% measurement innovation
y_err = ym - y;

Sk = Hk*Pk*Hk' + P.Rk;

% Kalman gain
Kk = Pk*Hk'/Sk;

% Modified Klaman gains for positive average solid concentrations
for j = P.bnd_pos_sep:P.nj-1
    Kk((j-1)*P.ne+2) = Kk(end);
end

for j = 2:P.bnd_sep_neg
    Kk((j-1)*P.ne+2) = Kk(2);
end

%% Check mass conservation
c = sol_est(P.idx_c:P.ns:P.ns*P.nj, :);
cl = x_est(P.idx_cl:P.nx:P.nx*P.nj, :);

xk = reshape([cl(:,k)'; c(:,k)'], [], 1);

xk_old = xk;

xk = xk + Kk*y_err;

xk = reshape(xk, P.ne, []);
xk_old = reshape(xk_old, P.ne, []);

ms = zeros(1,P.nj);
ms(1) = 0.5*P.epss_neg*P.L_neg/P.n_neg;
ms(2:P.bnd_sep_neg-1) = P.epss_neg*P.L_neg/P.n_neg;
ms(P.bnd_sep_neg) = 0.5*(P.epss_neg*P.L_neg/P.n_neg);
ms(P.bnd_sep_neg+1:P.bnd_pos_sep-1) = 0;
ms(P.bnd_pos_sep) = 0.5*(P.epss_pos*P.L_pos/P.n_pos);
ms(P.bnd_pos_sep+1:P.nj-1) = P.epss_pos*P.L_pos/P.n_pos;
ms(P.nj) = 0.5*P.epss_pos*P.L_pos/P.n_pos;

nLi_solid_new = ms*xk(2,:)';
nLi_solid_new_pos = ms(P.bnd_pos_sep:P.nj)*xk(2,P.bnd_pos_sep:P.nj)';
nLi_solid_new_neg = ms(1:P.bnd_sep_neg)*xk(2,1:P.bnd_sep_neg)';

xk(2,1:P.bnd_sep_neg) = xk(2,1:P.bnd_sep_neg)*(P.nLi_solid_matlab-nLi_solid_new_pos)/nLi_solid_new_neg;


%% Update error covariance matrix
Pk = (eye(size(Pk)) - Kk*Hk)*Pk;

P.Pk = Pk;

%% Store new corrected states
x_estnew = x_est;
x_estnew(P.idx_cl:P.nx:P.nx*P.nj, k) = xk(1,:)';
sol_est(P.idx_c:P.ns:P.ns*P.nj, k) = xk(2,:)';
% sol_est(P.idx_q:P.ns:P.ns*P.nj, k) = xk(3,:)';


c  = sol_est(P.idx_c:P.ns:P.ns*P.nj, k);
Q = sol_est;
Q(P.idx_c:P.ns:P.ns*P.nj, :) = [];
for j = 1:P.nj
    if j > P.bnd_sep_neg && j < P.bnd_pos_sep
        x_estnew((j-1)*P.nx+P.idx_css, k) = 0;
    else
        if j <= P.bnd_sep_neg
            Rp = P.Rp_neg;
            Ds = P.Ds_neg;
            csmax = P.csmax_neg;
        else
            Rp = P.Rp_pos;
            Ds = P.Ds_pos;
            csmax = P.csmax_pos;
        end
        
        lambda = P.lambda;
        N = length(lambda)-1;

        tau = Ds*ts(k)/Rp^2;
        dtau = Ds*P.dt/Rp^2;

        delta = x_estnew((j-1)*P.nx+P.idx_jn, k)*Rp/csmax/Ds;
        
        sumCol1 = sum((Q((j-1)*(P.ns-1)+(P.idx_q:P.idx_qN)-1,k-1)-2*dtau*delta)./(1+dtau*lambda(1:N).^2));

        sumCol2 = -2*delta*((1/10-sum(1./lambda(1:N).^2))*(1-exp(-lambda(N+1)^2*tau))+sqrt(tau/pi)*erfc(lambda(N+1)*sqrt(tau)));

        x_estnew((j-1)*P.nx+P.idx_css, k) = c(j) + csmax*(sumCol1+sumCol2);
    end
end

x_estnew(:,k) = shoeHorns(x_estnew(:,k),x_est(:,k),1,P);
x_est(:,k) = x_estnew(:,k);

%% recalculate consistent solution with updated states
iter = 0;
while 1
    iter = iter+1;

    [gz, Jz] = approxBattery_EKF(k,P,x_est,cur,ts,sol_est);

    deltaz = Jz\gz;

    for j = 1:P.nj
        z(1+(j-1)*(P.nx-2):j*(P.nx-2),1) = x_est(3+(j-1)*P.nx:j*P.nx,k);
    end
    zknew = z + deltaz;
    
    xknew = x_est(:,k);
    for j = 1:P.nj
        xknew(3+(j-1)*P.nx:j*P.nx,1) = zknew(1+(j-1)*(P.nx-2):j*(P.nx-2),1);
    end
    xknew = shoeHorns(xknew,x_est(:,k),1,P);
    
    deltax = xknew-x_est(:,k);
    x_est(:,k) = xknew;
    
    i = 30;
    if P.n_neg <= i
        i = 1;
    end

    errlim = P.tolabs*ones((P.nx-2),1);
    errlim(P.idx_jn-2) = errlim(P.idx_jn-2)*1e-6; % reduce absolute tolerance for jn
    if sum((abs(deltaz((P.nx-2)+1:(P.nx-2)*2)) < max(errlim, P.tolrel*abs(z((P.nx-2)+1:(P.nx-2)*2,1)))) ...
        | (abs(deltaz((P.nx-2)*i+1:(P.nx-2)*(i+1)))< max(errlim, P.tolrel*abs(z((P.nx-2)*i+1:(P.nx-2)*(i+1),1))))) == (P.nx-2)
        break;
    elseif iter > P.lim
        error('Solver could not converge: maximum iteration number reached');
    end

end
