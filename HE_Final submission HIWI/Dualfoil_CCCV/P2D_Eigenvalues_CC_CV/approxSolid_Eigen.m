function sol = approxSolid_Eigen(k,P,x,sol)

c = sol(P.idx_c:P.ns:P.ns*P.nj, :);
Q = sol;
Q(P.idx_c:P.ns:P.ns*P.nj, :) = [];
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
% css = x(P.idx_css:P.ns:P.ns*P.nj, :);

for j = 1:P.nj
    if j > P.bnd_sep_neg && j < P.bnd_pos_sep
        c(j,k) = 0;
        Q((j-1)*(P.ns-1)+(P.idx_q:P.idx_qN)-1,k) = 0;
    else
        if j <= P.bnd_sep_neg
            Rp = P.Rp_neg;
            Ds = P.Ds_neg;
            csmax = P.csmax_neg;
        elseif j >= P.bnd_pos_sep
            Rp = P.Rp_pos;
            Ds = P.Ds_pos;
            csmax = P.csmax_pos;
        end

        c(j,k) = c(j,k-1) - P.dt*3*jn(j,k)/Rp;

        dtau = Ds*P.dt/Rp^2;
    
        delta = jn(j,k)*Rp/csmax/Ds;
        
        lambda = P.lambda;
        N = P.idx_qN-1;
        
        Q((j-1)*(P.ns-1)+(P.idx_q:P.idx_qN)-1,k) = (Q((j-1)*(P.ns-1)+(P.idx_q:P.idx_qN)-1,k-1)-2*dtau*delta)./(1+dtau*lambda(1:N).^2);
    end

    sol(P.idx_c+(j-1)*P.ns, k) = c(j,k);
    sol((P.idx_q:P.idx_qN)+(j-1)*P.ns, k) = Q((j-1)*(P.ns-1)+(P.idx_q:P.idx_qN)-1,k);
end

% as(j,k) = 39/4*css(j,k) - 3*q(j,k)*Rp - 35/4*c(j,k);
% bs(j,k) =-35*css(j,k) + 10*q(j,k)*Rp + 35*c(j,k);
% ds(j,k) = 105/4*css(j,k) - 7*q(j,k)*Rp - 105/4*c(j,k);