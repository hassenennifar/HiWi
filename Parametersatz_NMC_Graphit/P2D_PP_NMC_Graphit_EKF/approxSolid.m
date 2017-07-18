function sol = approxSolid(k,P,x,sol)

c = sol(P.idx_c:P.ns:P.ns*P.nj, :);
q = sol(P.idx_q:P.ns:P.ns*P.nj, :);
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.ns:P.ns*P.nj, :);

for j = 1:P.nj
    if j > P.bnd_sep_neg && j < P.bnd_pos_sep
        c(j,k) = 0;
        q(j,k) = 0;
    else
        if j <= P.bnd_sep_neg
            Rp = P.Rp_neg;
            Ds = P.Ds_neg;
        elseif j >= P.bnd_pos_sep
            Rp = P.Rp_pos;
            Ds = P.Ds_pos;
        end

        c(j,k) = c(j,k-1) - P.dt*3*jn(j,k)/Rp;
        q(j,k) = 1/(1+P.dt*30*Ds/Rp/Rp)*(q(j,k-1) - P.dt*45/2*jn(j,k)/Rp/Rp);
%         c(j,k) = css(j,k)-Rp/35/Ds*8*Ds*q(j,k)+Rp/35/Ds*jn(j,k)
%         q(j,k) = 35/8/Rp*(css(j,k)-c(j,k-1))+jn(j,k)/Ds;
    end

    

    sol(P.idx_c+(j-1)*P.ns, k) = c(j,k);
    sol(P.idx_q+(j-1)*P.ns, k) = q(j,k);
end

% as(j,k) = 39/4*css(j,k) - 3*q(j,k)*Rp - 35/4*c(j,k);
% bs(j,k) =-35*css(j,k) + 10*q(j,k)*Rp + 35*c(j,k);
% ds(j,k) = 105/4*css(j,k) - 7*q(j,k)*Rp - 105/4*c(j,k);