function [g, Jrow] = eqn_outVoltage(k,P,x)

Jrow = zeros(1,P.nx*P.nj+2);
% if isequal(P.mode, 'CC')
    phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);

    g = -x(end,k) + (phis(end,k)-phis(1,k)- P.RG*x(end-1,k));

    Jrow(P.idx_phis) = 1;
    Jrow(P.idx_phis*P.nj) = -1;
    Jrow(P.nj*P.nx+1) = P.RG;
    Jrow(P.nj*P.nx+2) = 1;
% elseif isequal(P.mode, 'CV')
%     g = 0;
%     Jrow(P.nj*P.nx+2) = 1;
% end